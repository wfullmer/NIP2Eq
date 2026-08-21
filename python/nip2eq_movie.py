#!/usr/bin/env python3
"""
Parse NIP2Eq's amovout/umovout ASCII output (output1/2/5/6.dat) and
render it as a movie.

Each frame in one of these files looks like:

    VARIABLES = "x", "a"
    ZONE I=          24 , ZONETYPE=ORDERED,
    DATAPACKING=POINT, SOLUTIONTIME=   0.0000000000000000
         0.25000     0.80000
         0.75000     0.80000
         ...

with frames repeating back-to-back until EOF. This module reads that
directly - no Tecplot required.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

_ZONE_RE = re.compile(r'I\s*=\s*(\d+)')
_TIME_RE = re.compile(r'SOLUTIONTIME\s*=\s*([-+0-9.EeDd]+)')
_VAR_RE = re.compile(r'"([^"]+)"')


def read_tecplot_series(path):
    """Read a single amovout/umovout file into a dict of numpy arrays.

    Returns {"names": [xname, yname], "time": (nframes,),
             "data": (nframes, npts, 2)}
    """
    path = Path(path)
    lines = path.read_text().splitlines()

    names = None
    times = []
    frames = []

    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if not line.startswith('VARIABLES'):
            raise ValueError(f"{path}:{i + 1}: expected VARIABLES line, got {line!r}")
        frame_names = _VAR_RE.findall(line)
        if names is None:
            names = frame_names
        elif frame_names != names:
            raise ValueError(f"{path}:{i + 1}: variable list changed mid-file")

        npts = int(_ZONE_RE.search(lines[i + 1]).group(1))
        raw_t = _TIME_RE.search(lines[i + 2]).group(1)
        times.append(float(raw_t.replace('D', 'E').replace('d', 'e')))

        start = i + 3
        rows = [[float(v) for v in lines[start + k].split()] for k in range(npts)]
        frames.append(rows)
        i = start + npts

    return {"names": names, "time": np.array(times), "data": np.array(frames)}


def load_case(directory):
    """Load whichever of output1/2/5/6.dat exist in `directory`.

    Keys: a_num (output1.dat), u_num (output2.dat),
          a_ex  (output5.dat), u_ex  (output6.dat)
    """
    directory = Path(directory)
    files = {
        "a_num": "output1.dat",
        "u_num": "output2.dat",
        "a_ex": "output5.dat",
        "u_ex": "output6.dat",
    }
    case = {}
    for key, fname in files.items():
        fpath = directory / fname
        if fpath.exists():
            case[key] = read_tecplot_series(fpath)
    if "a_num" not in case:
        raise FileNotFoundError(f"no output1.dat (void fraction) found in {directory}")
    return case


def _ylim(data, pad_frac=0.05):
    lo, hi = data.min(), data.max()
    pad = pad_frac * (hi - lo) or 1.0
    return lo - pad, hi + pad


def animate_case(case, out=None, fps=10, stride=1, show_exact=True):
    a_num, u_num = case["a_num"], case.get("u_num")
    a_ex = case.get("a_ex") if show_exact else None
    u_ex = case.get("u_ex") if show_exact else None

    # many ICs have no exact solution and aex/uex are left as zero arrays
    if a_ex is not None and np.allclose(a_ex["data"][..., 1], 0.0):
        a_ex = None
    if u_ex is not None and np.allclose(u_ex["data"][..., 1], 0.0):
        u_ex = None

    frames = np.arange(0, a_num["time"].size, stride)

    nrows = 2 if u_num is not None else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(7, 3.5 * nrows), sharex=True, squeeze=False)
    ax_a = axes[0, 0]

    line_a_num, = ax_a.plot([], [], 'o-', ms=3, label='numerical')
    line_a_ex = None
    if a_ex is not None:
        line_a_ex, = ax_a.plot([], [], 'k--', label='exact')
    ax_a.set_ylabel(a_num["names"][1])
    ax_a.set_xlim(a_num["data"][..., 0].min(), a_num["data"][..., 0].max())
    ax_a.set_ylim(*_ylim(a_num["data"][..., 1]))
    ax_a.legend(loc='upper right')

    line_u_num = line_u_ex = None
    if u_num is not None:
        ax_u = axes[1, 0]
        line_u_num, = ax_u.plot([], [], 'o-', ms=3, color='C1', label='numerical')
        if u_ex is not None:
            line_u_ex, = ax_u.plot([], [], 'k--', label='exact')
        ax_u.set_ylabel(u_num["names"][1])
        ax_u.set_xlabel(u_num["names"][0])
        ax_u.set_xlim(u_num["data"][..., 0].min(), u_num["data"][..., 0].max())
        ax_u.set_ylim(*_ylim(u_num["data"][..., 1]))
        ax_u.legend(loc='upper right')
    else:
        ax_a.set_xlabel(a_num["names"][0])

    title = fig.suptitle('')

    def update(frame_idx):
        artists = [line_a_num]
        line_a_num.set_data(a_num["data"][frame_idx, :, 0], a_num["data"][frame_idx, :, 1])
        if line_a_ex is not None:
            line_a_ex.set_data(a_ex["data"][frame_idx, :, 0], a_ex["data"][frame_idx, :, 1])
            artists.append(line_a_ex)
        if line_u_num is not None:
            line_u_num.set_data(u_num["data"][frame_idx, :, 0], u_num["data"][frame_idx, :, 1])
            artists.append(line_u_num)
            if line_u_ex is not None:
                line_u_ex.set_data(u_ex["data"][frame_idx, :, 0], u_ex["data"][frame_idx, :, 1])
                artists.append(line_u_ex)
        title.set_text(f't = {a_num["time"][frame_idx]:.4f}')
        artists.append(title)
        return artists

    anim = animation.FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)

    if out:
        out = Path(out)
        if out.suffix.lower() == '.gif':
            anim.save(out, writer=animation.PillowWriter(fps=fps))
        else:
            anim.save(out, writer=animation.FFMpegWriter(fps=fps))
        print(f"wrote {out}")
    else:
        plt.show()

    return anim


def main():
    p = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    p.add_argument("directory", help="directory containing output1/2/5/6.dat")
    p.add_argument("-o", "--out", help="output file (.mp4 or .gif); shown interactively if omitted")
    p.add_argument("--fps", type=int, default=10)
    p.add_argument("--stride", type=int, default=1, help="use every Nth frame")
    p.add_argument("--no-exact", action="store_true", help="don't overlay the exact solution")
    args = p.parse_args()

    case = load_case(args.directory)
    animate_case(case, out=args.out, fps=args.fps, stride=args.stride, show_exact=not args.no_exact)


if __name__ == "__main__":
    main()
