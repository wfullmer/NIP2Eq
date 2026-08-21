# NIP2Eq


# Movie generation (python/nip2eq_movie.py)
The solver's `amovout`/`umovout` subroutines (in `src/mod_inout.f90`) write the
void fraction `a` and velocity `u` fields to `output1.dat`/`output2.dat`
(numerical) and `output5.dat`/`output6.dat` (exact solution, if any) every
`print freq` timesteps, in a plain-text, Tecplot-style ASCII format: each
frame is a `VARIABLES`/`ZONE`/`SOLUTIONTIME` header followed by the data
points for that frame, with frames appended back-to-back.

`python/nip2eq_movie.py` reads that format directly (no Tecplot required)
and renders it as a movie with matplotlib:

```
python3 python/nip2eq_movie.py <case_dir> -o movie.mp4 [--fps N] [--stride N] [--no-exact]
```

- `<case_dir>` is a directory containing `output1.dat`/`output2.dat` (and
  optionally `output5.dat`/`output6.dat`), e.g. a benchmark or verification
  case directory after running `nip2eq.x`.
- `-o/--out` sets the output file (`.mp4` needs ffmpeg, `.gif` needs
  Pillow); omit it to show the animation interactively instead.
- `--fps` sets the frame rate (default 10).
- `--stride` uses every Nth frame, for faster/coarser previews.
- `--no-exact` skips overlaying the exact solution even if present.

The exact-solution overlay is skipped automatically when `aex`/`uex` are
identically zero (i.e. no exact solution exists for the problem, e.g. the
water-faucet case).

Requires `numpy` and `matplotlib` (and `ffmpeg` for `.mp4` output).


# old note
This was an old note that I left to myself when I wrote this code back in 2012 or so:

need to add: if(bc = periodic) then call bndry     affter the source call
  then wrap all the "update" terms into a subroutine
  need to apply bndry after every exact call (see exact subroutine, all calls need to be removed)
  put all subroutines into a seperate file that will get cat'd into a large otherwise empty module
