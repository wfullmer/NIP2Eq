# Compiler and Flags
FC = gfortran
FFLAGS = -g

# Executable name
TARGET = ../nip2eq.x

# Separate modules
MODS = mod_prec.o mod_global.o mod_inout.o mod_subs.o

# Files that make up mod_subs.f90
# Note: Order must match buildit.sh exactly
SUBS_SOURCES = header.f90 ghost.f90 grid1d.f90 init.f90 update.f90 \
               bndry.f90 exact.f90 source.f90 fou.f90 mma.f90 \
               mmu.f90 edfou2.f90 getrhs.f90 rk3s1.f90 rk3s2.f90 rk3s3.f90

# Default target
all: $(TARGET)

# Linking the executable
$(TARGET): main.f90 $(MODS)
	$(FC) $(FFLAGS) -o $(TARGET) main.f90 $(MODS)

#Rules for standalone modules
mod_prec.o: mod_prec.f90
	$(FC) $(FFLAGS) -c mod_prec.f90

mod_global.o: mod_global.f90
	$(FC) $(FFLAGS) -c mod_global.f90

mod_inout.o: mod_inout.f90
	$(FC) $(FFLAGS) -c mod_inout.f90

#Special rule to generate mod_subs.f90 then compile it
mod_subs.o: mod_subs.f90
	$(FC) $(FFLAGS) -c mod_subs.f90

mod_subs.f90: $(SUBS_SOURCES)
	cat $(SUBS_SOURCES) > mod_subs.f90
	echo 'END MODULE subs ' >> mod_subs.f90

# Clean target
clean:
	rm -f *.o .mod $(TARGET) mod_subs.f90 output*.dat

.PHONY: all clean
