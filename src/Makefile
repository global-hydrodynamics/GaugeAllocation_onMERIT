FC	= ifort
FFLAGS = -O3 -warn all -assume byterecl -heap-arrays -nogen-interfaces

TARGETS=\
allocate_GRDC \
HydroWeb_allocate \


all: $(TARGETS)

clean:
	$(RM) -rf *.o *.s core *~ *trace *.mod *.dSYM  $(TARGETS)

.SUFFIXES : .f90
.f90:
	$(FC) ${FFLAGS} $^ -o $@
