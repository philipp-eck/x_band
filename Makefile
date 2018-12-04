F90 = ifort

OPS = -heap-arrays 1 -traceback

OBJS = x_band_PROCAR.o

SRCS = x_band_PROCAR.f90

EXE = x_band_PROCAR.x


x_band_PROCAR : $(OBJS)
		 $(F90) $(OPS) -o $(EXE) $(OBJS) 



x_band_PROCAR.o : $(SRCS)
		 $(F90) $(OPS) -c $(SRCS)

clean : 
		rm -f x_band_PROCAR.o x_band_PROCAR.x

rebuild : clean x_band_PROCAR
