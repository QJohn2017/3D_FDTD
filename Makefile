DEP = complex.o matrix.o mode.o timefield3d.o adjoint3d.o
CFLAGS = -O3 -fopenmp -lm
LIBS = 

%.o:%.c
	gcc $< -c -o $*.o $(CFLAGS) $(LIBS)

3D_FDTD: $(DEP) 3d_fdtd.c 
	gcc 3d_fdtd.c -o 3D_FDTD $(DEP) $(CFLAGS) $(LIBS)


clean: 
	rm -f  `find . -name "*.o"` 3D_FDTD
