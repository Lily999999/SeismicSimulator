nProcesses = 10
nCols =3
nRows = 3
nIterations = 5
all: makeall

makeall: runv2.c
	mpicc runv2.c -o wsn.o -lm

run:
	mpirun -np $(nProcesses) -oversubscribe wsn.o $(nCols) $(nRows) $(nIterations)

clean:
	rm wsn.o