

mpi_test: mpi_test.c mpi_lib.h
	mpicc mpi_test.c -o mpi_test -O0

clean:
	rm ./mpi_test

run:
	mpirun -np 4 ./mpi_test


