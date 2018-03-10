#include <stdio.h>
#include <mpi.h>
#include <math.h>


int main(int argc, char *argv[]){

	double* A;
	double* xPrevAns;
		double* ptrA;
		double* ptrPrev;

	double x;
	double f;
	double eps;
	double Norm, MaxNorm;
	int rank, size;



    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	A = new double [size];
	ptrA = A;
	xPrevAns = new double [size];
	ptrPrev = xPrevAns;

	MaxNorm = -1;
	eps = 0.001;
	for(int i = 0; i < size; i++){
		*ptrPrev++ = -1;
		*ptrA++ = 1;
	}
	*(A + rank) = 2 * size;
	f = size * (size + 1) / 2 + (rank + 1) * (2 * size - 1);

	do{
		ptrA = A;		
		ptrPrev = xPrevAns;
		x = f;
		for(int i = 0; i < size; i++)
			if(i != rank)
				x -= *ptrA++ * *ptrPrev++;
			else{
				ptrA++;
				ptrPrev++;
			}
			

		x /= *(A + rank);
		Norm = fabs(x - *(xPrevAns + rank));

		MPI_Allgather(&x, 1, MPI_DOUBLE, xPrevAns, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Reduce(&Norm, &MaxNorm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Bcast(&MaxNorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	}while(MaxNorm > eps);

	if(rank==0){
		ptrPrev = xPrevAns;
		for(int i = 0; i < size; i++)
			printf("x[%d] = %f\n", i, *ptrPrev++);
	}
	delete []A;
	delete []xPrevAns;


	MPI_Finalize();

	return 0;
}
