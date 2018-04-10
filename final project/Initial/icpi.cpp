#include <mpi.h>
#include <iostream>
using namespace std;

int main(int argc,char *argv[])
{
    int  namelen, numprocs, myid;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc,&argv);

	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);	

	MPI_Get_processor_name(processor_name,&namelen);
	
	MPI_Status status;
	
	
    MPI_Finalize();
    return 0;
}
