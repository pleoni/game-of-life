/*********************************/
/* mpi_getinfo.c                 */
/*********************************/

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num()  0
#endif

#ifdef OPEN_MPI
#include "mpi.h"
#endif

#include <iostream>
#include <stdio.h>
using namespace std;

#include "grid.hh"

int        nthreads;

int main(int argc, char **argv)
{


  // ---------------------------------
  // Initialize MPI communicator
  // ---------------------------------

#ifdef OPEN_MPI

  int        rank=0, np=1, tot_threads;
  int        resultlen;
  int        i;
  MPI_Status status;
  char       name[MPI_MAX_PROCESSOR_NAME];

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Get_processor_name(name, &resultlen);
#endif

  // ---------------------------------
  // Show # of process and locations
  // ---------------------------------
  #pragma omp parallel
  nthreads=omp_get_num_threads();

#ifdef _OPENMP
printf("-- OPENMP with %d nthreads \n",  nthreads);
#endif

#ifdef OPEN_MPI
  if(rank == 0)
    {
      printf("There are [%d] MPI Process on this job !\n",np);
      printf("-- MPI process [%2d of %d] with %d nthreads is executing on node: %s\n", rank,np, nthreads,name);
      tot_threads = nthreads;
      for(i = 1; i < np; i++)
        {
          MPI_Recv(&nthreads, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(&resultlen, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(name, resultlen+1, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
          printf("-- MPI process [%2d of %d] with %d nthreads is executing on node: %s\n",i,np, nthreads,name);
          tot_threads += nthreads;
        }
      printf("The total number of threads for this run is %d !\n",tot_threads);
    }
  else
    {
      MPI_Send(&nthreads, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&resultlen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(name, resultlen+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
#endif

  grid GRID(3,2,1,1,1,1,1);
  //GRID.vars[0][GRID.idx_ijk(0,0,0)] = 1.0;
  //GRID(0,0) = 10.0;
  for(int i=0;i<GRID.size();i++) {
    cout << "index " << i << " (" << GRID.idx_x(i) << "," << GRID.idx_y(i)  << "," << GRID.idx_z(i) <<") val is: " <<  GRID.vars[0][i] << endl;  
  }
  GRID.dump();

//  grid G1(10);
//  G1.dump();
//  G1.set_size(2,2,2);
//  G1.set_stencil(1);
//  G1.dump();
//
//  G1.set_size(2,2);
//  G1.set_stencil(1);
//  G1.dump();
//
//  G1.set_size(2);
//  G1.set_stencil(1);
//  G1.dump();
//  

// ---------------------------------
// Finalize MPI communicator
// ---------------------------------

#ifdef OPEN_MPI
  MPI_Finalize();
#endif

  return 0;
}
