/*********************************/
/* mpi_getinfo.c                 */
/*********************************/

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num()  0
#endif

//#ifdef MPI // <-- non va bene: e' il nome di un namespace in mpicxx.h
#include <mpi.h> // OPEN_MPI is defined here
//#endif

using namespace std;

#include "grid.hh"

/////// GLOBAL VARIABLES ///////

int num_threads;
int mygpu;

#ifdef OPEN_MPI
int        mpi_rank, mpi_size, tot_threads;
int        hostlen;
char       hostname[MPI_MAX_PROCESSOR_NAME];
MPI_Status status;
#endif

// user options
int nrows, ncols, nsteps, mygpu, ncomp, nthreads, DEBUG=1;
char version[] = "0.2";

////// FUNCTION DECLARATIONS //////

void init_MPI(int argc, char **argv);
int  init_omp();
void check_processes();
void parse_options(int argc, char * argv[]);
void usage(char * argv[]);
void finalize_all();

void do_step(grid &Grid, grid &Next_Grid);
void copy_borders(grid &Grid);
void do_display(grid &Grid);

void grid_tests();

/////////// MAIN ///////////

int main(int argc, char **argv)
{

  // grid_tests();
  parse_options(argc,argv);
  init_MPI(argc,argv); // Initialize MPI communicato
  init_omp();  // Get thread count in current process
  check_processes(); // Show # of processes and locations

  grid Grid( 10,10,10, 1,1,1, 3 );  /* The Grid */
  grid Next_Grid(Grid);        /* Auxiliary grid */

  if (DEBUG==1 && mpi_rank==0) {
    #pragma omp single
    {
    Grid.dump();
    Next_Grid.dump();
    }
  }

#if _OPENACC
	init_GPU(); //GPU setup
#endif

//#pragma acc data copyin(Grid,neighbors) create(Next_Grid)
for(int k=1; k<nsteps; k++) {    //**** MAIN LOOP *****//

		if (DEBUG==2) do_display(Grid);

		do_step(Grid,Next_Grid);

		swap_grids_pointers(Grid,Next_Grid);  // def. in grid.cc

    		copy_borders(Grid);

  }


  finalize_all();

  return 0;
}


/////////////// FUNCTION DEFINITIONS /////////////////

void parse_options(int argc, char * argv[]) {

  int i;
  while ( (i = getopt(argc, argv, "W:vc:r:s:d:ht:f:C:n:G:")) != -1) {
    switch (i) {
      case 'r':  nrows       = strtol(optarg, NULL, 10);  break;
      case 'c':  ncols       = strtol(optarg, NULL, 10);  break;
      case 's':  nsteps      = strtol(optarg, NULL, 10);  break;
      case 'G':  mygpu       = strtol(optarg, NULL, 10);  break;
      case 'n':  ncomp       = strtol(optarg, NULL, 10);  break;
      case 't':  nthreads    = strtol(optarg, NULL, 10);  break;
      case 'd':  grid::DEBUG = DEBUG = strtol(optarg, NULL, 10);  break;
      case 'h':  usage(argv); exit(1);
      case 'v':  printf("%s version %s\n",argv[0],version); exit(1);
      case '?':  usage(argv); exit(1);
      default:   usage(argv); exit (1);
        }
    }
}


void init_MPI(int argc, char **argv) {
#ifdef OPEN_MPI

  if ( MPI_Init(&argc, &argv) != MPI_SUCCESS ) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
    exit(-1);
    return;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(hostname,&hostlen);

  return;
#endif // OPEN_MPI
}


int init_omp() {

  #pragma omp parallel
  #pragma omp single
  num_threads=omp_get_num_threads();

#ifdef _OPENMP
  printf("-- OPENMP with %d nthreads \n",  num_threads);
#endif // _OPENMP

  return 0;

}


void do_step(grid &Grid, grid &Next_Grid) {

  int k,j,i;
  double neighbors=0.0;

	#if _OPENACC
	#pragma acc kernels present(Grid,Next_Grid,neighbors)
	#pragma acc loop independent
	#endif
  for (int m=0; m<Grid.Nvars; m++) {  // variables
    //nota: esite un algoritmo molto piu' efficiente di questo per il game of life.
	#if _OPENACC
	#pragma acc loop independent
	#endif
    for (i=0; i<Grid.N[0]; i++){
	#if _OPENACC
	#pragma acc loop independent
	#endif
      for (j=0; j<Grid.N[1]; j++) {
	#if _OPENACC
	#pragma acc loop independent
	#endif
        for (k=0; k<Grid.N[2]; k++) {
          // Domanda: quali vicini guardo? Se seguo lo stesso schema del game of life 2d, devo contare ben 26 vicini;
          // altrimenti posso evitare di andare in diagonale, e diventano solo 6, ma non si riduce al game of life se pongo dim=2;
          // esiste anche una via di mezzo, in cui ne conto 18, ma nemmeno questo si riduce al game of life in 2d...
          neighbors = Grid(m,i+1,j,k) + Grid(m,i-1,j,k) + Grid(m,i,j+1,k) + Grid(m,i,j-1,k) + Grid(m,i,j,k+1) + Grid(m,i,j,k-1); // Per ora ne conto solo 6.
          if ( neighbors > 3.0 || neighbors < 2.0 )
            Next_Grid(m,i,j,k) = 0.0;
          else if ( neighbors == 3.0 )
            Next_Grid(m,i,j,k) = 1.0;
          else
            Next_Grid(m,i,j,k) = Grid(m,i,j,k);
        } // k
      } // j
    } // i
  } // m

}


void check_processes() {  // Show # of processes and locations

#ifdef OPEN_MPI
  if(mpi_rank == 0)
    {
      printf("There are [%d] MPI Process on this job !\n",mpi_size);
      printf("-- MPI process [%2d of %d] with %d threads is executing on node: %s\n", mpi_rank, mpi_size, num_threads, hostname);
      tot_threads = num_threads;
      for(int i = 1; i < mpi_size; i++)
        {
          int r_numthreads;  // remote numthreads
          int r_hostlen;  // remote hostname length
          char r_hostname[MPI_MAX_PROCESSOR_NAME];  // remote hostname
          MPI_Recv(&r_numthreads, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(&r_hostlen, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
          MPI_Recv(r_hostname, r_hostlen+1, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
          printf("-- MPI process [%2d of %d] with %d threads is executing on node: %s\n", i, mpi_size, r_numthreads, r_hostname);
          tot_threads += num_threads;
        }
      printf("The total number of threads for this run is %d !\n",tot_threads);
    }
  else
    {
      MPI_Send(&num_threads, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&hostlen, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(hostname, hostlen+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
#endif // OPEN_MPI

}


void copy_borders(grid &Grid) {

  #pragma omp single
  {
  // TODO

 #if _OPENACC
 #pragma acc update host("new borders on Grid")
 #endif

  }

};


void do_display(grid &Grid) {

  #pragma omp single
  {
  // TODO
  }

};


void finalize_all() {

#ifdef OPEN_MPI
  MPI_Finalize(); // Finalize MPI communicator
#endif

}


void usage(char * argv[])  {

  printf ("\n%s [-c ncols] [-r nrows] [-t num_thr] [-s nsteps] [-G ngpu] [-d debug] [-v] [-h]",argv[0]);
  printf ("\n -d <0|1|2> : <no output | debug info (default) | display interactively> ");
  printf ("\n -s <int>   : steps  (step num. default=1000)");
  printf ("\n -n <int>   : Computation load  (default=1000)");
  printf ("\n -G <int>   : GPU (default 0)");
  printf ("\n -t <int>   : Threads num ( default = OMP_NUM_THREADS if set, otherwise = cores num");
  printf ("\n -v   : version ");
  printf ("\n -h   : help ");
  printf ("\n");

}


void grid_tests() {

  grid GRID(3,2,1,1,1,1,1); // (Nx,Ny,Nz,Sx,Sy,Sz,Nvars)
  GRID.vars[0][GRID.idx_ijk(0,0,0)] = 1.0;
  GRID(0,0) = 10.0;
  for(int i=0; i<GRID.size(); i++) {
    cout << "index " << i << " (" << GRID.idx_x(i) << "," << GRID.idx_y(i)  << "," << GRID.idx_z(i) <<") val is: " <<  GRID.vars[0][i] << endl;
  }
  GRID.dump();

  grid G1(10);
  G1.dump();
  G1.set_size(2,2,2);
  G1.set_stencil(1);
  G1.dump();

  G1.set_size(2,2);
  G1.set_stencil(1);
  G1.dump();

  G1.set_size(2);
  G1.set_stencil(1);
  G1.dump();

}

//////////////////////// init_GPU ///////////////////////////////////////////

#if _OPENACC
void init_GPU() {

	acc_init(acc_device_nvidia);
	int myrealgpu, num_devices;
	acc_device_t my_device_type;

	//my_device_type = acc_device_cuda; //uncomment this if you are using CAPS
	my_device_type = acc_device_nvidia; //comment this if you are using CAPS

	acc_set_device_type(my_device_type) ;
	
	num_devices = acc_get_num_devices(my_device_type) ;
	fprintf(stderr,"\nNumber of devices available: %d \n",num_devices);
	
	acc_set_device_num(mygpu,my_device_type);
	fprintf(stderr,"Trying to use GPU: %d\n",mygpu);
	
	myrealgpu = acc_get_device_num(my_device_type);
	fprintf(stderr,"Actually I am using GPU: %d\n\n",myrealgpu);

	if(mygpu != myrealgpu) {
		fprintf(stderr,"I cannot use the requested GPU: %d\n",mygpu);
	exit(1);
	}
}
#endif
