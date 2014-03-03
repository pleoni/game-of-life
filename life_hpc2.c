// Roberto Alfieri, Marco Borelli, Roberto De Pietri, Paolo Leoni
// University of Parma - INFN
// life_hpc2.c

char version[]="2014.03.03";
int DEBUG=1;

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>    //usleep
#include <sys/time.h>  //gettimeofday
#include <string.h>    //strcpy

// #include "grid.hh" //grid class

#include <mpi.h>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_num_threads() 0
	#define omp_get_thread_num()  1
#endif

#if _OPENACC
	#include <openacc.h>
#endif

#define MYSTRLEN 80

//////////// functions ////////////

void options(int argc, char * argv[]) ;
void usage(char * argv[]);
void init_grid(int nrows, int rmin, int rmax, int ncols, double ** grid, double prob);
//void allocate_grid(int,  int , double *** grid);
void randomize_grid(int, int , double ** grid, double base_life);
double rand_double();
void do_step(int rmin, int rmax, int cmin, int cmax, double ** grid, double ** next_grid);
void do_display(int rmin, int rmax, int cmin, int cmax, double ** grid);
void save_data(char filename[]);
void save_data_as_text(char filename[]);
void swap_grids();
//void grid_copy(int rmin, int rmax, int cmin, int cmax, double ** grid, double ** next_grid);
void random_initByTime(int rank) ;
void clearscreen();
void copy_border(int rmin, int rmax, int cmin, int cmax, double ** grid);
void log_initialize();
void log_start_main_loop(double ta, double tb);
void log_finalize(double ta, double tb, double tc);

void init_GPU();

//////////// global vars ////////////

int mygpu=0; //default GPU

int nsteps=1000;       //!< Number of Steps
int ncols=80;          //!< Number of Columns
int nrows=40;          //!< Number of rRows
double base_life=0.2;  //!< Base probability
char datafile[MYSTRLEN]="";
char hostname[MYSTRLEN];
long datasize;
int ncomp=1000;            //!< Computational load
double sum=0.0;

double ** grid;
double ** next_grid;
double ** temp_grid;

double *A, *B; // array for computation

double *col_send_l, *col_send_r, *col_recv_l, *col_recv_r; //border buffers

////// omp + MPI //////

int num_threads=1;  // default omp threads
int omp_rank=0, omp_size=1; // actual omp threads

int mpi_rank=0, mpi_size=0;
int prev_rank=0, next_rank=0;
//char mpirank_name[4];
//char mpisize_name[4];

MPI_Request request[4];
MPI_Status  status[4];

/////////////////////////////////////////
///////////////// main //////////////////
/////////////////////////////////////////

int main(int argc, char ** argv) {

  double ta, tb, tc;
  struct timeval tempo ;

  options(argc, argv);         /* optarg management */

  gethostname(hostname, MYSTRLEN);  /* get hostname */
  datasize=nrows*ncols*sizeof(double) ;     /* tot number of bytes, excluding stencils */
  //if (!strcmp(datafile,""))  sprintf(datafile,"life_%ld_%d_%d.dat",datasize,mpi_rank,mpi_size);  /* file name */

  //double *ptr;
  grid      = (double **)  malloc ( sizeof(double*) * (nrows+2)  );  // init grid
  next_grid = (double **)  malloc ( sizeof(double*) * (nrows+2)  );  // init next_grid

  /*  comp  */
  posix_memalign((void*)&(A), 64, ncomp*sizeof(double)); //allocates aligned memory
  posix_memalign((void*)&(B), 64, ncomp*sizeof(double)); //allocates aligned memory

  int i;

  for (i=0; i< ncomp; i++) A[i]=rand_double();
  for (i=0; i< ncomp; i++) B[i]=rand_double();

  /*  Buffers for borders  */
  col_send_l = (double *)  malloc ( sizeof (double) * (nrows+2) ) ;
  col_send_r = (double *)  malloc ( sizeof (double) * (nrows+2) ) ;
  col_recv_l = (double *)  malloc ( sizeof (double) * (nrows+2) ) ;
  col_recv_r = (double *)  malloc ( sizeof (double) * (nrows+2) ) ;

  /*  Init MPI  */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  prev_rank = (mpi_rank-1+mpi_size) % mpi_size;
  next_rank = (mpi_rank+1) % mpi_size;

  mygpu = mpi_rank % 2;

  //sprintf(mpirank_name,"%d",mpi_rank);  /* convert integer to string */
  //sprintf(mpisize_name,"%d",mpi_size);  /* convert integer to string */

  /*  Init OMP  */
  #ifdef _OPENMP
  omp_set_num_threads(num_threads);  // default is num_threads=1, unless the user provided a different value
  #endif
  omp_rank = omp_get_thread_num();
  omp_size = omp_get_num_threads();

  /* Initialize */

  log_initialize();

  gettimeofday(&tempo,0);  ta=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TA

  init_grid( nrows, 1, nrows+1, ncols, grid, base_life);  // -- init grids --
  init_grid( nrows, 1, nrows+1, ncols, next_grid, 0);     // ----------------

  gettimeofday(&tempo,0);  tb=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TB

  log_start_main_loop(ta,tb);

  int rmin=1;
  int cmin=1;

  int rmax=nrows;
  int cmax=ncols;

  int k;

  #if _OPENACC
  	init_GPU();
  #endif


  #pragma acc data copy(A[0:ncomp],B[0:ncomp],grid[0:nrows+2][0:ncols+2]) create(col_send_l[0:nrows+2],col_send_r[0:nrows+2],col_recv_l[0:nrows+2],col_recv_r[0:nrows+2],next_grid[0:nrows+2][0:ncols+2],sum)
  for(k=1; k<nsteps; k++) {    /* MAIN LOOP */
  
    do_step(rmin,rmax,cmin,cmax, grid, next_grid);

    #pragma acc update host (col_send_r[0:nrows+2], col_send_l[0:nrows+2], grid[0:nrows+2][0:ncols+2])

    if (DEBUG==2) do_display (1, nrows, 1, ncols,  grid);

    #pragma acc wait

    swap_grids();

    copy_border(1, nrows, 1, ncols,  grid);

    #pragma acc update device (col_recv_r[0:nrows+2], col_recv_l[0:nrows+2])

  }

  // TODO: update grid on host

  //#pragma acc data copyout(grid[nrows+2][ncols+2],next_grid[nrows+2][ncols+2])

  gettimeofday(&tempo,0); tc=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TC

  log_finalize(ta,tb,tc);

  if (mpi_rank==0 && strcmp(datafile,"")) {
    save_data(datafile);
    save_data_as_text(strcat(datafile,".txt"));
  }

  MPI_Finalize();

  return 0;

}  /* main */



////////////////////////// copy_border //////////////////////////


void  copy_border(int rmin, int rmax, int cmin, int cmax, double ** grid) {

  int i;

  // copy rows (top-bottom)
  #pragma acc kernels present(grid[nrows+2][ncols+2])
  #pragma acc loop gang independent
  for (i = cmin - 1; i <= cmax + 1; ++i) {
    grid[rmin-1][i] = grid[rmax][i];
	  grid[rmax+1][i] = grid[rmin][i];
	}

  if (mpi_size==1) {

    // copy cols (left-right)
    #pragma acc kernels present(grid[nrows+2][ncols+2])
    #pragma acc loop gang independent
      for (i = rmin - 1; i <= rmax + 1; ++i) {
        grid[i][cmin-1] = grid[i][cmax];
        grid[i][cmax+1] = grid[i][cmin];
      }

  }

  else

  {

    int tag = 999;

    MPI_Sendrecv(col_send_r, nrows + 2, MPI_DOUBLE, next_rank, tag, //send
                col_recv_l, nrows + 2, MPI_DOUBLE, prev_rank, tag, //recv
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(col_send_l, nrows + 2, MPI_DOUBLE, prev_rank, tag, // send
                col_recv_r, nrows + 2, MPI_DOUBLE, next_rank, tag, // recv
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //free(col_send);
    //free(col_recv);

    MPI_Barrier(MPI_COMM_WORLD);

  }

/////////////////////////////////////////////////

 return ;

}


////////////////////////// options //////////////////////////////////////

void options(int argc, char * argv[]) {

  int i;
  while ( (i = getopt(argc, argv, "W:vc:r:s:d:ht:f:C:n:G:")) != -1) {
    switch (i) {
    case 'c':  ncols       = strtol(optarg, NULL, 10);  break;
    case 'G':  mygpu       = strtol(optarg, NULL, 10);  break;
    case 'r':  nrows       = strtol(optarg, NULL, 10);  break;
    case 's':  nsteps      = strtol(optarg, NULL, 10);  break;
    case 'n':  ncomp       = strtol(optarg, NULL, 10);  break;
    case 't':  num_threads = strtol(optarg, NULL, 10);  break;
    case 'd':  DEBUG       = strtol(optarg, NULL, 10);  break;
    case 'f':  strncpy(datafile,optarg,MYSTRLEN-1); datafile[MYSTRLEN-1]='\0'; break;
    case 'h':  usage(argv); exit(1);
    case 'v':  printf("%s version %s\n",argv[0],version); exit(1);
    case '?':  usage(argv); exit(1);
    default:   usage(argv); exit (1);
    }
  }

}

////////////////////////// usage //////////////////////////////////////

void usage(char * argv[])  {

  printf ("\n%s [-c ncols] [-r nrows] [-t num_thr] [-s nsteps] [-G ngpu] [-d debug] [-v] [-h]",argv[0]);
  printf ("\n -d <0|1|2> : <no output | debug info (default) | display interactively> ");
  printf ("\n -s <int>   : steps  (step num. default=1000)");
  printf ("\n -n <int>   : Computation load  (default=1000)");
  printf ("\n -G <int>   : GPU (default 0)");
  printf ("\n -t <int>   : Threads num ( default = OMP_NUM_THREADS if set, otherwise = cores num");
  printf ("\n -f <filename> : output data file");
  printf ("\n -v   : version ");
  printf ("\n -h   : help ");
  printf ("\n");

}


void swap_grids() {

  // swap pointers
  temp_grid=grid;
  grid=next_grid;
  next_grid=temp_grid;

}


////////////////////////// grid_copy //////////////////////////////////////
/* // UNUSED
void grid_copy(int rmin, int rmax, int cmin, int cmax, double ** grid, double ** next_grid)
{
 int i,j;
 #pragma omp parallel for private(i,j)
	#pragma acc kernels async(2) present(grid[nrows+2][ncols+2],next_grid[nrows+2][ncols+2])
	#pragma acc loop gang independent
 for (i=rmin;i<=rmax;i++)
 {
     #pragma acc loop vector independent
     for (j=cmin;j<=cmax;j++)
	{
         grid[i][j]=next_grid[i][j];
	}
 }

} */

////////////////////////////// do_step //////////////////////////////

void do_step(int rmin, int rmax, int cmin, int cmax, double ** grid, double ** next_grid)
{
  int k,l,j,i;
  double neighbors=0.0;
  #pragma acc kernels present(grid[0:nrows+2][0:ncols+2],col_recv_l[0:nrows+2], col_recv_r[0:nrows+2])
  {
    #pragma acc loop vector independent
    for (i=0; i<nrows+2; i++) grid[i][0]=col_recv_l[i] ;  //Copy recv buff to Col 0
    #pragma acc loop vector independent
    for (i=0; i<nrows+2; i++) grid[i][ncols+1]=col_recv_r[i];  //copy recv buff to Col n+1
  }
  #pragma acc kernels present(grid[nrows+2][ncols+2],next_grid[nrows+2][ncols+2],sum,A[0:ncomp],B[0:ncomp])
  {
    #pragma acc loop gang independent
    #pragma omp parallel for private(i,j,k)
    for (i=rmin; i<=rmax; i++) {  // righe
      #pragma acc loop vector independent
      for (j=cmin; j<=cmax; j++) {  // colonne
        #pragma ivdep
        #pragma vector aligned
        #pragma acc loop independent reduction(+: sum)
        for (k=0; k < ncomp; k++)  sum += A[k] + B[k]; // COMP

        // LIFE
        neighbors = grid[i+1][j+1] + grid[i+1][j] + grid[i+1][j-1] + grid[i][j+1] + grid[i][j-1] + grid[i-1][j+1]+grid[i-1][j]+grid[i-1][j-1];
        if ( ( neighbors > 3.0 ) || ( neighbors < 2.0 ) )
          next_grid[i][j] = 0.0;
        else if ( neighbors == 3.0 )
          next_grid[i][j] = 1.0;
        else
          next_grid[i][j] =  grid[i][j];
        }
      }
    }

    // prepare buffers
    #pragma acc kernels present(grid[nrows+2][ncols+2],col_send_l[0:nrows+2],col_send_r[0:nrows+2])
    {
    #pragma acc loop vector independent
    for (i=0; i<nrows+2; i++) col_send_l[i]=grid[i][1];  // Copy Col 1 to send buff
    #pragma acc loop vector independent
    for (i=0; i<nrows+2; i++) col_send_r[i]=grid[i][ncols];  //Copy Col n to send buff
    }
}

/////////////////////////// do_display ////////////////////////////////////

void do_display(int rmin, int rmax, int cmin, int cmax, double ** grid)

{
  int i,j;
  int delay=200000;       /* usec sleep in do_display */

  clearscreen();
  for(i=cmin;i<=cmax;i++) printf("-"); printf ("\n");
  for (i=rmin;i<=rmax;i++) {
    for (j=cmin;j<=cmax;j++)
      if (grid[i][j]==0) printf(" ");
      else  printf ("x");
    printf ("\n");
  }

  usleep(delay);
}


/////////////////////////// allocate_grid ////////////////////////////////////

/* // UNUSED
void allocate_grid(int nrows, int ncols, double *** grid){

    int *a,i,j;
    (*grid) =  (double **)    malloc ( sizeof(double *) * (nrows+2) );

    for (i=0; i<nrows+2;i++) {
        (*grid)[i] =  (double *)  malloc (  sizeof (double) * (ncols+2) ) ;
        for (j=0;j<ncols+2;j++) {
            (*grid)[i][j]=0;
        }
    }
} */

/////////////////////////// init_grid ////////////////////////////////////

void init_grid(int nrows, int rmin, int rmax, int ncols, double ** grid, double prob){


  int i,j;

  random_initByTime(rand()%100);

  if (rmin==1) rmin=0;
  if (rmax==nrows) rmax=nrows+1;

  if (DEBUG == 1) printf("init-grid %d-%d %d %f\n",  rmin,rmax,ncols,prob);

  for (i=rmin; i<rmax+1;i++) {

    grid[i] =  (double *)  malloc ( sizeof (double) * (ncols+2) );
    if (prob)
      for ( j=0;j<=ncols+1;j++)
        if ( rand_double() <prob ) grid[i][j]=1;
        else  grid[i][j]=0;
  //  for (j=0;j<ncols+2;j++) (*grid)[i][j]=0;

  }

}


/////////////////////////// randomize_grid ////////////////////////////////////

void randomize_grid(int nrows, int ncols, double ** grid, double prob){

  int i,j;

  random_initByTime(rand()%100) ;
  for ( i=1;i<=nrows;i++)
    for ( j=1;j<=ncols;j++)
      if (rand_double() < prob)
        grid[i][j] = 1.0;
      else
        grid[i][j] = 0.0;

}

/////////////////////////// rand_double ////////////////////////////////////

double rand_double() {

    return (double)rand()/(double)RAND_MAX;

}

/////////////////////////// random_iniyByTime ////////////////////////////////////

void random_initByTime(int rank) {

    time_t ltime;
    time(&ltime);

    //srand((unsigned) ltime + 100*rank);

    srand((unsigned) 123);

}

/////////////////////////// clearscreen ////////////////////////////////////

void clearscreen() {

if (system( "clear" )) system( "cls" );

}


//////////////////////// init_GPU //////////////////////////////////////

#if _OPENACC
void init_GPU() {

	acc_init(acc_device_nvidia);
	int myrealgpu, num_devices;
	acc_device_t my_device_type;

	//my_device_type = acc_device_cuda; //uncomment this if you are using CAPS
	my_device_type = acc_device_nvidia; //comment this if you are using CAPS

	acc_set_device_type(my_device_type) ;

	num_devices = acc_get_num_devices(my_device_type) ;
	if (DEBUG==1) fprintf(stderr,"\nNumber of devices available: %d \n",num_devices);

	acc_set_device_num(mygpu,my_device_type);
	if (DEBUG==1) fprintf(stderr,"Trying to use GPU: %d\n",mygpu);

	myrealgpu = acc_get_device_num(my_device_type);
	if (DEBUG==1) fprintf(stderr,"Actually I am using GPU: %d\n\n",myrealgpu);

	if(mygpu != myrealgpu) {
		if (DEBUG==1) fprintf(stderr,"I cannot use the requested GPU: %d\n",mygpu);
		exit(1);
	}
}
#endif

/////////////////// logging functions /////////////////////////

void log_initialize() {

  if (DEBUG==1) fprintf(stdout,"\n%s-%d MPI_INIT mpi_size:%d omp_size:%d ncols:%d nrows:%d nsteps:%d file:%s debug:%d\n", hostname, mpi_rank, mpi_size, num_threads, ncols,nrows,nsteps,datafile, DEBUG);
  if (DEBUG==1) fprintf(stdout,"\nComp load: %d\n",ncomp);
  if (DEBUG==1) fprintf(stderr,"\n%s-%d ALLOCATE MEMORY  (%ld grid + %ld new_grid = %ld bytes ) \n", hostname,mpi_rank, datasize, datasize, datasize*2);
  if (DEBUG==1) fprintf(stderr,"%s-%d  %d iterations - Start timer\n",hostname ,mpi_rank, nsteps);
  if (DEBUG==1) fprintf(stderr,"%s-%d  START_MEM_ALLOC \n", hostname,mpi_rank);

}

void log_start_main_loop(double ta, double tb) {

  if (DEBUG==1) fprintf (stderr, "%s-%d  END_MEM_ALLOC %f sec - START_DO_STEPS %d rows:%d-%d cols:%d-%d\n", hostname,mpi_rank, tb-ta, nsteps, 0,nrows+1,0,ncols+1);

}

void log_finalize(double ta, double tb, double tc) {

  if (DEBUG==1) fprintf(stderr,"%s-%d %d/%d OMP-PARALLEL STOP\n", hostname,mpi_rank,omp_rank,omp_size);

  if (DEBUG >0) fprintf(stderr,"%s-%d - Finalize  - %f sec  \n" , hostname,mpi_rank, tc-ta);

  if (DEBUG==0) {
    if (mpi_rank==0)
      fprintf(stderr,"%d %d %d %d %d %d %f %f %f # %s \n" ,  mpi_size, omp_size, ncols, nrows, nsteps, ncomp,  tb-ta, tc-tb, tc-ta, hostname );
    else
      fprintf(stderr,"#%d %d %d %d %d %d %f %f %f # %s \n" ,  mpi_size, omp_size, ncols, nrows, nsteps, ncomp,  tb-ta, tc-tb, tc-ta, hostname );
  }

}


void save_data(char filename[]) {  // todo: make this function endian-independent

  FILE * fp_outfile = fopen(filename,"wb");
  if (fp_outfile==NULL) {
    printf("Error: couldn't open outuput file %s.\n",filename);
    return;
  }

  //#pragma acc update host (grid[nrows+2][ncols+2])

  int tot_rows = nrows+2;
  int tot_cols = ncols+2;
  int errors=0;

  int i,j;
  fwrite(&tot_rows,sizeof(tot_rows),1,fp_outfile);
  fwrite(&tot_cols,sizeof(tot_cols),1,fp_outfile);
  for (i=0; i<tot_rows; i++) {
    if ( fwrite(grid[i],sizeof(double),tot_cols,fp_outfile) != tot_cols )
      errors++;
  }

  fclose(fp_outfile);

  if (errors)
    printf("There were %d errors saving binary data to file: %s.\n",errors,filename);
  else
    printf("Data saved succesfully to file: %s.\n",filename);

}

void save_data_as_text(char filename[]) {

  FILE * fp_outfile = fopen(filename,"w");
  if (fp_outfile==NULL) {
    printf("Error: couldn't open outuput file %s.\n",filename);
    return;
  }

  int tot_rows = nrows+2;
  int tot_cols = ncols+2;
  int errors=0;

  int i,j;
  fprintf(fp_outfile,"tot_rows: %d, tot_cols: %d.\n\n",tot_rows,tot_cols);
  for (i=0; i<tot_rows; i++) {
    fprintf(fp_outfile,"[%d] ",i);
    for (j=0; j<tot_cols; j++) {
      if ( fprintf(fp_outfile,"%g ",grid[i][j]) < 0 )
        errors++;
    }
    fprintf(fp_outfile,"\n");
  }

  fclose(fp_outfile);

  if (errors)
    printf("There were %d errors saving text data to file: %s.\n",errors,filename);
  else
    printf("Data saved succesfully to file: %s.\n",filename);

}
