/*!
\mainpage LIFE
\section intro Introduzione
Game of life\n

\version   12.09.14
\author Roberto Alfieri - University of Parma - INFN


/*
  Roberto Alfieri - 

  Version History:
  12.09.14   Adapted for serial version
  12.06.13   added ifdef MPI
  11.07.15   support of OMP_NUM_THREADS
  11.05.25   openMP splitted by rows, MPI by columns; parallized (openMP) version of grid_init
  11.03.20   split main and grid routines (life_common.h)
  10.10.04   Add options -f -v
  
Serial: gcc   life.c -o life 
OMP:    gcc   life.c --define OMP -fopenmp -o life_omp 
MPI:    mpicc life.c --define MPI -o life_mpi
Hybrid: mpicc life.c --define MPI --define OMP -fopenmp -o life_hybrid

*/

char version[]="12.09.14";

int DEBUG=1;

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>    //usleep
#include <sys/time.h>  //gettimeofday
#include <string.h>    //strcpy
 
#ifdef MPI
#include <mpi.h>
#endif /* MPI */

#ifdef OMP
#include <omp.h>
#endif /* OMP */
 
///////////////////////////

void options(int argc, char * argv[]) ;
void usage(char * argv[]);
void  copy_border(int rmin, int rmax, int cmin, int cmax, short ** grid);
void init_grid(int nrows, int rmin, int rmax, int ncols, short *** grid, double prob);
void allocate_grid(int,  int , short *** grid);
void save_grid(int nrows, int ncols, short *** grid, char * file, int rank, int rank_tot);
int  load_grid(int nrows, int ncols, short *** grid, char * file, int rank, int rank_tot);
void randomize_grid(int, int , short ** grid, double base_life);
double rand_double();
void  do_step(int rmin, int rmax, int cmin, int cmax, short ** grid, short ** next_grid);
void  do_display(int rmin, int rmax, int cmin, int cmax, short ** grid);
void grid_copy(int rmin, int rmax, int cmin, int cmax, short ** grid, short ** next_grid);
void random_initByTime(int rank) ;
void clearscreen();
void copy_border(int rmin, int rmax, int cmin, int cmax, short ** grid);

 int nsteps=1000;       //!< Number of Steps 
 int ncols=80;          //!< Number of Columns 
 int nrows=40;          //!< Number of rRows 
 double base_life=0.2;  //!< Base probability
 int CKPTsteps=500;     //!< CheckPoint intervals (steps) 
 int Write=0;           //!< if Write=1  save the grid to datafile after execution 
 char datafile[80]=""; 
 char hostname[80];
 long datasize;
 
// OMP
 int num_threads=1;

 //MPI
 char mpirank_name[4]; 
 char mpisize_name[4]; 

 short ** grid;              
 short ** next_grid;

/////////// MPI /////////

int mpi_rank=0, mpi_size=0;  
int prev_rank=0, next_rank=0;  

#ifdef MPI
MPI_Request request[4]; 
MPI_Status  status[4];
#endif /* MPI */

/////////// OMP //////////
int omp_rank=0, omp_size=1; 

/////////////////////////////////////////
int main(int argc, char ** argv) {

 double  ta, tb;
 struct timeval tempo ;

 options(argc, argv);         /* optarg management */ 
 gethostname(hostname, 100);  /* get hostname */
 datasize=nrows*ncols*sizeof(short) ;     /* tot number of bytes */ 
	
 short *ptr;
 grid      = (short **)  malloc ( sizeof(ptr) * (nrows+2)  );  // init grid
 next_grid = (short **)  malloc ( sizeof(ptr) * (nrows+2)  );  // init next_grid
	

#ifdef OMP
  char *nt=getenv("OMP_NUM_THREADS");
  printf ("valore OMP_NUM_THREADS=%s \n", nt);
  if ( nt != NULL )  num_threads=strtol(nt,NULL,0);
  omp_set_num_threads(num_threads);  // if is set (>0) keep this value else use  OMP_NUM_THREADS  
#endif /* OMP */


#ifdef MPI
   MPI_Init(&argc, &argv); 
   MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
   prev_rank = (mpi_rank-1+mpi_size) % mpi_size;
   next_rank = (mpi_rank+1) % mpi_size; 
   sprintf(mpirank_name,"%d",mpi_rank);  /* convert integer to string */
   sprintf(mpisize_name,"%d",mpi_size);  /* convert integer to string */ 
#endif /* MPI */

	if (!strcmp(datafile,""))
        	sprintf(datafile,"life_%ld_%d_%d.dat",datasize,mpi_rank,mpi_size);  /* file name */	

//initialize
 
    if (DEBUG==1) fprintf(stdout,"\n%s-%d MPI_INIT mpi_size:%d omp_size:%d ncols:%d nrows:%d nsteps:%d file:%s debug:%d\n", 
                                     hostname, mpi_rank, mpi_size,num_threads,  ncols,nrows,nsteps,datafile, DEBUG);
    if (DEBUG==1) fprintf(stderr, "\n%s-%d ALLOCATE MEMORY  (%ld grid + %ld new_grid = %ld bytes ) \n",
	                                                   hostname,mpi_rank, datasize, datasize, datasize*2);
												   													   
    if (DEBUG==1) fprintf(stderr, "%s-%d  %d iterations - Start timer\n",hostname ,mpi_rank, nsteps);
		
    gettimeofday(&tempo,0);  ta=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TA
    
  #pragma omp parallel  private(omp_rank)
  {
 #ifdef OMP 
    omp_rank=omp_get_thread_num();
    omp_size=omp_get_num_threads();
#endif /* OMP */	
	
    if (DEBUG==1) fprintf(stderr,"%s-%d %d/%d OMP-PARALLEL START_MEM_ALLOC \n", hostname,mpi_rank,omp_rank,omp_size); 
  
    int rmin=(nrows*omp_rank)/omp_size+1;    //first row
    int rmax=(nrows*(omp_rank+1))/omp_size;  //last row 
	int cmin=1;      // first col
    int cmax=ncols;  // last col
	
    init_grid( nrows, rmin, rmax, ncols, &grid, base_life);
    init_grid( nrows, rmin, rmax, ncols, &next_grid, 0);
    #pragma omp barrier

    gettimeofday(&tempo,0);  tb=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TB
   
   if (DEBUG==1) fprintf (stderr, "%s-%d %d/%d OMP_PARALLEL END_MEM_ALLOC %f sec - START_DO_STEPS %d rows:%d-%d cols:%d-%d\n",
                                    hostname,mpi_rank,omp_rank,omp_size, tb-ta, nsteps, rmin,rmax,cmin,cmax);
 
									
	int k;		

 for(k=1; k<nsteps; k++)
     {
 //     if (DEBUG==1) printf ("%d \n",k);
      do_step(rmin,rmax,cmin,cmax, grid, next_grid);
     #pragma omp barrier
     grid_copy(rmin,rmax,cmin,cmax, grid, next_grid);
     #pragma omp barrier
	 
//	if (DEBUG==1) fprintf(stderr,"%s-%d %d/%d OMP-PARALLEL STEP:%d\n", hostname,mpi_rank,omp_rank,omp_size,k); 
  
  #pragma omp single
      {
      copy_border(1, nrows, 1, ncols,  grid);
      if (DEBUG==2) do_display (1, nrows, 1, ncols,  grid);
	  
      // checkpointing - begin /////////////////////////////////////
       
       if (!(k%CKPTsteps)) {
         if (Write==1) {
           gettimeofday(&tempo,0);  tb=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TB
           if (DEBUG==1) fprintf(stderr, "\n%s Checkpoint writing started at %d steps after %8.2f sec from start \n", hostname, k, (tb-ta));
            save_grid( nrows, ncols,  &grid, datafile,0,1); 
           gettimeofday(&tempo,0);  tb=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TB
           if (DEBUG==1) fprintf(stderr, "\n%s Checkpoint writing finished at %d steps after %8.2f sec from start \n", hostname, k, (tb-ta));

            }
    	}
	   // checkpointing - end /////////////////////////////////////
      }
    }

	if (DEBUG==1) fprintf(stderr,"%s-%d %d/%d OMP-PARALLEL STOP\n", hostname,mpi_rank,omp_rank,omp_size); 
	
  } // //end openMP parallel
  
    gettimeofday(&tempo,0); tb=tempo.tv_sec+(tempo.tv_usec/1000000.0); // Save current time in TB
	
    if (DEBUG>0) fprintf(stderr,"%s-%d Finalize  - %f sec  \n" , hostname,mpi_rank, tb-ta);

	if (Write==1) {
	  if (DEBUG==1) fprintf(stderr, "\n%s-%d PROGRAM END, saving on file %s..  \n",hostname,mpi_rank, datafile);
      save_grid( nrows, ncols,  &grid, datafile, mpi_rank, mpi_size);
    }

#ifdef MPI
    MPI_Finalize();
#endif /* MPI */
 
  return 0;
}

////////////////////////////// border //////////////////////////////

void  copy_border(int rmin, int rmax, int cmin, int cmax, short ** grid) {

  int i;

// copy rows
  for (i = cmin - 1; i <= cmax + 1; ++i) {
    grid[rmin-1][i] = grid[rmax][i];
    grid[rmax+1][i] = grid[rmin][i];
  }

#ifndef MPI

// copy cols
  for (i = rmin - 1; i <= rmax + 1; ++i) {
    grid[rmin-1][i] = grid[rmax][i];
    grid[rmax+1][i] = grid[rmin][i];
  }

#else

////////////////////////////// MPI 
int tag = 999;

// create 2 temp columns 
short *col_send, *col_recv; 
col_send= (short *)  malloc (  sizeof (short) * (nrows+2) ) ;
col_recv= (short *)  malloc (  sizeof (short) * (nrows+2) ) ;

for (i=0; i<nrows+2; i++) col_send[i]=grid[i][ncols];  //Copy Col n to send buff

 MPI_Sendrecv(col_send, nrows + 2, MPI_SHORT, next_rank, tag, //send
              col_recv, nrows + 2, MPI_SHORT, prev_rank, tag, //recv
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//printf("MPIsendrecv send to %d - recv from %d\n", next_rank, prev_rank); 

for (i=0; i<nrows+2; i++)  grid[i][0]=col_recv[i] ;  //Copy recv buff to Col 0


for (i=0; i<nrows+2; i++) col_send[i]=grid[i][1];  // Copy Col 1 to send buff

 MPI_Sendrecv(col_send, nrows + 2, MPI_SHORT, prev_rank, tag, // send
              col_recv, nrows + 2, MPI_SHORT, next_rank, tag, // recv 
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//printf("MPIsendrecv send to %d - recv from %d\n", prev_rank, next_rank);

for (i=0; i<nrows+2; i++) grid[i][ncols+1]=col_recv[i];  //copy recv buff to Col n+1
 
free(col_send);
free(col_recv);

MPI_Barrier(MPI_COMM_WORLD);

#endif /* MPI */

/////////////////////////////////////////////////  

 return ;
}

////////////////////////// options //////////////////////////////////////

void options(int argc, char * argv[]) {

  int i;
   while ( (i = getopt(argc, argv, "W:vc:r:s:d:ht:f:C:")) != -1) {
        switch (i) {
        case 'c':  ncols       = strtol(optarg, NULL, 10);  break;
        case 'r':  nrows       = strtol(optarg, NULL, 10);  break;
        case 's':  nsteps      = strtol(optarg, NULL, 10);  break;
        case 't':  num_threads = strtol(optarg, NULL, 10);  break;
        case 'd':  DEBUG       = strtol(optarg, NULL, 10);  break;
        case 'W':  Write       = strtol(optarg, NULL, 10);  break;			
		case 'C':  CKPTsteps   = strtol(optarg, NULL, 10);  break;	
        case 'f':  strcpy(datafile,optarg)               ;  break;
        case 'h':  usage(argv); exit(1);
        case 'v':  printf("%s version %s\n",argv[0],version); exit(1);
        case '?':  usage(argv); exit(1);
        default:   usage(argv); exit (1);
        }
    }
}

////////////////////////// usage //////////////////////////////////////

void usage(char * argv[])  {

  printf ("\n%s [-c ncols] [-r nrows] [-t num_thr] [-s nsteps] [-d debug] [-f filename] [-W 1|0] [-C CKPTsteps] [-v] [-h]",argv[0]); 
  printf ("\n -d <0|1|2> : <no output | debug info (default) | display interactively> ");   
  printf ("\n -s <int>   : steps  (step num. default=1000)"); 
  printf ("\n -C <int>   : CheckPoint freq. (step num. default=500)");
  printf ("\n -t <int>   : Threads num ( default = OMP_NUM_THREADS if set, otherwise = cores num"); 
  printf ("\n -W <1|0>   : Write yes|no(default) the grid config into a file  "); 
  printf ("\n -f <string>: filename  (default: life_ <filesize>.dat) "); 
  printf ("\n -v   : version "); 
  printf ("\n -h   : help "); 
  printf ("\n"); 

}

////////////////////////// grid_copy //////////////////////////////////////

void grid_copy(int rmin, int rmax, int cmin, int cmax, short ** grid, short ** next_grid)
{
int i,j;
 for (i=rmin;i<=rmax;i++)
                for (j=cmin;j<=cmax;j++) grid[i][j]=next_grid[i][j];  
}

////////////////////////////// do_step //////////////////////////////

void do_step(int rmin, int rmax, int cmin, int cmax, short ** grid, short ** next_grid)
 { 
  int i,j,k,l;
  for (i=rmin;i<=rmax;i++) {
                for (j=cmin;j<=cmax;j++) {
                        int neighbors=0;
                        for (k=i-1;k<=i+1; k++) {
                                for (l=j-1;l<=j+1; l++) {
                                        if (!(k==i&&l==j)&&grid[k][l]>0) {
                                                neighbors++;
                                        }
                                }
                                if (neighbors>3) continue;
                        }
                        if (neighbors<2||neighbors>3) {
                                next_grid[i][j]=0;
                        } else if (grid[i][j]>0||neighbors==3) {
                                next_grid[i][j]=grid[i][j]+1;
                        }
                }
        }

}

/////////////////////////// do_display ////////////////////////////////////

void do_display(int rmin, int rmax, int cmin, int cmax, short ** grid)
{
  int i,j; 
  int delay=800000;       /* usec sleep in do_display */

 
  clearscreen();

  for(i=cmin;i<=cmax;i++) printf("-"); printf ("\n");

    for (i=rmin;i<=rmax;i++) 
         { for (j=cmin;j<=cmax;j++)  
	    if (grid[i][j]==0) printf(" "); else printf ("x"); printf ("\n"); }
		
   //for(i=cmin;i<=cmax;i++) printf("-"); printf ("\n");
   usleep(delay);

 }

/////////////////////////// allocate_grid ////////////////////////////////////

void allocate_grid(int nrows, int ncols, short *** grid){

    int *a,i,j;

    (*grid) =  (short **)    malloc ( sizeof(a) * (nrows+2) );
    for (i=0; i<nrows+2;i++) {
        (*grid)[i] =  (short *)  malloc (  sizeof (short) * (ncols+2) ) ;

        for (j=0;j<ncols+2;j++) {
            (*grid)[i][j]=0;    
        }
    }
}

/////////////////////////// init_grid ////////////////////////////////////

void init_grid(int nrows, int rmin, int rmax, int ncols, short *** grid, double prob){

    int i,j; 

    random_initByTime(rand()%100);
    if (rmin==1) rmin=0;
    if (rmax==nrows) rmax=nrows+1;  

    printf("init-grid %d-%d %d %f\n",  rmin,rmax,ncols,prob);
	
    for (i=rmin; i<rmax+1;i++) 
  	{
	(*grid)[i] =  (short *)  malloc (  sizeof (short) * (ncols+2) ) ;
    if (prob)  for ( j=0;j<=ncols+1;j++)  if ( rand_double() <prob ) (*grid)[i][j]=1; else 	 (*grid)[i][j]=0;	
//        for (j=0;j<ncols+2;j++) (*grid)[i][j]=0;       
    }
}

/////////////////////////// save_grid ////////////////////////////////////


void save_grid(int nrows, int ncols, short *** grid, char *datafile, int rank_num , int rank_tot){
FILE *file;
int i,n=1;
double t1, t2, dt12;
struct timeval tempo ;
long datasize=nrows*ncols*sizeof(short) ;     /* tot number of bytes */
char hostname[100];
gethostname(hostname, 100);   

gettimeofday(&tempo,0);  t1=tempo.tv_sec+(tempo.tv_usec/1000000.0); // start timer
file = fopen(datafile,"w+");

fwrite(&n,sizeof(int),1,file);     // chunk num
fwrite(&n,sizeof(int),1,file);     // chunk tot
fwrite(&nrows,sizeof(int),1,file); // num rows
fwrite(&ncols,sizeof(int),1,file); // num cols

   for(i=1;i<=nrows;i++)
          { fwrite((*grid)[i],sizeof(int),ncols,file); }

fclose(file);

gettimeofday(&tempo,0); t2=tempo.tv_sec+(tempo.tv_usec/1000000.0); // stop timer
dt12=t2-t1;

printf("#W-Bytes \t time(sec)  \t MBytes/s \t ExecutionHost\n");
printf("%ld   \t %2.6f  \t %f \t %s \n", datasize , dt12, 1.0e-6 * datasize /dt12, hostname);

}

/////////////////////////// load_grid ////////////////////////////////////

int load_grid(int nrows, int ncols, short *** grid, char *datafile ,int rank_num, int rank_tot){
int num,tot, rows, cols;
char s[100]; int n;
int i,j;
double t1, t2, dt12;
struct timeval tempo ;
long datasize=nrows*ncols*sizeof(short) ;     /* tot number of bytes */
char hostname[100];
gethostname(hostname, 100);

gettimeofday(&tempo,0);  t1=tempo.tv_sec+(tempo.tv_usec/1000000.0); // start timer

FILE *file;
if ( (file = fopen(datafile,"rb") ) == NULL ) return 0;

fread( &num,sizeof(int),1,file);
fread( &tot,sizeof(int),1,file);
fread( &rows,sizeof(int),1,file);
fread( &cols,sizeof(int),1,file);
if (DEBUG==1) fprintf(stderr, "#Reading data n:%d t:%d r:%d c:%d \n", num, tot, rows, cols);

if ( (rows!=nrows) || (cols!=ncols) ) { fclose(file); printf ("dimension mismatch\n"); return 0;}

    for (i=1;i<=nrows;i++)
          { fread( (*grid)[i],sizeof(int),ncols,file);
                  if (DEBUG>1) {for(j=1;j<=ncols;j++) printf ("%d ", (*grid)[i][j]);  printf("\n");}
          }
                  fclose(file);

gettimeofday(&tempo,0); t2=tempo.tv_sec+(tempo.tv_usec/1000000.0); // stop timer
dt12=t2-t1;
printf("#R-Bytes \t time(sec)  \t MBytes/s \t ExecutionHost\n");
printf("%ld   \t %2.6f \t %f \t %s\n", datasize , dt12, 1.0e-6 * datasize /dt12, hostname);

return 1;

}

/////////////////////////// randomize_grid ////////////////////////////////////

void randomize_grid(int nrows, int ncols, short ** grid, double prob){
   int i,j; 
   random_initByTime(rand()%100) ;
    
    for ( i=1;i<=nrows;i++) 
        for ( j=1;j<=ncols;j++) 
            if ( rand_double() <prob ) grid[i][j]=1; 
}

/////////////////////////// rand_double ////////////////////////////////////

double rand_double() {
    return (double)rand()/(double)RAND_MAX;
}

/////////////////////////// random_iniyByTime ////////////////////////////////////

void random_initByTime(int rank) {
    time_t ltime;
    time(&ltime);
    srand((unsigned) ltime + 100*rank);
}

/////////////////////////// clerascreen ////////////////////////////////////

void clearscreen() 
{
if (system( "clear" )) system( "cls" );
}


