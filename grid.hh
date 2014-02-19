using namespace std;

class grid {
public:
  int dim;              // Dimension of the grid (can be 1,2,3)
  int N[3];             // Computational grid size in the diriection (x,y,z)
  int S[3];             // Size of the stencil  (# of points) in the diriection (x,y,z)
  int Nvars;            // Number of grid variables       
  double **vars;         // array of pointer to the allocated grid variables
public: 
  double& operator()(int var,int idx)            {return vars[var][idx];}
  double& operator()(int var,int i,int j)        {return vars[var][origin + i + STRIDEy * j ];}
  double& operator()(int var,int i,int j, int k) {return vars[var][origin + i + STRIDEy * j + STRIDEz * k];}
public: 
  // Function to read private member
  int size() const {return SIZE;}; 
  // -----------------------------------------------------------------------
  // Inline function that returns the correct index given the linear one
  // -----------------------------------------------------------------------
  int idx_x(int i)  const { return  ( i % D[0] ) - S[0]  ;};
  int idx_y(int i)  const { return  ( i / D[0] ) % D[1] -S[1]  ;};
  int idx_z(int i)  const { return  ( i / D[0] ) / D[1] -S[2] ;};
  int idx_ijk(int i,int j,int k)  const {return origin + IDX_ijk(i,j,k);};
  int IDX_ijk(int i,int j,int k)  const {return (i + STRIDEy * j + STRIDEz * k  ); }; 
private:
  // ------------------ Local Grid dimensions ---------------------
  int V[3];             // Vectorization size in direction (x,y,z) 
                        //  --->  This allow that the allocated dimension in the direction i
                        //        is a multiple of V[i]. Default is V[_]=1 (scalar)
  int D[3];             // Real allocated grid size: D[i] = N[i] + 2 * S[i]
                        // If vettorization is request .... ((N[i]-1)/V[i]+1)*V[i]
                        //                                  instead of N[i] 
  int SIZE;             // Total number of point to be allocate = D[0]*D[1]*D[2]
  // ----------------------------------------------------------------
  // Useful data for index conversion that it is beter to store 
  // ----------------------------------------------------------------
  int origin;           // Origin of the first point to be computed
  int STRIDEx;          // Next in direction x => +1
  int STRIDEy;          // Next in direction y => + D[0]
  int STRIDEz;          // Next in direction z => + D[0]*D[1]
public: 
  void dump();
  // *** --- SET SIZE --- ***
  grid& set_size(int Nx,int Ny,int Nz) {dim=3;N[0]=Nx;N[1]=Ny;N[2]=Nz; this->set_stencil(1);return *this;}
  grid& set_size(int Nx,int Ny)        {dim=2;N[0]=Nx;N[1]=Ny;N[2]=1;  this->set_stencil(1);return *this;}
  grid& set_size(int Nx)               {dim=1;N[0]=Nx;N[1]=N[2]=1;     this->set_stencil(1);return *this;}
  // *** --- vectorization setting --- ***
  grid& set_vectorization(int Vx,int Vy,int Vz) {V[0]=Vx;V[1]=Vy;V[2]=Vz; this->compute_val(); return *this;}
  // *** --- SET STENCIL --- ***
  grid& set_stencil(int Sx,int Sy,int Sz) {S[0]=Sx;S[1]=Sy;S[2]=Sz;    this->compute_val(); return *this;}
  grid& set_stencil(int Stencil) {S[0]=S[1]=S[2]=0;for (int i=0;i<dim;i++) S[i]=Stencil; this->compute_val(); return *this;}
  // *** --- ENSURE THAT PRIVATE VARIABLES ARE CONSISTENT --- ***
  grid& compute_val() {
    D[0] = ( ( ( N[0] -1 ) / V[0] ) + 1 ) * V[0] +2*S[0];
    D[1] = ( ( ( N[1] -1 ) / V[1] ) + 1 ) * V[1] +2*S[1];
    D[2] = ( ( ( N[2] -1 ) / V[2] ) + 1 ) * V[2] +2*S[2];
    // Setting total size
    SIZE = D[0]*D[1]*D[2];
    origin = S[0]+ D[0] * S[1] + D[0]*D[1] *S[2] ; 
    STRIDEx = 1;
    STRIDEy = D[0];
    STRIDEz = D[0]*D[1];
    return *this;
  };
public:  // CONSTRUCTOR for the class
  grid(int NVARS);
  grid(int Nx,int Ny,int Nz,int Sx,int Sy,int Sz,int NVARS,double initval=1.0,double initval2=2.0);

};


// ****************************************************************************
//
// To be moved in grid.cc
//
// ****************************************************************************

#include <iostream>


// -----------------------------------------------------------------------
// Constructor for a VOID structure for a given number of vars
// -----------------------------------------------------------------------
grid::grid(int NVARS) {
  dim = 0;
  N[0]=N[1]=N[2]=1;
  V[0]=V[1]=V[2]=1;
  S[0]=S[1]=S[2]=0;
  this->compute_val();
  Nvars = NVARS;
  vars = new double*[Nvars];
};

// ------------------------------------------------------------------------------
// General constructor for a 3d GRID structure (It also allocated memory)
//
// It also initialize values with two different values depending if the variables
// are on the compuational grid or in teh stancil part of teh array 
// --------------------------------------------------------------------------------
grid::grid(int Nx,int Ny,int Nz,int Sx,int Sy,int Sz,int NVARS,double initval,double initval2) {

    dim = 3;
    // Setting x-direction
    V[0] = 1;
    N[0] = Nx;
    S[0] = Sx;
    // Setting y-direction
    V[1] = 1;
    N[1] = Ny;
    S[1] = Sy;
    // Setting x-direction
    V[2] = 1;
    N[2] = Nz;
    S[2] = Sz;
    this->compute_val();
    Nvars = NVARS;
    vars = new double*[Nvars];
    for(int m=0;m< Nvars; m++) {
      // vars[i]=malloc(sizeof(double));
      printf("Variable %d: Allocated size %d\n",m,SIZE);
      vars[m]=new double[SIZE];
      for(int i=0;i<SIZE;i++) {
        vars[m][i] = initval*(m+1.0);
      }
      for(int i=0;i<N[0];i++) {
        for(int j=0;j<N[1];j++) {
          for(int k=0;k<N[2];k++) {
            (*this)(m,i,j,k) = initval2;
          };
        };
      };
    };
};


void grid::dump() {
  cout << "Dumping data about grid class" << endl ;
  cout << "dimension is: " << dim  << endl ;
  cout << " # of points:" << N[0] <<  "," << N[1] << "," << N[2] << endl ;
  cout << " stencil:    " << S[0] <<  "," << S[1] << "," << S[2] << endl ;
  cout <<  "allocated:  " << D[0] <<  "," << D[1] << "," << D[2] << endl ;
  cout << "STRIDEs      " << STRIDEx <<  "," << STRIDEy <<  "," << STRIDEz << endl; 
  cout << "SIZE   = " << SIZE   << endl ;
  cout << "Origin = " << origin << endl ;
}

