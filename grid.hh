using namespace std;

class grid {

/* variables */

public:
  int dim;              // Dimension of the grid (can be 1,2,3)
  int N[3];             // Computational grid size in the direction (x,y,z)
  int S[3];             // Size of the stencil  (# of points) in the diriection (x,y,z)
  int Nvars;            // Number of grid variables
  double **vars;        // array of pointers to the allocated grid variables
  static int DEBUG;     // debug flag for class-specific operations

private:
  // ------------------ Local Grid dimensions ---------------------
  int V[3];             // Vectorization size in direction (x,y,z)
                        //  --->  This allow that the allocated dimension in the direction i
                        //        is a multiple of V[i]. Default is V[_]=1 (scalar)
  int D[3];             // Real allocated grid size: D[i] = N[i] + 2 * S[i]
                        // If vettorization is request .... ((N[i]-1)/V[i]+1)*V[i]
                        //                                  instead of N[i]
  int SIZE;             // Total number of points to be allocated = D[0]*D[1]*D[2]
  // ----------------------------------------------------------------
  // Useful data for index conversion that it is better to store
  // ----------------------------------------------------------------
  int origin;           // Origin of the first point to be computed
  int STRIDEx;          // Next in direction x => +1
  int STRIDEy;          // Next in direction y => + D[0]
  int STRIDEz;          // Next in direction z => + D[0]*D[1]

/* operators */

public:
  double& operator()(int var,int idx)            {return vars[var][idx];}  // get element by linear index
  //double& operator()(int var,int i,int j)        {return vars[var][origin + i + STRIDEy * j ];}               // get element by coordinates
  double& operator()(int var,int i,int j, int k) {return vars[var][origin + i + STRIDEy * j + STRIDEz * k];}  // get element by coordinates

/* contructors */

public:
  grid(int NVARS);
  grid(int Nx,int Ny,int Nz,int Sx,int Sy,int Sz,int NVARS,double initval=1.0,double initval2=2.0);

  grid(const grid &grid2); // copy constr. (by const ref.)

/* other functions */

public:
  // Function to read private member
  int size() const {return SIZE;};
  // -----------------------------------------------------------------------
  // Inline function that returns the correct index given the linear one
  // -----------------------------------------------------------------------
  int idx_x(int i)  const { return  ( i % D[0] ) - S[0] ;};
  int idx_y(int i)  const { return  ( i / D[0] ) % D[1] -S[1] ;};
  int idx_z(int i)  const { return  ( i / D[0] ) / D[1] -S[2] ;};
  int idx_ijk(int i,int j,int k)  const {return origin + IDX_ijk(i,j,k);};
  int IDX_ijk(int i,int j,int k)  const {return (i + STRIDEy * j + STRIDEz * k  ); };

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
    origin = S[0] + D[0] * S[1] + D[0]*D[1] *S[2] ;
    STRIDEx = 1;
    STRIDEy = D[0];
    STRIDEz = D[0]*D[1];
    return *this;
  };

  friend void swap_grids(grid &G1, grid &G2);
  friend void swap_grids_pointers(grid &G1, grid &G2);

};
