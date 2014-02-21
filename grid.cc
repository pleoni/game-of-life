#include <iostream>
#include <stdio.h>
#include <cstdlib>

#include "grid.hh"


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
// General constructor for a 3d GRID structure (It also allocates memory)
//
// It also initializes values with two different values depending if the variables
// are on the compuational grid or in the stencil part of the array
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
    vars = new double*[Nvars]; // allocate the array
    for(int m=0; m<Nvars; m++) {
      // vars[i]=malloc(sizeof(double)*SIZE);
      vars[m]=new double[SIZE];  // allocate the actual memory for data
      if (DEBUG==1) printf("Variable %d: Allocated size %d\n",m,SIZE);
      for(int i=0;i<SIZE;i++) {
        vars[m][i] = initval*(m+1.0);  // set all points
      }
      for(int i=0;i<N[0];i++) {
        for(int j=0;j<N[1];j++) {
          for(int k=0;k<N[2];k++) {
            (*this)(m,i,j,k) = initval2;  // only set points in the "real" grid
          };
        };
      };
    }; // for m=0:Nvars
};


grid::grid(const grid &g) { // copy constructor (by const reference)

  dim = g.dim;
  Nvars = g.Nvars;
  for (int i=0; i<3; i++) {
    N[i] = g.N[i];
    S[i] = g.S[i];
    V[i] = g.V[i];
    D[i] = g.D[i];
  }
  SIZE = g.SIZE;
  STRIDEx = g.STRIDEx;
  STRIDEy = g.STRIDEy;
  STRIDEz = g.STRIDEz;
  origin = g.origin; // ???

  vars = new double*[Nvars]; // allocate the array
  for(int m=0; m<Nvars; m++) {
    vars[m]=new double[SIZE];  // allocate the actual memory for data
    if (DEBUG==1) printf("(Copy constructor) Variable %d: Allocated size %d\n",m,SIZE);
    for(int i=0; i<SIZE; i++) {
      vars[m][i] = g.vars[m][i];  // copy data
    }
  }

}


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

int grid::DEBUG=0;

void swap_grids(grid G1, grid G2) {
  // TODO
  printf("WARNING: swap_grids() is not implemented yet!");
}

void swap_grids_pointers(grid G1, grid G2) {
  // Watch out: this function only swaps the pointers to the data.
  // All other variables must be equal for this to make sense!
  double ** temp;
  temp = G1.vars;
  G1.vars = G2.vars;
  G2.vars = temp;
}
