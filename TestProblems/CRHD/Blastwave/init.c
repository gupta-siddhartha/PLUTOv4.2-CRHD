/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * Reference : Chevalier, R. A. 1983, ApJ, 272,765
 *
 *********************************************************************** */
{ 
  double dr, vol, r;

   dr = g_domBeg[IDIR] + 1.0*(g_domEnd[IDIR]-g_domBeg[IDIR])/(double)NX1;

  double rinj = dr;//0.01;//2.*dr;

  vol = 4.0/3.0*CONST_PI*rinj*rinj*rinj;

#if GEOMETRY == CARTESIAN
  r = sqrt(x1*x1 + x2*x2 + x3*x3);
#elif GEOMETRY  == SPHERICAL
  r = x1; // Spherical geometry
#endif

  v[RHO] = g_inputParam[RHO_AMB]/UNIT_DENSITY;
  EXPAND(v[VX1] = 0.0;, v[VX2] = 0.0;, v[VX3] = 0.0;)
  v[PRS] = (g_inputParam[RHO_AMB]/CONST_mH/0.6)*CONST_kB*g_inputParam[TEMP_AMB]/unitPRS;
  
  #if CR_FLUID != NO
  v[PCR] =v[PRS]; // or any small value (e.g, 1.e-10);
  #endif
//printf("%e", rinj); exit(1);
  if (r <= rinj ) {
   v[PRS] += (g_gamma - 1.0)*(g_inputParam[E_SN]/unitENERGY)/vol;
  }         
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;


  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){

      if(d->Vc[PRS][k][j][i] < 1e-10){d->Vc[PRS][k][j][i] = 1e-10;}
  #if CR_FLUID != NO
      if(d->Vc[PCR][k][j][i] < 1e-10){d->Vc[PCR][k][j][i] = 1e-10;}
  #endif
    }
  }
}


