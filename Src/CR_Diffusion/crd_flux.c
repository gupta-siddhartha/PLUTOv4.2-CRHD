/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the Cosmic ray diffusion flux.

  Compute the Cosmic ray diffusion flux
 
  \date    Aug 24, 2015
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CRD_Flux (double ***Ecr, const State_1D *state,  
              double **dcoeff, int beg, int end, Grid *grid)
/*! 
 * Compute the Cosmic ray diffusion flux, state->par_flx.
 *
 * \param [in]     Ecr       3D array containing the dimensionless 
 *                         temperature
 * \param [in,out] state   pointer to a State_1D structure
 * \param [out]    dcoeff  the diffusion coefficient needed for computing
 *                         the time step.
 * \param [in]      beg   initial index of computation
 * \param [in]      end   final   index of computation
 * \param [in]      grid  pointer to an array of Grid structures
 *
 * \return This function has no return value.                       
 *********************************************************************** */
{
  int  i, j, k, nv;
  double bgradEcr, Bmag, dEcrmag;
  double Fc, Fcmag, Fsat, g1 = g_gammacr - 1.0;
  double x1, x2, x3;
  double alpha, uL, uR, suL, suR, bn;
  double vi[NVAR], kpar=0.0, knor=0.0, phi;
  double bck_fld[3];
  static double **gradEcr;

/* -----------------------------------------------------------
   1. Allocate memory, compute temperature gradient in the
      required direction and set 2 out of the 3 coordinates.
   ----------------------------------------------------------- */

  if (gradEcr == NULL) {
    gradEcr = ARRAY_2D(NMAX_POINT, 3, double);
  }

  GetGradientECR (Ecr, gradEcr, beg, end, grid);
  if (g_dir == JDIR || g_dir == KDIR) x1 = grid[IDIR].x[g_i];
  if (g_dir == IDIR || g_dir == KDIR) x2 = grid[JDIR].x[g_j];
  if (g_dir == IDIR || g_dir == JDIR) x3 = grid[KDIR].x[g_k];

/* ----------------------------------------------- 
   2. Compute Cosmic ray diffusion Flux (crdflx).
   ----------------------------------------------- */

  for (i = beg; i <= end; i++){
    
    for (nv = 0; nv < NVAR; nv++) {
      vi[nv]  = 0.5*(state->vh[i][nv] + state->vh[i+1][nv]);
    }

  /* ------------------------------------------------
     2 Calculate the Cosmic ray diffusion coefficient
     ------------------------------------------------ */

    if (g_dir == IDIR) x1 = grid[IDIR].xr[i];
    if (g_dir == JDIR) x2 = grid[JDIR].xr[i];
    if (g_dir == KDIR) x3 = grid[KDIR].xr[i];

    
    CRD_kappa(vi, x1, x2, x3, &kpar, &knor, &phi);
    dEcrmag  = D_EXPAND(  gradEcr[i][0]*gradEcr[i][0], 
                      + gradEcr[i][1]*gradEcr[i][1], 
                      + gradEcr[i][2]*gradEcr[i][2]);
    dEcrmag = sqrt(dEcrmag) + 1.e-12;
 
#if PHYSICS == HD 

    Fc = kpar*gradEcr[i][g_dir];  
  
    state->crd_flux[i][ECR] = Fc;
    dcoeff[i][ECR]         = fabs(kpar);

#endif
  }
}


