/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLL Riemann solver for HD.

  Solve the Riemann problem for the HD equations using the
  single-state HLL solver by Toro.

  On input, this function takes left and right primitive state vectors
  \c state->vL and \c state->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface
  \c i+1/2 (note that the \c i refers to \c i+1/2).

  Also during this step, compute maximum wave propagation speed (cmax)
  for  explicit time step computation.

  \b Reference:
   - "Riemann Solver and Numerical Methods for Fluid Dynamics"
      by E.F. Toro (Chapter 10)
   - "On Godunov-Type Method near Low Densities"
     by B. Einfeldt, C.D. Munz, P.L. Roe, JCP 92, 273-295 (1991)

  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 23, 2012

  \Changes: HLL function includes cosmic rays
  \Reference: Gupta, Sharma, Mignone 2021.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLL_Solver (const State_1D *state, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic/isothermal MHD equations 
 * using the HLL Riemann solver.
 *
 * \param[in,out] state   pointer to State_1D structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int    nv, i;
  double   scrh;
  static double *pL, *pR, *SL, *SR, *a2L, *a2R;
  static double **fL, **fR;
  double *uR, *uL;
  double bmax, bmin, *vL, *vR, aL, aR;

  static double *pcrL, *pcrR; // NEW

/* -- Allocate memory -- */

  if (fL == NULL){
    fL = ARRAY_2D(NMAX_POINT, NFLX, double);
    fR = ARRAY_2D(NMAX_POINT, NFLX, double);

    pR = ARRAY_1D(NMAX_POINT, double);
    pL = ARRAY_1D(NMAX_POINT, double);
    SR = ARRAY_1D(NMAX_POINT, double);
    SL = ARRAY_1D(NMAX_POINT, double);

    a2R = ARRAY_1D(NMAX_POINT, double);
    a2L = ARRAY_1D(NMAX_POINT, double);

    pcrR = ARRAY_1D(NMAX_POINT, double); // NEW
    pcrL = ARRAY_1D(NMAX_POINT, double); // NEW
  }

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, pcrL, beg, end); // NEW
  Flux (state->uR, state->vR, a2R, fR, pR, pcrR, beg, end); // NEW

  HLL_Speed (state->vL, state->vR, a2L, a2R, SL, SR, beg, end);
  
  double uhll_MXn, uhll_RHO; //NEW
  
  for (i = beg; i <= end; i++) {

    scrh = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i]  =scrh;

    if (SL[i] > 0.0){
    
      for (nv = NFLX; nv--; ) state->flux[i][nv] = fL[i][nv];
      state->press[i] = pL[i];

      #if CR_FLUID != NO // NEW
      state->presscr[i] = pcrL[i];
      state->u_riemann[i][RHO] = state->vL[i][RHO];
      state->u_riemann[i][VXn] = state->vL[i][VXn];
      state->u_riemann[i][PCR] = state->vL[i][PCR];
      #endif

    }else if (SR[i] < 0.0){

      for (nv = NFLX; nv--; ) state->flux[i][nv] = fR[i][nv];
      state->press[i] = pR[i];

      #if CR_FLUID != NO // NEW
      state->presscr[i] = pcrR[i];
      state->u_riemann[i][RHO] = state->vR[i][RHO];
      state->u_riemann[i][VXn] = state->vR[i][VXn]; 
      state->u_riemann[i][PCR] = state->vR[i][PCR];
      #endif

    }else{

      uR = state->uR[i];
      uL = state->uL[i];

      scrh = 1.0 / (SR[i] - SL[i]);
      for (nv = NFLX; nv--;  ) {
        state->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                             SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
        state->flux[i][nv] *= scrh;
      }
      state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;

      #if CR_FLUID != NO // NEW
      state->presscr[i] = (SR[i]*pcrL[i] - SL[i]*pcrR[i])*scrh;  //not used               
     
      uhll_MXn = scrh*( SR[i]*uR[MXn] - SL[i]*uL[MXn] + (fL[i][MXn]-fR[i][MXn])  \
                                               +  (pL[i] - pR[i]) ) ;
                                                              
      uhll_RHO = scrh*(SR[i]*uR[RHO] - SL[i]*uL[RHO] + (fL[i][RHO]-fR[i][RHO]) );
          
      state->u_riemann[i][RHO] = uhll_RHO;
      state->u_riemann[i][VXn] = uhll_MXn / uhll_RHO;
      state->u_riemann[i][PCR] = scrh*( SR[i]*uR[ECR] - SL[i]*uL[ECR] + (fL[i][ECR]-fR[i][ECR]) )*(g_gammacr-1.);
     #endif 
     
    }
  } /* end loops on points */
}
