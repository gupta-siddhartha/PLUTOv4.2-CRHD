/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief HLLC Riemann solver for MHD.

  Solve the Riemann problem for the HD equations using the 
  two-state HLLC solver by Toro.
  
  On input, this function takes left and right primitive state vectors 
  \c state->vL and \c state->vR at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
   -   "Riemann Solver and Numerical Methods for Fluid Dynamics"
        by E.F. Toro (Chapter 10)
       
  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 23, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/* ********************************************************************* */
void HLLC_Solver (const State_1D *state, int beg, int end, 
                  real *cmax, Grid *grid)
/*!
 * Solve Riemann problem using the HLLC Riemann solver.
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
  double scrh, vxr, vxl;
  double usL[NFLX], usR[NFLX], vs;
  double qL, qR, wL, wR;
  double *vL, *vR, *uL, *uR;
  #if EOS == ISOTHERMAL
   double rho, mx;
  #endif
  static double *pL, *pR, *SL, *SR, *a2L, *a2R;
  static double **fL, **fR;

  static double *pcrL, *pcrR; // NEWi
  double uhll_MXn, uhll_RHO; //NEW
  #if CR_FLUID != NO && CR_FLUID == NC_PdV_TOTENG
  print1(" >> Error: HLLC solver does not fit well with CRs.");
  print1(" >> Choose solver 'hll\n'");
  #endif
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
    
    // WHERE TO FREE ?
    pcrR = ARRAY_1D(NMAX_POINT, double); // NEW
    pcrL = ARRAY_1D(NMAX_POINT, double); // NEW
  }

/* ----------------------------------------------------
    Compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (state->vL, a2L, NULL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (state->vR, a2R, NULL, beg, end, FACE_CENTER, grid);

  Flux (state->uL, state->vL, a2L, fL, pL, pcrL, beg, end);  // NEW
  Flux (state->uR, state->vR, a2R, fR, pR, pcrR, beg, end);  // NEW

  HLL_Speed (state->vL, state->vR, a2L, a2R, SL, SR, beg, end);

  for (i = beg; i <= end; i++) {

    scrh = MAX(fabs(SL[i]), fabs(SR[i]));
    cmax[i]  = scrh;

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

      vR = state->vR[i]; uR = state->uR[i];
      vL = state->vL[i]; uL = state->uL[i];

      vxr = vR[VXn];
      vxl = vL[VXn];

     /*#if CR_FLUID != NO
     vR[PRS] += vR[PCR];
     vL[PRS] += vL[PCR];
     #endif*/
 
#if SHOCK_FLATTENING == MULTID   
      if ((state->flag[i] & FLAG_HLL) || (state->flag[i+1] & FLAG_HLL)){        
         scrh  = 1.0/(SR[i] - SL[i]);
         for (nv = NFLX; nv--; ){
           state->flux[i][nv]  = SL[i]*SR[i]*(uR[nv] - uL[nv])
                              +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
           state->flux[i][nv] *= scrh;
         }
         state->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
         
         #if CR_FLUID != NO // NEW
         state->presscr[i] = (SR[i]*pcrL[i] - SL[i]*pcrR[i])*scrh;
         #endif
         
         continue;
       }
#endif

  /* ---------------------------------------
                   get u* 
     --------------------------------------- */    

      #if HAVE_ENERGY
       qL = vL[PRS] + uL[MXn]*(vL[VXn] - SL[i]);
       qR = vR[PRS] + uR[MXn]*(vR[VXn] - SR[i]);

       wL = vL[RHO]*(vL[VXn] - SL[i]);
       wR = vR[RHO]*(vR[VXn] - SR[i]);

       #if CR_FLUID != NO //NEW
       qL += vL[PCR]; 
       qR += vR[PCR];
       #endif

       vs = (qR - qL)/(wR - wL); /* wR - wL > 0 since SL < 0, SR > 0 */
/*
      vs = vR[PRS] - vL[PRS] + uL[MXn]*(SL[i] - vxl) 
                           - uR[MXn]*(SR[i] - vxr);
      vs /= vL[RHO]*(SL[i] - vxl) - vR[RHO]*(SR[i] - vxr);
*/

       usL[RHO] = uL[RHO]*(SL[i] - vxl)/(SL[i] - vs);
       usR[RHO] = uR[RHO]*(SR[i] - vxr)/(SR[i] - vs);
       EXPAND(usL[MXn] = usL[RHO]*vs;     usR[MXn] = usR[RHO]*vs;      ,
              usL[MXt] = usL[RHO]*vL[VXt]; usR[MXt] = usR[RHO]*vR[VXt];  ,
              usL[MXb] = usL[RHO]*vL[VXb]; usR[MXb] = usR[RHO]*vR[VXb];)
           
       usL[ENG] =    uL[ENG]/vL[RHO] 
                  + (vs - vxl)*(vs + vL[PRS]/(vL[RHO]*(SL[i] - vxl)));
       usR[ENG] =    uR[ENG]/vR[RHO] 
                  + (vs - vxr)*(vs + vR[PRS]/(vR[RHO]*(SR[i] - vxr)));

       #if CR_FLUID != NO // NEW
       usL[ENG] += (vs - vxl)*(vL[PCR]/(vL[RHO]*(SL[i] - vxl)));
       usR[ENG] += (vs - vxr)*(vR[PCR]/(vR[RHO]*(SR[i] - vxr)));
       
       usL[ECR] = uL[ECR]*(SL[i] - vxl)/(SL[i] - vs);
       usR[ECR] = uR[ECR]*(SR[i] - vxr)/(SR[i] - vs);       
       #endif
       
       usL[ENG] *= usL[RHO];
       usR[ENG] *= usR[RHO]; 
       
      #elif EOS == ISOTHERMAL
       scrh = 1.0/(SR[i] - SL[i]);
       rho  = (SR[i]*uR[RHO] - SL[i]*uL[RHO] - fR[i][RHO] + fL[i][RHO])*scrh;
       mx   = (SR[i]*uR[MXn] - SL[i]*uL[MXn] - fR[i][MXn] + fL[i][MXn])*scrh;
       
       usL[RHO] = usR[RHO] = rho;
       usL[MXn] = usR[MXn] = mx;
       vs  = (  SR[i]*fL[i][RHO] - SL[i]*fR[i][RHO] 
              + SR[i]*SL[i]*(uR[RHO] - uL[RHO]));
       vs *= scrh;
       vs /= rho;
       EXPAND(                                            ,
              usL[MXt] = rho*vL[VXt]; usR[MXt] = rho*vR[VXt]; ,
              usL[MXb] = rho*vL[VXb]; usR[MXb] = rho*vR[VXb];)
      #endif

        
 /* ---- verify consistency condition ------ */
 /*
      for (nv = 0; nv < NFLX; nv++) {
        scrh  = (vs - SL[i])*usL[nv] + (SR[i] - vs)*usR[nv];
        scrh -= SR[i]*uR[nv] - SL[i]*uL[nv] + 
                fL[i][nv] - fR[i][nv];
          if (fabs(scrh) > 1.e-9){
          printf (" Consistency condition violated\n");
          printf ("%d %d  %12.6e\n",i,nv, scrh);
        }
      }
 */          
/*  ----  Compute HLLC flux  ----  */

      if (vs >= 0.0){
        for (nv = NFLX; nv--;   ) {
          state->flux[i][nv] = fL[i][nv] + SL[i]*(usL[nv] - uL[nv]);
        }
        state->press[i] = pL[i];
        #if CR_FLUID != NO // NEW
        state->presscr[i] = pcrL[i];
        state->u_riemann[i][PCR] = usL[ECR]*(g_gammacr-1.); //NEW
        #endif
                  
      } else {
        for (nv = NFLX; nv--;   ) {
          state->flux[i][nv] = fR[i][nv] + SR[i]*(usR[nv] - uR[nv]);
        }
        state->press[i] = pR[i];
        #if CR_FLUID != NO // NEW
        state->presscr[i] = pcrR[i];        
        state->u_riemann[i][PCR] = usR[ECR]*(g_gammacr-1.); //NEW
        #endif
      }
      
      // Normal component of velocity is identian to star state.
      

     #if CR_FLUID != NO // NEW
     scrh = 1.0 / (SR[i] - SL[i]);
    
      //state->presscr[i] = (SR[i]*pcrL[i] - SL[i]*pcrR[i])*scrh;  //not used               
     
      //uhll_MXn = scrh*( SR[i]*uR[MXn] - SL[i]*uL[MXn] + (fL[i][MXn]-fR[i][MXn])  \
                                               +  (pL[i] - pR[i]) ) ;
                                                              
      //uhll_RHO = scrh*(SR[i]*uR[RHO] - SL[i]*uL[RHO] + (fL[i][RHO]-fR[i][RHO]) );
          
      //state->u_riemann[i][RHO] = uhll_RHO;
      state->u_riemann[i][VXn] = vs; //uhll_MXn / uhll_RHO;
      //printf("%0.4e %0.4e\n", uhll_MXn / uhll_RHO, vs);
      state->u_riemann[i][PCR] = scrh*( SR[i]*uR[ECR] - SL[i]*uL[ECR] + (fL[i][ECR]-fR[i][ECR])  )*(g_gammacr-1.);
      //state->flux[i][ECR] = scrh*(SL[i]*SR[i]*(uR[ECR] - uL[ECR]) + SR[i]*fL[i][ECR] - SL[i]*fR[i][ECR]);
     #endif 
     
    }
      
  } /* end loops on points */
}
