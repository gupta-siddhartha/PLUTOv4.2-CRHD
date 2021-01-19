/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Compute the outermost wave speeds for HLL-based solvers.

  HLL_Speed() computes an estimate to the leftmost and rightmost
  wave signal speeds bounding the Riemann fan based on the input states
  ::vR and ::vL.
  Depending on the estimate, several variants are possible.
 
  \authors A. Mignone (mignone@ph.unito.it)
  \date    June 6, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

/*   ****  switches, see below for description  ****  */

#define DAVIS_ESTIMATE     YES //NO
#define EINFELDT_ESTIMATE  NO
#define ROE_ESTIMATE       NO

/* ********************************************************************* */
void HLL_Speed (double **vL, double **vR, double *a2L, double *a2R,
                double *SL, double *SR, int beg, int end)
/*!
 * Compute leftmost (SL) and rightmost (SR) speed for the Riemann fan.
 * 
 * \param [in]  vL   left  state for the Riemann solver
 * \param [in]  vR   right state for the Riemann solver
 * \param [in] a2L   1-D array containing the square of the sound speed
 *                   for the left state
 * \param [in] a2R   1-D array containing the square of the sound speed
 *                   for the right state
 * \param [out] SL   the (estimated) leftmost speed of the Riemann fan
 * \param [out] SR   the (estimated) rightmost speed of the Riemann fan
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 * Switches:
 *
 *    ROE_ESTIMATE (YES/NO), DAVIS_ESTIMATE (YES/NO). TVD_ESTIMATE (YES/NO)
 *    JAN_HLL (YES/NO) 
 *
 *    These switches set how the wave speed estimates are
 *    going to be computed. Only one can be set to 'YES', and
 *    the rest of them must be set to 'NO'  
 *
 *     ROE_ESTIMATE:    b_m = \min(0, \min(u_R - c_R, u_L - c_L, u_{roe} - c_{roe}))     
 *                      b_m = \min(0, \min(u_R + c_R, u_L + c_L, u_{roe} + c_{roe}))
 * 
 *                      where u_{roe} and c_{roe} are computed using Roe averages.
 *  
 *     DAVIS_ESTIMATE:  b_m = \min(0, \min(u_R - c_R, u_L - c_L))     
 *                      b_m = \min(0, \min(u_R + c_R, u_L + c_L))  
 *
 *********************************************************************** */
{
  int    i;
  double scrh, s, c;
  double aL, dL;
  double aR, dR;
  double dvx, dvy, dvz;
  double a_av, du, vx;
  static double *sl_min, *sl_max;
  static double *sr_min, *sr_max;

  if (sl_min == NULL){
    sl_min = ARRAY_1D(NMAX_POINT, double);
    sl_max = ARRAY_1D(NMAX_POINT, double);

    sr_min = ARRAY_1D(NMAX_POINT, double);
    sr_max = ARRAY_1D(NMAX_POINT, double);
  }

/* ----------------------------------------------
              DAVIS Estimate  
   ---------------------------------------------- */

  #if DAVIS_ESTIMATE == YES

   MaxSignalSpeed (vL, a2L, sl_min, sl_max, beg, end);
   MaxSignalSpeed (vR, a2R, sr_min, sr_max, beg, end);
   for (i = beg; i <= end; i++) {

     SL[i] = MIN(sl_min[i], sr_min[i]);
     SR[i] = MAX(sl_max[i], sr_max[i]);

     aL = sqrt(a2L[i]);
     aR = sqrt(a2R[i]);

     scrh  = fabs(vL[i][VXn]) + fabs(vR[i][VXn]);    
     scrh /= aL + aR;

     g_maxMach = MAX(scrh, g_maxMach); 
   }

  #endif

/* ----------------------------------------------
        Einfeldt Estimate for wave speed 
   ---------------------------------------------- */

  #if EINFELDT_ESTIMATE == YES

   for (i = beg; i <= end; i++) {

     aL = sqrt(a2L[i]);
     aR = sqrt(a2R[i]);

     dL   = sqrt(vL[i][RHO]);
     dR   = sqrt(vR[i][RHO]);
     a_av = 0.5*dL*dR/( (dL + dR)*(dL + dR));
     du   = vR[i][VXn] - vL[i][VXn];
     scrh = (dR*aL*aL + dR*aR*aR)/(dL + dR);
     scrh += a_av*du*du;
      
     //SL[i] = (dl*ql[VXn] + dr*qr[VXn])/(dl + dr);
  
     //bmin = MIN(0.0, um[VXn] - sqrt(scrh));
     //bmax = MAX(0.0, um[VXn] + sqrt(scrh));      
   }
  #endif

/* ----------------------------------------------
              Roe-like Estimate  
   ---------------------------------------------- */

  #if ROE_ESTIMATE  == YES
   

  double rhoL_half, rhoR_half;
  double H_av, v_av;
  double EL, ER, HL, HR;
  double vl, vr;
  double wl, wr, gl, gr, g_av;
  

   for (i = beg; i <= end; i++) {

     aL = sqrt(a2L[i]);
     aR = sqrt(a2R[i]);


     rhoL_half = sqrt(vL[i][RHO]);
     rhoR_half = sqrt(vR[i][RHO]);
     

     // v_av
     v_av  = rhoL_half * vL[i][VXn] + rhoR_half * vR[i][VXn];
     v_av /= (rhoL_half + rhoR_half);


     // a_av;
     vl = EXPAND(vL[i][VXn]*vL[i][VXn],, );
     vr = EXPAND(vR[i][VXn]*vR[i][VXn],, );
  
     EL = 0.5 * vL[i][RHO] * vl * vl + vL[i][PRS]/(g_gamma-1.) + vL[i][PCR]/(g_gammacr-1.);
     ER = 0.5 * vR[i][RHO] * vr * vr + vR[i][PRS]/(g_gamma-1.) + vR[i][PCR]/(g_gammacr-1.);
    

 
     HL = (EL + vL[i][PRS] + vL[i][PCR])/vL[i][RHO];
     HR = (ER + vR[i][PRS] + vR[i][PCR])/vR[i][RHO];


     H_av = (rhoL_half * HL +  rhoR_half * HR)  / (rhoR_half + rhoL_half); 


    // g_Av;
    wl = vL[i][PCR]/(vL[i][PRS] + vL[i][PCR]);  wr = vR[i][PCR]/(vR[i][PRS] + vR[i][PCR]);

    gl = (5. + 3 * wl)/(3. * (1. + wl)); gr = (5. + 3 * wr)/(3.*(1. + wr)) ;

    g_av = (rhoL_half * gl + rhoR_half * gr)/(rhoR_half + rhoL_half);

     a_av = sqrt( (g_av - 1.)*(H_av  - .5*v_av) );
 

     SL[i] =  v_av -  a_av;
     SR[i] =  v_av +  a_av;

     scrh = (fabs(vL[i][VXn]) + fabs(vR[i][VXn]))/(aL + aR);
     g_maxMach = MAX(scrh, g_maxMach);

   }
  #endif
   
}
#undef DAVIS_ESTIMATE     
#undef EINFELDT_ESTIMATE  
#undef ROE_ESTIMATE       
