/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the hydro (HD) flux.                                             

  Compute the flux of the conservative HD equations in the direction 
  given by ::g_dir.
  This function defines the component of the hyperbolic flux tensor 
  of the standard HD equations.\n
  In what follows:
  - \c VXn, \c MXn are the velocity, momentum components in the direction 
    given by ::g_dir (normal, \c "n")
  - \c VXt, \c MXt and \c VXb, \c MXb are the transverse components
    (tangent \c "t" and bi-tangent \c "b").

 \author A. Mignone (mignone@ph.unito.it)
 \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (double **u, double **w, double *a2, double **fx, double *p, double *pcr, 
           int beg, int end)
/*!
 * \param [in]    u    1D array of conserved quantities
 * \param [in]    w    1D array of primitive quantities
 * \param [in]   a2    1D array of sound speeds
 * \param [out]  fx    1D array of fluxes (total pressure excluded)
 * \param [out]   p    1D array of pressure values
 * \param [out]   pcr  1D array of cosmic ray pressure // NEW
 * \param [in]   beg   initial index of computation 
 * \param [in]   end   final   index of computation
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  int   nv, ii;

  for (ii = beg; ii <= end; ii++) {
    fx[ii][RHO] = u[ii][MXn];
    EXPAND(fx[ii][MX1] = u[ii][MX1]*w[ii][VXn]; ,
           fx[ii][MX2] = u[ii][MX2]*w[ii][VXn]; ,
           fx[ii][MX3] = u[ii][MX3]*w[ii][VXn];)
#if HAVE_ENERGY
    p[ii] = w[ii][PRS];
    fx[ii][ENG] = (u[ii][ENG] + w[ii][PRS])*w[ii][VXn];
#elif EOS == ISOTHERMAL
    p[ii] = a2[ii]*w[ii][RHO];
#endif


#if CR_FLUID != NO // NEW
     p[ii] += w[ii][PCR]; // term includes CR pressure
     pcr[ii]      = w[ii][PCR];
     
    // Adds CR work done term in the energy flux i.e.,\
       fx_{eng} = ( e_{internal} + (p_{th} + p_{cr}) ) v 

     fx[ii][ENG] += w[ii][PCR]*w[ii][VXn]; 
     
     if (CR_FLUID == NC_VdP_TOTENG) { 
     //option 1: (VdP method) fx_{ecr}=(e_{cr} + p_{cr}) v; 
     
       fx[ii][ECR] = (u[ii][ECR] + w[ii][PCR])*w[ii][VXn];
       
     } else {
      //option != 1 (other methods) fx_{ecr} =  e_{cr} v;
     
       fx[ii][ECR] = u[ii][ECR]*w[ii][VXn]; 
       
     } 
#endif

   // FOR HLLI solver: Momentum flux includes pressure term 
   // which is not the case for other HLL solver.
   // fx[ii][MX1] += w[ii][PRS];
   // #if CR_FLUID != NO 
   //fx[ii][MX1] += w[ii][PCR]; //NEW1
   //#endif
    
    
/*
#if DUST == YES
    fx[ii][RHO_D] = u[ii][MXn_D];
    EXPAND(fx[ii][MX1_D] = u[ii][MX1_D]*w[ii][VXn_D]; ,
           fx[ii][MX2_D] = u[ii][MX2_D]*w[ii][VXn_D]; ,
           fx[ii][MX3_D] = u[ii][MX3_D]*w[ii][VXn_D];)
#endif
*/
  }
}
