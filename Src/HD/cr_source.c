#include "pluto.h"


double CR_Source(int swp, const State_1D *state, double **vp,\
             double **vm, double dtdx, double dtdV, double *A)
/*!
 * Purpose : Adding coupling term for CR-HD
 *
 * \param [in] swp     Runs along the direction
 * \param [in] state   To access CR pressure, and hll speed.
 * \param [in] vp, vm  To access cell centered variables
 * \param [in] dtdx    dt . gradient
 * \param [in] dtdV    dt. divergence
 * \param [in] A       Interface area
 *
 * \method             NC_VdP_TOTENG => ENG eq. (eth+eke+ecr) and
 *                                       source term is vdp in ECR eq.
 *                     NC_PdV => ENG eq. (eth+eke) and
 *                               source term is pdv in both ENG and ECR eq.
 *
 * \g_inputParam[OPT]: Options are given to check different implementation.
 *                     (To be removed in final version)
 *
 * \returns            CR source term 
 *
 * \Author: Siddhartha Gupta (gsiddhartha@uchicago.edu)
 * Reference: Eq. 29 and 30 in Gupta, Sharma, Mignone 2021.
 *********************************************************************** */
{
 double CR_HD_Source;
 double *pcr  = state->presscr;
    
#if (CR_FLUID == NC_PdV_TOTENG)
   double vl, vr, Pcr;
   
      vl  = state->u_riemann[swp-1][VXn];
      vr  = state->u_riemann[swp][VXn];
   
      Pcr = 0.5*(state->u_riemann[swp-1][PCR]+state->u_riemann[swp][PCR]);          
    
   #if GEOMETRY == CARTESIAN
     CR_HD_Source = (vr-vl) * dtdx *  Pcr;
   #else 
     CR_HD_Source = (A[swp]*vr-A[swp-1]*vl)* dtdV * Pcr ;
   #endif
   
#endif
 
 
  return(CR_HD_Source);
}


