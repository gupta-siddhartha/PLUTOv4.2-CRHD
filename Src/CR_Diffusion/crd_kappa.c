/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the cosmic ray diffusion coefficient.
  \date    June 04, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void CRD_kappa(double *v, double x1, double x2, double x3,
              double *kpar, double *knor, double *phi)
/*!
 * Compute thermal conduction coefficients.
 *
 * \param [in] v   array of primitive variables
 * \param [in] x1  coordinate in the X1 direction
 * \param [in] x2  coordinate in the X2 direction
 * \param [in] x3  coordinate in the X3 direction
 * \param [out] kpar pointer to the conduction coefficient 
 *                   \f$ \kappa_\parallel \f$ in the direction of magnetic
 *                     field
 * \param [out] knor pointer to the conduction coefficient 
 *                   \f$ \kappa_\perp \f$ perpendicular to magnetic field 
 * \param [out] phi    pointer to the parameter \f$ \phi \f$ controlling the
 *                     magnitude of the saturated flux.
 *    
 *********************************************************************** */
{
 
}
