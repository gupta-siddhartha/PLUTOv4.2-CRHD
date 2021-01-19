/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Cosmic Ray diffusion (CRD) module header file
*/
/* ///////////////////////////////////////////////////////////////////// */

void   CRD_Flux (double ***, const State_1D *, double **, int, int, Grid *);
void   CRD_kappa (double *, double, double, double, double *, double *, double *);
void   GetGradientECR (double ***, double **, int, int, Grid *);
