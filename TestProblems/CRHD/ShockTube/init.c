
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 * Reference : Table 1 in  Kudoh, Y. & Hanawa, T. 2016, MNRAS, 462, 4517
 *********************************************************************** */
{
 #if GEOMETRY  == CARTESIAN
  if(x1 < 0.0){
 #elif GEOMETRY  == SPHERICAL
  if(x1 < 0.4){
 #endif

    v[RHO] = 1.0;
    v[PRS] = 2.0;

    #if CR_FLUID != NO
     v[PCR] = 1.0;
    #endif

    v[TRC] = 1.0;
  }
  else
 {
    v[RHO] = 0.2;
    v[PRS] = 0.02;

    #if CR_FLUID != NO
      v[PCR] = 0.1;
    #endif

    v[TRC] = 0.0;

  }
  v[VX1] = 0.0;

  // g_gamma = 5./3.; g_gammacr = 1.39;
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
FILE *consrv;

   if (g_stepNumber == 0){
     consrv   = fopen("consrv.txt", "w");   
     fprintf(consrv,"# [1]Time, [2]BOXMass, [3]MASStest,\
                       [4]TE_term, [5]KE_term, [6]CR_term, [7]ENGYtest\n\n" );
   } else if (g_stepNumber > 0){
     consrv   = fopen("consrv.txt","a"); 
   }

   //  ----------- CHECKING CONSERVATION ------------------------ //
   int i, j, k;
   double *x1, *dx1, *dV1, dvolume, dmass;
   double BOXMass;
   double TE_term, KE_term, CR_term;

   x1 = grid[IDIR].x; dx1 = grid[IDIR].dx;
   dV1 = grid[IDIR].dV;

   BOXMass = 0.0; 
   dvolume = 0.0; dmass   = 0.0;
   TE_term = 0.0; KE_term = 0.0; CR_term = 0.0;
   DOM_LOOP(k,j,i){

     dvolume = dV1[i];
     dmass = (d->Vc[RHO][k][j][i])*dvolume;

     BOXMass  += dmass;
     TE_term  += ( d->Vc[PRS][k][j][i]/(g_gamma-1.0) )*dvolume;
     KE_term  += 0.5*dmass*pow( d->Vc[VX1][k][j][i], 2.0);
     
     #if CR_FLUID != NO
     CR_term += ( d->Vc[PCR][k][j][i]/(g_gammacr-1.0) )*dvolume;
     #endif

   }

   static double MASS0, TE0, KE0, CR0;
   if (g_stepNumber == 0){
     MASS0 = BOXMass;
     TE0 = TE_term;
     KE0 = KE_term;
     CR0 = CR_term;
   }
   double MASStest =  (BOXMass)/(MASS0);
   double ENGYtest = ( TE_term + KE_term + CR_term)/(TE0 + KE0 + CR0) ;


 // "# [1]Time(sec), [2]BOXMass, [3]MASStest, [4]TE_term, [5]KE_term, [6]CR_term, [7]ENGYtest"

  fprintf(consrv,"%0.7e\t%0.7e %0.9e\t%0.7e %0.7e %0.7e %0.9e\n\n",\
       g_time, MASS0, MASStest, TE_term, KE_term, CR_term, ENGYtest); 
   

   fclose(consrv); 
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *********************************************************************** */
{ 
}
