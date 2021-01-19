#include "pluto.h"
/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  

  double ***Test0;  Test0 = GetUserVar("test0");
  double ***Test1;  Test1 = GetUserVar("test1");

  double ***p, ***rho;
    p   = d->Vc[PRS];
    rho = d->Vc[RHO];
  
  DOM_LOOP(k,j,i){ // To check the shocked-zone
    Test0[k][j][i]  = d->flag[k][j][i];
    Test1[k][j][i]  = d->UDA[3][k][j][i];
  }
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





