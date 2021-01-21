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

  /* 
  DOM_LOOP(k,j,i){ // To check the shocked-zone
    Test1[k][j][i]  = d->UDA[3][k][j][i];
  }*/

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





