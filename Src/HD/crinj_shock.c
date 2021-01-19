/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shock finding algorithm.
  
  Search and flag computational zones lying in a shock wave.
  The flagging strategy is based on two switches designed to detect 
  the presence of compressive motion or shock waves in the fluid:
  \f[
     \nabla\cdot\vec{v} < 0 \qquad{\rm and}\qquad
     \Delta x\frac{|\nabla p|}{p} > \epsilon_p
  \f]  
  where \f$\epsilon_p\f$ sets the shock strength.
  At the discrete level we replace the two conditions by 
  \f[
    \sum_d \frac{ A_{\vec{i}+\HALF\hvec{e}_d}v_{d,\vec{i}+\HALF\hvec{e}_d}
                 -A_{\vec{i}-\HALF\hvec{e}_d}v_{d,\vec{i}-\HALF\hvec{e}_d} }
                {\Delta{\cal V}_{d,\vec{i}}}  < 0
                \qquad{\rm and}\qquad
    \sum_{d} \left|p_{\vec{i}+\hvec{e}_d} - p_{\vec{i}-\hvec{e}_d}\right|
             < 
     \epsilon_p \min_d\left(p_{\vec{i}+\hvec{e}_d},
                            p_{\vec{i}-\hvec{e}_d},p_{\vec{i}}\right)
  \f]  
  where \f$\hvec{i} = (i,j,k)\f$ is a vector of integer numbers 
  giving the position of a computational zone, while \f$\hvec{e}_d =
  (\delta_{1d},\delta_{2d},\delta_{3d})\f$ is a unit vector in the direction
  given by \c d.
  Once a zone has been tagged as lying in a shock, different flags may be
  switched on or off to control the update strategy in these critical regions.

  This function can be called called when:
  - \c SHOCK_FLATTENING has been set to \c MULTID: in this case shocked zones
    are tagged with \c FLAG_MINMOD and \c FLAG_HLL that will later
    be used to force the reconstruction with the minmod limiter
    and the Riemann solver to HLL.
  - \c ENTROPY_SWITCH has been turned on: this flag will be checked later in
    the ConsToPrim() functions in order to recover pressure from the
    entropy density rather than from the total energy density.
    The update process is:
 
    - start with a vector of primitive variables  <tt> {V,s} </tt>
      where \c s is the entropy;
    - set boundary condition on \c {V}; 
    -  compute \c {s} from \c {V};
    - flag zones where entropy may be used (flag = 1);
    - evolve equations for one time step;
    - convert <tt> {U,S} </tt> to primitive:
      \code
        if (flag == 0) {  // Use energy
          p = p(E)
          s = s(p)
        }else{            // use entropy
          p = p(S)
          E = E(p)
        }  
        \endcode

  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    July 28, 2015

  \Changes, File is renamed as pluto_flag_shock.c
            Added Mach number calculation.
 
  \Changes: Function is renamed and made compatible for cosmic ray injection at shocks
  \Method: Shocked zones are detected using condition i-iii section 4.4
           of Gupta, Sharma, Mignone 2021.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN 0.5
#endif


// Following function are defined by SG
// Start
double GetMachNumber(int k, int j, int i, const Data *d); 
// Caution: Function determines local Mach number not the actual shock Mach number.
// end

//#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH // Removed by SG
/* *************************************************************** */
void DetectShock_CR (const Data *d, Grid *grid)
/*!
 * Detect shocks in presence of CR-fluid
 * Method is identical to flag_shock.c file
 * Some extra filters are introduced for a more robust detection.
 *  
 ***************************************************************** */
{
  int  i, j, k, nv;
  int  ip, jp, kp;
  double divv, gradp, pt_min;
  double dpx1, *dx1, *dV1, ***vx1, pt_min1, dvx1;
  double dpx2, *dx2, *dV2, ***vx2, pt_min2, dvx2;
  double dpx3, *dx3, *dV3, ***vx3, pt_min3, dvx3;
  
  double *dVx, *dVy, *dVz;
  double *Ar, *Ath, *r, *th, s;
  static double ***pt;

  static double ***tmp,***rho;   // NEW
  double gradrho_gradtmp;        // NEW
  double drhox1, drhox2, drhox3; // NEW
  double dtmpx1, dtmpx2, dtmpx3; // NEW

  double *x1, *x2, *x3; // NEW

/* ----------------------------------------------------
   0. Define pointers to variables and allocate memory
   ---------------------------------------------------- */
  x1 = grid[IDIR].x; x2 = grid[JDIR].x; x3 = grid[KDIR].x;

  dx1 = grid[IDIR].dx; dV1 = grid[IDIR].dV;
  dx2 = grid[JDIR].dx; dV2 = grid[JDIR].dV; 
  dx3 = grid[KDIR].dx; dV3 = grid[KDIR].dV; 
  
  Ar  = grid[IDIR].A;
  r   = grid[IDIR].x;
  Ath = grid[JDIR].A;
  th  = grid[JDIR].x;

  EXPAND(vx1 = d->Vc[VX1];  ,
         vx2 = d->Vc[VX2];  ,
         vx3 = d->Vc[VX3];)

  if (pt == NULL)   pt = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double);
  if (rho == NULL) rho = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double); // NEW
  if (tmp == NULL) tmp = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, double); // NEW

  // By siddhartha
  #if CR_FLUID != NO
  double Etot, e_cr, w_cr = g_inputParam[Pcr_Shock];
  e_cr = 2. * w_cr / (1. + w_cr);
  #endif
    
/* --------------------------------------------------------
   1. Compute total pressure and flag, initially all zones
      with ENTROPY_SWITCH.
   -------------------------------------------------------- */

  TOT_LOOP(k,j,i){  
#if EOS == ISOTHERMAL
    pt[k][j][i] = d->Vc[RHO][k][j][i]*g_isoSoundSpeed*g_isoSoundSpeed;
#else
  #if HAVE_ENERGY 
     pt[k][j][i] = d->Vc[PRS][k][j][i];
  #endif    
#endif
  
   rho[k][j][i] = d->Vc[RHO][k][j][i];
   tmp[k][j][i] = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
   
   #if CR_FLUID != NO
    pt[k][j][i] += d->Vc[PCR][k][j][i];
   #endif   

   // Following variables are used to 
   // check different shock detection filters
   // See e.g., section 4.4 in Gupta et al 2021.
   // (Initialized with zero.) 
   d->UDA[0][k][j][i] = 0.0; 
   d->UDA[1][k][j][i] = 0.0; 
   d->UDA[2][k][j][i] = 0.0; 
   d->UDA[3][k][j][i] = 0.0; 

  }
 
/* ----------------------------------------------
   2. Track zones lying in a shock
   ---------------------------------------------- */

  for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++){ 
  for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){ 
  for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){

  /* -- Compute divergence of velocity -- */
     
    #if GEOMETRY == CARTESIAN

     D_EXPAND(dvx1 = vx1[k][j][i + 1] - vx1[k][j][i - 1];   ,
              dvx2 = vx2[k][j + 1][i] - vx2[k][j - 1][i];   ,
              dvx3 = vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == CYLINDRICAL

     D_EXPAND(dvx1 =   Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                     - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);   ,
              dvx2 = vx2[k][j + 1][i] - vx2[k][j - 1][i];           , 
              dvx3 = vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == POLAR

     D_EXPAND(dvx1 =  Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                    - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);  ,
              dvx2 = (vx2[k][j + 1][i] - vx2[k][j - 1][i])/r[i];      ,
              dvx3 =  vx3[k + 1][j][i] - vx3[k - 1][j][i];)

    #elif GEOMETRY == SPHERICAL

     D_EXPAND(dvx1 =  Ar[i]  *(vx1[k][j][i + 1] + vx1[k][j][i])
                    - Ar[i-1]*(vx1[k][j][i - 1] + vx1[k][j][i]);               ,
              dvx2 = (  Ath[j] *(vx2[k][j + 1][i] + vx2[k][j][i])
                     - Ath[j-1]*(vx2[k][j - 1][i] + vx2[k][j][i]))/fabs(r[i]); ,
              s    = th[j];
              dvx3 = (vx3[k + 1][j][i] - vx3[k - 1][j][i])/(r[i]*sin(s));)

    #endif

    divv = D_EXPAND(dvx1/dV1[i], + dvx2/dV2[j], + dvx3/dV3[k]);
 
    if (divv < 0.0){

    d->UDA[0][k][j][i] = 1.0; // Set yes to first filter
    /* -----------------------------------------------
        Compute undivided difference of the total 
        pressure and minimum value in neighbour zones
       ----------------------------------------------- */
       
      pt_min = pt[k][j][i];
      D_EXPAND(pt_min1 = MIN(pt[k][j][i+1], pt[k][j][i-1]); ,
               pt_min2 = MIN(pt[k][j+1][i], pt[k][j-1][i]);  ,
               pt_min3 = MIN(pt[k+1][j][i], pt[k-1][j][i]); )

      D_EXPAND(pt_min = MIN(pt_min, pt_min1);  ,
               pt_min = MIN(pt_min, pt_min2);  ,
               pt_min = MIN(pt_min, pt_min3);)
      
      D_EXPAND(dpx1 = fabs(pt[k][j][i+1] - pt[k][j][i-1]);  ,  
               dpx2 = fabs(pt[k][j+1][i] - pt[k][j-1][i]);  , 
               dpx3 = fabs(pt[k+1][j][i] - pt[k-1][j][i]);)   
                
      gradp = D_EXPAND(dpx1, + dpx2, + dpx3);

       if (gradp > EPS_PSHOCK_FLATTEN*pt_min) {

         d->UDA[1][k][j][i]=1.0;  // Set yes to second filter
         
         d->flag[k][j][i]   |= FLAG_HLL;
        
         d->flag[k][j][i]   |= FLAG_HLL;
         D_EXPAND(
           d->flag[k][j][i+1] |= FLAG_HLL;
           d->flag[k][j][i-1] |= FLAG_HLL;  ,
           d->flag[k][j-1][i] |= FLAG_HLL;  
           d->flag[k][j+1][i] |= FLAG_HLL;  ,
           d->flag[k-1][j][i] |= FLAG_HLL;
           d->flag[k+1][j][i] |= FLAG_HLL;)
           

     /* -----------------------------------------------
         NEW: 
         Compute Grad(rho).Grad(T) ;
       ----------------------------------------------- */       

         D_EXPAND(drhox1 = (rho[k][j][i+1] - rho[k][j][i-1]);  ,
                  drhox2 = (rho[k][j+1][i] - rho[k][j-1][i]);  ,
                  drhox3 = (rho[k+1][j][i] - rho[k-1][j][i]);)

         D_EXPAND(dtmpx1 = (tmp[k][j][i+1] - tmp[k][j][i-1]);  ,
                  dtmpx2 = (tmp[k][j+1][i] - tmp[k][j-1][i]);  ,
                  dtmpx3 = (tmp[k+1][j][i] - tmp[k-1][j][i]);)

       #if GEOMETRY == CARTESIAN 
       gradrho_gradtmp = D_EXPAND(drhox1*dtmpx1/(dx1[i]*dx1[i]),\
                                + drhox2*dtmpx2/(dx2[j]*dx2[j]),\
                                + drhox3*dtmpx3/(dx3[k]*dx3[k]));
       #elif GEOMETRY == SPHERICAL
       gradrho_gradtmp = D_EXPAND(drhox1*dtmpx1/(dx1[i]*dx1[i]),\
                                + drhox2*dtmpx2/(r[i]*r[i])/(dx2[j]*dx2[j]),\
                                + drhox3*dtmpx3/(r[i]*r[i]*sin(th[j])*sin(th[j]))/(dx3[k]*dx3[k]));
       #endif

       d->UDA[2][k][j][i] = GetMachNumber(k, j, i, d);
       D_EXPAND(\
             d->UDA[2][k][j][i-1] = GetMachNumber(k, j, i-1, d); 
             d->UDA[2][k][j][i+1] = GetMachNumber(k, j, i+1, d); ,\
             d->UDA[2][k][j-1][i] = GetMachNumber(k, j-1, i, d); 
             d->UDA[2][k][j+1][i] = GetMachNumber(k, j+1, i, d); ,\
             d->UDA[2][k-1][j][i] = GetMachNumber(k-1, j, i, d); 
             d->UDA[2][k+1][j][i] = GetMachNumber(k+1, j, i, d); );

       
      if(gradrho_gradtmp > 0.0){

       d->UDA[3][k][j][i] = GetMachNumber(k, j, i, d);
       D_EXPAND(\
             d->UDA[3][k][j][i-1] = GetMachNumber(k, j, i-1, d);
             d->UDA[3][k][j][i+1] = GetMachNumber(k, j, i+1, d); ,\
             d->UDA[3][k][j-1][i] = GetMachNumber(k, j-1, i, d); 
             d->UDA[3][k][j+1][i] = GetMachNumber(k, j+1, i, d); ,\
             d->UDA[3][k-1][j][i] = GetMachNumber(k-1, j, i, d); 
             d->UDA[3][k+1][j][i] = GetMachNumber(k+1, j, i, d); );

         }
      }
    } 
      
  }}}
  

#if CR_FLUID != NO  
 //====================================================//
 //                 CR injection                       //
 //====================================================//
  DOM_LOOP(k,j,i){ 
    if(d->UDA[3][k][j][i] >= 1.){
      if(w_cr != 0.0){
        Etot = d->Vc[PRS][k][j][i]/(g_gamma-1.) + d->Vc[PCR][k][j][i]/(g_gammacr-1.);
        d->Vc[PRS][k][j][i] = (1.- e_cr) * Etot * (g_gamma   - 1.);
        d->Vc[PCR][k][j][i] = e_cr       * Etot * (g_gammacr - 1.);
      }
    }
  }
#endif


#ifdef PARALLEL
   AL_Exchange (d->flag[0][0], SZ_char);
#endif
}

#undef  EPS_PSHOCK_FLATTEN 

double GetMachNumber(int k, int j, int i, const Data *d)
/* ****************************************************
 *  Input : [1] grid number i, j, k
 *          [2] Data structure to access primitive variables
 *
 *  Output : [1] returns Mach number.
 *
 *  Method : 
 *           Step 1. Finds the minimum sound speed amoung 
 *                   neighbouring cells (a_min). 
 *           Step 2. Speedof the flaged zone (v).
 *           Step 3. M = v/a_min
 *
 * ****************************************************/
{
  double Mth;
  
  if(i < IOFFSET || i > NX1_TOT-IOFFSET) {return 0;}
  if(j < JOFFSET || j > NX2_TOT-JOFFSET) {return 0;}
  if(k < KOFFSET || k > NX3_TOT-KOFFSET) {return 0;}

  double vx1, vx2, vx3;
  EXPAND(vx1 = d->Vc[VX1][k][j][i];, 
         vx2 = d->Vc[VX2][k][j][i];, 
         vx3 = d->Vc[VX3][k][j][i];)

  double ***rho; rho = d->Vc[RHO];
  double ***prs; prs = d->Vc[PRS];
        
 // Get the minimum sound speed in neighbouring zones
  double a_min;
  double a2, a2_ip1, a2_im1, a2_jp1, a2_jm1, a2_kp1, a2_km1;
 
  a2 = g_gamma*prs[k][j][i]/rho[k][j][i]; // Reference cell
 
  EXPAND( a2_ip1 = g_gamma*prs[k][j][i+1]/rho[k][j][i+1]; 
          a2_im1 = g_gamma*prs[k][j][i-1]/rho[k][j][i-1];,
          a2_jp1 = g_gamma*prs[k][j+1][i]/rho[k][j+1][i];
          a2_jm1 = g_gamma*prs[k][j-1][i]/rho[k][j-1][i];,
          a2_kp1 = g_gamma*prs[k+1][j][i]/rho[k+1][j][i]; 
          a2_km1 = g_gamma*prs[k-1][j][i]/rho[k-1][j][i];  )
 
  EXPAND( a2 = MIN(a2, a2_ip1); a2 = MIN(a2, a2_im1); ,
          a2 = MIN(a2, a2_jp1); a2 = MIN(a2, a2_jm1); ,
          a2 = MIN(a2, a2_kp1); a2 = MIN(a2, a2_km1); )
         
  a_min = sqrt(a2);
  
  // Velocity of the flagged zone
  double v2_shock;
  v2_shock =  EXPAND(vx1*vx1, + vx2*vx2, + vx3*vx3);

  //Mach number;
  Mth = sqrt(v2_shock)/a_min;         
 
  return(Mth);
}
