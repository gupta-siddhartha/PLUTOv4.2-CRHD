#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          WENO3
#define  TIME_STEPPING           RK3
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     1

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  CR_FLUID                NC_PdV_TOTENG
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  CR_DIFFUSION            NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  Pcr_Shock               0

/* [Beg] user-defined constants (do not change this line) */

#define  NC_HYBRID               NO
#define  EPS_PSHOCK_FLATTEN      0.1

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    NO
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   NO
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             DEFAULT
