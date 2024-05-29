#define  PHYSICS                 HD
#define  DIMENSIONS              1
#define  COMPONENTS              1
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     7

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  CR_FLUID                NC_PdV_TOTENG
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  CR_DIFFUSION            NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  E_SN                    0
#define  M_SN                    1
#define  RHO_AMB                 2
#define  TEMP_AMB                3
#define  mu_AMB                  4
#define  Pcr_Shock               5
#define  OPT                     6

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH             CONST_pc
#define  UNIT_VELOCITY           1.e5
#define  UNIT_DENSITY            CONST_mH
#define  unitPRS                 (UNIT_DENSITY*pow(UNIT_VELOCITY,2.0))
#define  unitMASS                (UNIT_DENSITY*pow(UNIT_LENGTH,3.0))
#define  unitTIME                (UNIT_LENGTH/UNIT_VELOCITY)
#define  unitENERGY              (unitMASS*pow(UNIT_LENGTH,2.0)*pow(unitTIME,-2.0))

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    NO
#define  PRINT_TO_FILE       YES
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             DEFAULT
