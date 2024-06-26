# PLUTOv4.2-CRHD 
PLUTO + Cosmic rays in two-fluid approximation.

Copyright (C) Gupta, Sharma, Mignone 2021 MNRAS

------------------------------------------------

This version includes cosmic ray physics under 
two-fluid approximation. Currently, it supports the HD module
in Cartesian and spherical geometries in 1D, 2D, and 3D. 

==================================

Reference

================================== 

Gupta, Sharma, Mignone 2021, MNRAS 
(https://ui.adsabs.harvard.edu/abs/2019arXiv190607200G/). 

==================================

Downloading

==================================

Open terminal and write the following command:
git clone git@github.com:gupta-siddhartha/PLUTOv4.2-CRHD.git <enter>
You can also download it directly in a zip format.

==================================

Setting environment

==================================

Set the path to the PLUTO source directory 
following the PLUTO user's manual, which can be found at Doc/.
e.g., if you are using MAC or Linux system, you can use 
>> export PLUTO_DIR=/Users/siddhartha/PLUTOv4.2-CRHD/

==================================

How to turn on CR-HD module

(or you can jump into the Test Run section)

==================================

CR-HD module is as easy as other modules in PLUTO.
Follow steps 1 - 9.

>>Step 1: python $PLUTO_DIR/setup.py 
          
>>Step 2: Setup problem 
          
>>Step 3: 

PHYSICS HD        and 

USER_DEF_PARAMETERS   1 
          
It will show  

-----------------

EOS                           IDEAL

ENTROPY_SWITCH                NO

CR_FLUID                      NO   <== New

THERMAL_CONDUCTION            NO

VISCOSITY                     NO

CR_DIFFUSION                  NO   <== New

ROTATING_FRAME                NO

-------------

>>Step 4:

EOS must be  IDEAL

>>Step 5: 

In CR_FLUID section, you may find three options:

a. NO
b. NC_PdV_TOTENG
c. NC_DCR_TOTENG

Set CR_FLUID to "NC_PdV_TOTENG", which is reported as "unsplit-pdv (Et+Ecr)" method in Gupta, Sharma, Mignone 2021.
NC_DCR_TOTENG option is NOT recommended, to be removed in future releases.

Once you choose "EOS   IDEAL" and "CR_FLUID  NC_PdV_TOTENG".

>>Step 6: 

It will ask "User-defined Parameters". Here you must name one parameter as "Pcr_Shock".

>>Step 7: 

Other steps such as makefile selection are identical to without CR case. 

>> Step 8: 

Once you click on "Quit", open pluto.ini.
Here you will find "Pcr_Shock            0" at the end.
This set CR pressure fraction "w_{cr}" at shock (see Eq 18 of Gupta, Sharma, Mignone 2021).
If you do not wish to use this parameter leave it as '0'.

>> Step 9: 

Open init.c.
Add the following line
v[PCR] == ...
Adiabatic index of CRs can be set as 
g_gammacr = 4/3. Default is 4/3 which stands for the relativistic gas.

You are all set. Have fun with PLUTO+CR-HD.

==================================

Test run

==================================

We have given a set up of a shock tube test problem
which is discussed in section 5.2 of Gupta, Sharma, Mignone 2021.

To generate output use the lines below.

cd TestProblems/CRHD/ShockTube/

python $PLUTO_DIR/setup.py

<enter, enter,.., Quit>

make

./pluto

This will generate data.0001.tab. To plot CR pressure
on 'gnuplot'
plot 'data.0001.tab u 1:6 w lp
Result should be identical to Figure 6.

Note that, column '6' stands for CRs in 1D [see column specification in tab.out].

==================================

If you need further help, contact www.siddharthagupta.com.

==================================

Acknowledgment

==================================

This code is distributed freely with the hope that
it can be useful to your research problem. 
We will be happy if you find the code useful and appreciate our effort. 
Journal details can be found at https://ui.adsabs.harvard.edu/abs/2019arXiv190607200G/abstract.


Have fun with PLUTO+CR-HD.

==================================


 
