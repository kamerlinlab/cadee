#!/bin/bash

#
# Retrieve temperatures (solvent, solute) and energy contributions (vdw, el, ...) from Qdyn log file
# 


logfile="$(pwd)/$1"
dat_folder=dynplot
if [[ -e "$dat_folder" ]]; then
    echo "Please remove the folder '$dat_folder' manually
"
    exit 1
fi
mkdir $dat_folder
cd $dat_folder

# energy contributions
#
#  Example output:  
# 
## =======================  Energy summary at step   1500 ========================
##                         el       vdW      bond     angle   torsion  improper
## solute                0.00      0.00     13.31      9.56      4.26      0.70
## solvent           -5623.35    806.80      0.00      0.00      0.00      0.00
## solute-solvent        0.00      0.00
## LRF                 -30.79
## Q-atom             -224.52     -4.85    -73.72      5.66      0.00      0.00
## 
##                      total       fix slvnt_rad slvnt_pol     shell    solute
## restraints         -218.28      0.00   -241.91     18.14      0.00      5.49
## 
##                      total potential   kinetic
## SUM               -4261.11  -5335.21   1074.10
#
# Save all this data to:
# E_solute.dat
# E_solvent.dat
# E_solute-solvent.dat
# E_LRF.dat
# EQ.dat
# E_restraints.dat
# E_total.dat

grep "Energy summary at step" -A15 ${logfile} | awk  '{ 

if ($2 == "Energy") step=$6;
else if ($1 == "solute") print step, $2, $3, $4, $5, $6, $7 >> "E_solute.dat";
else if ($1 == "solvent") print step, $2, $3, $4, $5, $6, $7 >> "E_solvent.dat";
else if ($1 == "solute-solvent") print step, $2,$3 >> "E_solute-solvent.dat";
else if ($1 == "LRF") print step, $2 >> "E_LRF.dat";
else if ($1 == "Q-atom") print step, $2, $3, $4, $5, $6, $7 >> "EQ.dat";
else if ($1 == "restraints") print step, $2, $3, $4, $5, $6, $7 >> "E_restraints.dat"; 
else if ($1 == "SUM") print step,$2,$3,$4 >> "E_total.dat"; 
}' 

# Q Energy contributions in each state
#
#
## ======================= Q-atom energies at step    500 ========================
## type   st lambda        el       vdW      bond     angle   torsion  improper
## Q-Q     1 1.0000    -22.88      3.10
## Q-Q     2 0.0000    -70.69     26.38
## 
## Q-wat   1 1.0000   -235.81     -3.69
## Q-wat   2 0.0000    -84.07    -18.18
## 
## Q-surr. 1 1.0000   -235.81     -3.69
## Q-surr. 2 0.0000    -84.07    -18.18
## 
## Q-any   1 1.0000   -258.70     -0.59    -93.59      1.65      0.00      0.00
## Q-any   2 0.0000   -154.76      8.20     -3.42     45.02      0.00      0.00
## 
## type   st lambda     total restraint
## Q-SUM   1 1.0000   -347.10      4.12
## Q-SUM   2 0.0000   -104.50      0.45
#
# Save all this data to:
# EQ_qq_1.dat
# EQ_qq_2.dat
# EQ_qwat_1.dat
# EQ_qwat_2.dat
# EQ_qsurr_1.dat
# EQ_qsurr_2.dat
# EQ_bonded_1.dat
# EQ_bonded_2.dat
# EQ_total_1.dat
# EQ_total_2.dat

grep "Q-atom energies at step" -A17 ${logfile} | awk  '{ 

if ($2 == "Q-atom") step=$6;
else if ($1 == "Q-Q" && $2 == "1") print step, $3, $4, $5 >> "EQ_qq_1.dat";
else if ($1 == "Q-Q" && $2 == "2") print step, $3, $4, $5 >> "EQ_qq_2.dat";
else if ($1 == "Q-wat" && $2 == "1") print step, $3, $4, $5 >> "EQ_qwat_1.dat";
else if ($1 == "Q-wat" && $2 == "2") print step, $3, $4, $5 >> "EQ_qwat_2.dat";
else if ($1 == "Q-surr." && $2 == "1") print step, $3, $4, $5 >> "EQ_qsurr_1.dat";
else if ($1 == "Q-surr." && $2 == "2") print step, $3, $4, $5 >> "EQ_qsurr_2.dat";
else if ($1 == "Q-any" && $2 == "1") print step, $3, $6, $7, $8, $9 >> "EQ_bonded_1.dat";
else if ($1 == "Q-any" && $2 == "2") print step, $3, $6, $7, $8, $9 >> "EQ_bonded_2.dat";
else if ($1 == "Q-SUM" && $2 == "1") print step, $3, $4, $5 >> "EQ_total_1.dat";
else if ($1 == "Q-SUM" && $2 == "2") print step, $3, $4, $5 >> "EQ_total_2.dat";
}' 





# MAKE GNUPLOT INPUTS

echo '
set term wxt persist
set key out

set title "Solute energy contributions"
plot       "E_solute.dat"         u :2 w lp ps 0.5 lw 0.5  title "El",\
           "E_solute.dat"         u :3 w lp ps 0.5 lw 0.5  title "vdW",\
           "E_solute.dat"         u :4 w lp ps 0.5 lw 0.5  title "Bond",\
           "E_solute.dat"         u :5 w lp ps 0.5 lw 0.5  title "Angle",\
           "E_solute.dat"         u :6 w lp ps 0.5 lw 0.5  title "Torsion",\
           "E_solute.dat"         u :7 w lp ps 0.5 lw 0.5  title "Improper"
' > E_solute.plot

echo '
set term wxt persist
set key out

set title "Solvent energy contributions"
plot       "E_solvent.dat"        u :2 w lp ps 0.5 lw 0.5  title "El",\
           "E_solvent.dat"        u :3 w lp ps 0.5 lw 0.5  title "vdW",\
           "E_solvent.dat"        u :4 w lp ps 0.5 lw 0.5  title "Bond",\
           "E_solvent.dat"        u :5 w lp ps 0.5 lw 0.5  title "Angle",\
           "E_solvent.dat"        u :6 w lp ps 0.5 lw 0.5  title "Torsion",\
           "E_solvent.dat"        u :7 w lp ps 0.5 lw 0.5  title "Improper"
' > E_solvent.plot

echo '
set term wxt persist
set key out

set title "Solute-solvent energy contributions"
plot       "E_solute-solvent.dat" u :2 w lp ps 0.5 lw 0.5   title "El",\
           "E_solute-solvent.dat" u :3 w lp ps 0.5 lw 0.5   title "vdW"

' > E_solute-solvent.plot

echo '
set term wxt persist
set key out

set title "LRF energy contribution"
plot "E_LRF.dat"                  u :2 w lp ps 0.5 lw 0.5   title "LRF El"

' > E_LRF.plot

echo '
set term wxt persist
set key out

set title "Q energy contributions (both states)"
plot       "EQ.dat"        u :2 w lp ps 0.5 lw 0.5  title "El",\
           "EQ.dat"        u :3 w lp ps 0.5 lw 0.5  title "vdW",\
           "EQ.dat"        u :4 w lp ps 0.5 lw 0.5  title "Bond",\
           "EQ.dat"        u :5 w lp ps 0.5 lw 0.5  title "Angle",\
           "EQ.dat"        u :6 w lp ps 0.5 lw 0.5  title "Torsion",\
           "EQ.dat"        u :7 w lp ps 0.5 lw 0.5  title "Improper"
' > EQ.plot



echo ' 
set term wxt persist
set key out

set title "Restraints energy contributions"
plot       "E_restraints.dat"     u :2 w lp ps 0.5 lw 0.5   title "Total",\
           "E_restraints.dat"     u :3 w lp ps 0.5 lw 0.5   title "Fix",\
           "E_restraints.dat"     u :4 w lp ps 0.5 lw 0.5   title "Solvent_rad",\
           "E_restraints.dat"     u :5 w lp ps 0.5 lw 0.5   title "Solvent_pol",\
           "E_restraints.dat"     u :6 w lp ps 0.5 lw 0.5   title "Shell",\
           "E_restraints.dat"     u :7 w lp ps 0.5 lw 0.5   title "Solute"
' > E_restraints.plot

echo '
set term wxt persist
set key out

set title "SUM of energies"
plot       "E_total.dat"          u :2 w lp ps 0.5 lw 0.5   title "Total",\
           "E_total.dat"          u :3 w lp ps 0.5 lw 0.5   title "Potential",\
           "E_total.dat"          u :4 w lp ps 0.5 lw 0.5   title "Kinetic"
' > E_total.plot

echo '
set term wxt persist
set key out

set title "Q energies: Nonbonding interactions"
plot "EQ_qq_1.dat"      u :3 w lp ps 0.5 lw 0.5  title "(1) El_qq",\
     "EQ_qq_2.dat"      u :3 w lp ps 0.5 lw 0.5  title "(2) El_qq",\
     "EQ_qq_1.dat"      u :4 w lp ps 0.5 lw 0.5  title "(1) vdW_qq",\
     "EQ_qq_2.dat"      u :4 w lp ps 0.5 lw 0.5  title "(2) vdW_qq",\
     "EQ_qwat_1.dat"    u :3 w lp ps 0.5 lw 0.5  title "(1) El_qwat",\
     "EQ_qwat_2.dat"    u :3 w lp ps 0.5 lw 0.5  title "(2) El_qwat",\
     "EQ_qwat_1.dat"    u :4 w lp ps 0.5 lw 0.5  title "(1) vdW_qwat",\
     "EQ_qwat_2.dat"    u :4 w lp ps 0.5 lw 0.5  title "(2) vdW_qwat",\
     "EQ_qsurr_1.dat"   u :3 w lp ps 0.5 lw 0.5  title "(1) El_qsurr",\
     "EQ_qsurr_2.dat"   u :3 w lp ps 0.5 lw 0.5  title "(2) El_qsurr",\
     "EQ_qsurr_1.dat"   u :4 w lp ps 0.5 lw 0.5  title "(1) vdW_qsurr",\
     "EQ_qsurr_2.dat"   u :4 w lp ps 0.5 lw 0.5  title "(2) vdW_qsurr"
' > EQ_nonbond.plot

echo '
set term wxt persist
set key out

set title "Q energies: bonding" 
plot "EQ_bonded_1.dat"   u :3 w lp ps 0.5 lw 0.5  title "(1) Bond",\
     "EQ_bonded_1.dat"   u :4 w lp ps 0.5 lw 0.5  title "(1) Angle",\
     "EQ_bonded_1.dat"   u :5 w lp ps 0.5 lw 0.5  title "(1) Torsion",\
     "EQ_bonded_1.dat"   u :6 w lp ps 0.5 lw 0.5  title "(1) Improper",\
     "EQ_bonded_2.dat"   u :3 w lp ps 0.5 lw 0.5  title "(2) Bond",\
     "EQ_bonded_2.dat"   u :4 w lp ps 0.5 lw 0.5  title "(2) Angle",\
     "EQ_bonded_2.dat"   u :5 w lp ps 0.5 lw 0.5  title "(2) Torsion",\
     "EQ_bonded_2.dat"   u :6 w lp ps 0.5 lw 0.5  title "(2) Improper"

' > EQ_bonding.plot

echo '
set term wxt persist
set key out

set title "Q energies: Total" 
plot "EQ_total_1.dat"   u :3 w lp ps 0.5 lw 0.5  title "State (1)",\
     "EQ_total_2.dat"   u :3 w lp ps 0.5 lw 0.5  title "State (2)"

' > EQ_total.plot


gnuplot E_solute.plot
gnuplot E_solvent.plot
gnuplot E_solute-solvent.plot
gnuplot E_LRF.plot
gnuplot EQ.plot
gnuplot E_restraints.plot
gnuplot E_total.plot
gnuplot EQ_nonbond.plot
gnuplot EQ_bonding.plot
gnuplot EQ_total.plot
