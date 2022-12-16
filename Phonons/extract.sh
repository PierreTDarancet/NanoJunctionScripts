#!/bin/bash
echo "Use: $0 Natoms Nconfigurations"
Natoms=$1 
Nconfigurations=$2

ii=0
grep "%block AtomicCoordinatesAndAtomicSpecies" -A $Natoms  minus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $1, $2, $3)}' > tmp1

grep "Atomic forces" -A $Natoms  minus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $2, $3, $4)}' > tmp2

echo $ii > forces_minus
paste tmp1 tmp2 >> forces_minus



grep "%block AtomicCoordinatesAndAtomicSpecies" -A $Natoms  plus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $1, $2, $3)}' > tmp3

grep "Atomic forces" -A $Natoms  plus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $2, $3, $4)}' > tmp4

echo $ii > forces_plus
paste tmp3 tmp4 >> forces_plus





while [ $ii -lt $Nconfigurations ]; do
ii=`expr $ii + 1`
echo $ii

grep "%block AtomicCoordinatesAndAtomicSpecies" -A $Natoms  minus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $1, $2, $3)}' > tmp1

grep "Atomic forces" -A $Natoms  minus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $3, $4, $5)}' > tmp2

echo $ii >> forces_minus
paste tmp1 tmp2 >> forces_minus



grep "%block AtomicCoordinatesAndAtomicSpecies" -A $Natoms  plus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $1, $2, $3)}' > tmp3

grep "Atomic forces" -A $Natoms  plus.out$ii | tail -n $Natoms \
      | awk '{printf("%10.5f   %10.5f    %10.5f\n", $3, $4, $5)}' > tmp4

echo $ii >> forces_plus
paste tmp3 tmp4 >> forces_plus




done
