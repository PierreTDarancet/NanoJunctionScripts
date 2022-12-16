#!/bin/bash 

MYNAME="GofR" 
MYVERSION="0.1"
MACHINE=`hostname`
WORKDIR=`pwd`
if [ -e /etc/profile ]; then
	source /etc/profile
fi
if [ -e /etc/profile ]; then
	source ~/.profile
fi
SCRIPTBEGINTIME=`date`

##  \ENV #####
# Copyright (C) 2016 Argonne National Lab
#
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# Contributors  : Pierre Darancet
# ***********************************************
#
# This program compute the distance between two atoms
#   of a cell with periodic boundary conditions
echo '# *********************************************************#'
echo '# *  COMPUTE G(R): Version 0.1                            *#'
echo '# *********************************************************#'
echo '# Copyright (C) 2016 Argonne National Lab                  #'
echo '#                                                          #'
echo '# This file is distributed under the terms of the          #'
echo '# GNU General Public License. See the file License         #'
echo '# in the root directory of the present distribution,       #'
echo '# or http://www.gnu.org/copyleft/gpl.txt .                 #'
echo '#                                                          #'
echo '# Contributors  : Pierre Darancet                          #'
echo '# *********************************************************#'
echo '#                                                          #'
echo '# # This program compute the distance between two atoms    #'
echo '#   of a cell with periodic boundary conditions            #'
echo '# ARGUMENTS                                                #'
NumberofInputFiles="$1"
echo '# $1 = Numbers of files:' "$NumberofInputFiles"
TotalNumberofAtoms="$2"
echo '# $2 = Total Number of Atoms:'  "$TotalNumberofAtoms"
echo '# 1st Atom Information                                     #'
Element_A="${3}"
echo '# $3 = 1st Atom Element Number:' "$Element_A"
echo '# 2nd Atom Information                                     #'
Element_B="${4}"
echo '# $4 = 2nd Atom Element Number:' "$Element_B"
OutputFile="${5}"
echo '# $5 =  Outputfile:' "$OutputFile"
InputFile_1="${6}"
echo '# $6 =  Input file 1:' "$InputFile_1"
InputFile_2="${7}"
echo '# $7 =  Input file 2:' "$InputFile_2"  
echo '# *********************************************************#'
echo '# This shell is running in the following directory:        #'
echo "$WORKDIR"

##################################################
# EXECUTION
##################################################
#echo "$MYNAME $MYVERSION beginning at" `date`  >> $WORKDIR/GofR.README
echo "$MYNAME $MYVERSION beginning at" `date` 
##################################################
FOUNDISSUE=.FALSE.
if [  "$NumberofInputFiles" == ""  ]; then
FOUNDISSUE=.TRUE.
elif [ "$TotalNumberofAtoms" == ""  ]; then
FOUNDISSUE=.TRUE.
elif [ "$Element_A" == ""  ]; then
FOUNDISSUE=.TRUE.
elif [ "$Element_B" == ""  ]; then
FOUNDISSUE=.TRUE.
elif [ "$OutputFile" == ""  ]; then
FOUNDISSUE=.TRUE.
elif [ "$InputFile_1" == ""  ]; then
FOUNDISSUE=.TRUE.
fi
echo $FOUNDISSUE

if [[ $FOUNDISSUE == .TRUE. ]]; then
echo "Problem with arguments: $@"
exit
fi


#for inputnum in $(seq 1 $NumberofInputFiles); do
#  nametmp=`echo 

echo '#Distance between elements:' $Element_A ' and ' $Element_B > $OutputFile
echo '#Format: Distance in (a,b,c) (frac coord), in x,y,z (Ang) Distance (Ang)' >> $OutputFile

# enter cell info

InputFile="$InputFile_1"


A1=`awk -v i=3 -v j=1 'FNR == i {print $j}' $InputFile`
A2=`awk -v i=3 -v j=2 'FNR == i {print $j}' $InputFile`
A3=`awk -v i=3 -v j=3 'FNR == i {print $j}' $InputFile`

B1=`awk -v i=4 -v j=1 'FNR == i {print $j}' $InputFile`
B2=`awk -v i=4 -v j=2 'FNR == i {print $j}' $InputFile`
B3=`awk -v i=4 -v j=3 'FNR == i {print $j}' $InputFile`

C1=`awk -v i=5 -v j=1 'FNR == i {print $j}' $InputFile`
C2=`awk -v i=5 -v j=2 'FNR == i {print $j}' $InputFile`
C3=`awk -v i=5 -v j=3 'FNR == i {print $j}' $InputFile`

echo '# Primitive Cell (Ang)' >> $OutputFile

echo "# ( $A1 ; $A2  ; $A3 )" >> $OutputFile
echo "# ( $B1 ; $B2  ; $B3 )" >> $OutputFile
echo "# ( $C1 ; $C2  ; $C3 )" >> $OutputFile


numberoflines=`cat $InputFile | wc -l`
echo "Total Number of lines $numberoflines" 
# number of iterations computation
tmpnumberA=`echo "scale=0; ( $numberoflines - 7 ) / ( $TotalNumberofAtoms + 1 ) " | bc`
tmpnumberB=`echo "scale=2; ( $numberoflines - 7 ) / ( $TotalNumberofAtoms + 1 ) " | bc`
#test for integer number of iterations

if (( $(bc <<< "$tmpnumberA == $tmpnumberB")  )) 
then 
  IterationNumber=$tmpnumberA
  echo "Number of iterations found=" $tmpnumberA
else  
  echo 'ERROR Check number of atoms' $tmpnumberA $tmpnumberB
  exit
fi


# Transfer the interesting part of input to tmpfile
tmpnumberA=`echo "scale=0; ( $numberoflines - 7 ) " | bc`
tail -n $tmpnumberA $InputFile > tmpfileA.POSCAR

# Printout only the two atoms of interest
# account for first line
tmpnumberA=`echo "scale=0; ( $TotalNumberofAtoms + 1 ) " | bc`
tmpnumberB=`echo "scale=0; ( $Element_A + 1 ) " | bc`
echo $tmpnumberB
if [ "$tmpnumberB" -eq "$tmpnumberA" ]; then
    tmpnumberB=0;
fi

awk -v TotAtom=${tmpnumberA} -v NAtom=${tmpnumberB} 'NR % TotAtom == NAtom' tmpfileA.POSCAR > ListElA.FracCoord
echo $tmpnumberB
tmpnumberB=`echo "scale=0; ( $Element_B + 1 ) " | bc`
echo $tmpnumberB
if [ "$tmpnumberB" -eq "$tmpnumberA" ]; then
    tmpnumberB=0;
fi
#awk -v TotAtom=${tmpnumberA} -v NAtom=0 'NR % TotAtom == NAtom' tmpfileA.POSCAR 
echo ${tmpnumberA} ${tmpnumberB}
awk -v TotAtom=${tmpnumberA} -v NAtom=${tmpnumberB} 'NR % TotAtom == NAtom' tmpfileA.POSCAR > ListElB.FracCoord
paste ListElA.FracCoord ListElB.FracCoord | awk '{printf "%4.8f %4.8f %4.8f \n" , ( ( ( $1 - $4 + 10.5) % 1.0 ) - 0.5 )  , ( ( ( $2 - $5 + 10.5) % 1.0 ) - 0.5 )  ,   ( ( ( $3 - $6 + 10.5) % 1.0 ) - 0.5 )  }' | awk -v aa=$A1 -v ab=$A2 -v ac=$A3 -v ba=$B1 -v bb=$B2 -v bc=$B3 -v ca=$C1 -v cb=$C2 -v cc=$C3 '{printf "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f \n", $1,$2,$3,  ((aa*$1)+(ba*$2)+(ca*$3)), ((ab*$1)+(bb*$2)+(cb*$3)), ((ac*$1)+(bc*$2)+(cc*$3)),  sqrt(   (((aa*$1)+(ba*$2)+(ca*$3))**2) + (((ab*$1)+(bb*$2)+(cb*$3))**2) +  (((ac*$1)+(bc*$2)+(cc*$3))**2)  ) }' >> $OutputFile

#if [[ $NumberofInputFiles <= 1 ]]; 
#then 
#   echo 'done'
#   exit
#else
   
   InputFile="$InputFile_2"
   
A1=`awk -v i=3 -v j=1 'FNR == i {print $j}' $InputFile`
A2=`awk -v i=3 -v j=2 'FNR == i {print $j}' $InputFile`
A3=`awk -v i=3 -v j=3 'FNR == i {print $j}' $InputFile`

B1=`awk -v i=4 -v j=1 'FNR == i {print $j}' $InputFile`
B2=`awk -v i=4 -v j=2 'FNR == i {print $j}' $InputFile`
B3=`awk -v i=4 -v j=3 'FNR == i {print $j}' $InputFile`

C1=`awk -v i=5 -v j=1 'FNR == i {print $j}' $InputFile`
C2=`awk -v i=5 -v j=2 'FNR == i {print $j}' $InputFile`
C3=`awk -v i=5 -v j=3 'FNR == i {print $j}' $InputFile`
echo '# Second file' >> $OutputFile
echo '# Primitive Cell (Ang)' >> $OutputFile

echo "# ( $A1 ; $A2  ; $A3 )" >> $OutputFile
echo "# ( $B1 ; $B2  ; $B3 )" >> $OutputFile
echo "# ( $C1 ; $C2  ; $C3 )" >> $OutputFile


numberoflines=`cat $InputFile | wc -l`
echo "Total Number of lines $numberoflines" 
# number of iterations computation
tmpnumberA=`echo "scale=0; ( $numberoflines - 7 ) / ( $TotalNumberofAtoms + 1 ) " | bc`
tmpnumberB=`echo "scale=2; ( $numberoflines - 7 ) / ( $TotalNumberofAtoms + 1 ) " | bc`
#test for integer number of iterations

if (( $(bc <<< "$tmpnumberA == $tmpnumberB")  )) 
then 
  IterationNumber=$tmpnumberA
  echo "Number of iterations found=" $tmpnumberA
else  
  echo 'ERROR Check number of atoms' $tmpnumberA $tmpnumberB
  exit
fi


# Transfer the interesting part of input to tmpfile
tmpnumberA=`echo "scale=0; ( $numberoflines - 7 ) " | bc`
tail -n $tmpnumberA $InputFile > tmpfileA.POSCAR

# Printout only the two atoms of interest
# account for first line
tmpnumberA=`echo "scale=0; ( $TotalNumberofAtoms + 1 ) " | bc`
tmpnumberB=`echo "scale=0; ( $Element_A + 1 ) " | bc`

echo $tmpnumberB
if [ "$tmpnumberB" -eq "$tmpnumberA" ]; then
    tmpnumberB=0;
fi
awk -v TotAtom=${tmpnumberA} -v NAtom=${tmpnumberB} 'NR % TotAtom == NAtom' tmpfileA.POSCAR > ListElA.FracCoord

tmpnumberB=`echo "scale=0; ( $Element_B + 1 ) " | bc`
echo $tmpnumberB
if [ "$tmpnumberB" -eq "$tmpnumberA" ]; then
    tmpnumberB=0;
fi
awk -v TotAtom=${tmpnumberA} -v NAtom=${tmpnumberB} 'NR % TotAtom == NAtom' tmpfileA.POSCAR > ListElB.FracCoord

paste ListElA.FracCoord ListElB.FracCoord | awk '{printf "%4.8f %4.8f %4.8f \n" , ( ( ( $1 - $4 + 10.5) % 1.0 ) - 0.5 )  , ( ( ( $2 - $5 + 10.5) % 1.0 ) - 0.5 )  ,   ( ( ( $3 - $6 + 10.5) % 1.0 ) - 0.5 )  }' | awk -v aa=$A1 -v ab=$A2 -v ac=$A3 -v ba=$B1 -v bb=$B2 -v bc=$B3 -v ca=$C1 -v cb=$C2 -v cc=$C3 '{printf "%4.8f %4.8f %4.8f %4.8f %4.8f %4.8f %4.8f \n", $1,$2,$3,  ((aa*$1)+(ba*$2)+(ca*$3)), ((ab*$1)+(bb*$2)+(cb*$3)), ((ac*$1)+(bc*$2)+(cc*$3)),  sqrt(   (((aa*$1)+(ba*$2)+(ca*$3))**2) + (((ab*$1)+(bb*$2)+(cb*$3))**2) +  (((ac*$1)+(bc*$2)+(cc*$3))**2)  ) }' >> $OutputFile


#fi



