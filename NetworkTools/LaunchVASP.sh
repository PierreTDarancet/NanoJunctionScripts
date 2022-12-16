#!/bin/bash
echo '# *********************************************************#'
echo '# *  Vasp Job Launcher: Version 0.1                       *#'
echo '# *********************************************************#'
echo '# Copyright (C) 2014 Argonne National Laboratory           #'
echo '#                                                          #'
echo '# This file is distributed under the terms of the          #'
echo '# GNU General Public License. See the file License         #'
echo '# in the root directory of the present distribution,       #'
echo '# or http://www.gnu.org/copyleft/gpl.txt .                 #'
echo '#                                                          #'
echo '# Contributors  : Pierre Darancet  (pdarancet@anl.gov)     #'
echo '# *********************************************************#'
echo '#                                                          #'
echo '# This program generates a pbs scrip for vasp in the       #'
echo '# directory it is launched from and launches it            #'
echo '# ARGUMENTS                                                #'
echo '# $1 = NUMBER OF NODES                                     #'
echo '# $2 = Number of PPN PER node                              #'
echo '# $3 = Walltime requested                                  #'
echo '# $4 = READMEINFO                                          #'
echo '# $5 = NPAR (optional)                                     #'
echo '# *********************************************************#'
NOW=`date "+%Y-%m-%dT-%H%M%S"`
ACCOUNT="staff"
LOCALDIR=`pwd`
NUMBEROFNODES=${1}
NUMBEROFPPN=${2}
TIMEREQUEST="${3}"
READMEINFO="${4}"
echo '# SCRIPT WAS CALLED BY:                                    #'
echo "# ${0} $@"
if [ "$NUMBEROFNODES" = '' ]; then 
	echo "Number of nodes not found/read"
	exit
elif [ "$NUMBEROFPPN" = '' ] ; then
	echo "Number of ppn/node not found/read"
	exit
elif [ "$TIMEREQUEST" = '' ] ; then
	echo "Wall time not found/read"
	exit
elif [ "$READMEINFO" = '' ] ; then
	echo "Readme info not found/read"
	exit
fi
echo '# Nodes:    ' ${NUMBEROFNODES}
echo '# PPN/Node: ' ${NUMBEROFPPN}
echo '# Walltime Requested: ' ${TIMEREQUEST}
echo '# READMEINFO:' "${READMEINFO}" 
echo '# ********************************************************#'
if [ "$5" == "" ] ; then 
	round() { echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc)); };
        tmp=`echo "scale=2; sqrt( $NUMBEROFNODES ) " |bc -l`
        NPAR=`echo $(round $tmp 0)`
        echo '# NPAR not present in call '
	echo "# NPAR set to $NPAR (sqrt(Nodes))"
	echo '# ********************************************************#'
else 
	NPAR=${5}
	echo "# NPAR = ${NPAR} (User provided)" 
	echo '# ********************************************************#'
fi

echo '# ********************************************************#'
echo '# Checking ~/.bashrc for modules... '
module=`cat ~/.bashrc  | grep -v "#" | grep "module" | grep "load" | grep "vasp"`
if [ "$module" = '' ]; then
	echo 'module not loaded in ~/.bashr'
	exit
else
	echo '# The following modules will be loaded:'
	echo "${module}"

fi
echo '# ********************************************************#'
echo '# Checking for vasp INPUTs'

if [ !  -e ./KPOINTS ] ; then
	echo 'KPOINTS file not present'
	exit
elif [ ! -e ./INCAR ] ; then
	echo ' INCAR file not present'
	exit
elif [ ! -e ./POSCAR ] ; then
	echo 'POSCAR file not present'
	exit
elif [ ! -e ./POTCAR ] ; then
	echo 'POTCAR file not present'
	exit
fi
echo '# ... OK'
echo '# ********************************************************#'
echo '# Correcting/Adding NPAR in INCAR'
grep -v "NPAR" ./INCAR  > ./tmp.file 
mv ./tmp.file ./INCAR
echo "NPAR = $NPAR " >> ./INCAR

echo '# ********************************************************#'
echo '# Generating pbscript'
if [ -e ./pbsscript ] ; then 
	echo 'Saving already present pbsscript'
	tmp=`echo ./pbsscript.${NOW} | tr -d ''`
	mv ./pbsscript $tmp
	echo "script moved to $tmp"
fi

echo '# ********************************************************#'

echo "
#!/bin/bash
# This script has been generated on $NOW
# by ${0}
# Launched on `hostname`
# in ${LOCALDIR}
# with the following arguments:
# ${@}
# Modules loaded in ~/.bashrc: ${module}
#  Basics: Number of nodes, processors per node (ppn), and walltime (hhh:mm:ss)
#PBS -l nodes=${NUMBEROFNODES}:ppn=${NUMBEROFPPN}
#PBS -l walltime=${TIMEREQUEST}
#PBS -N `echo ${PWD##*/}`
#PBS -A ${ACCOUNT}
#
#  File names for stdout and stderr.  If not set here, the defaults
# are <JOBNAME>.o<JOBNUM> and <JOBNAME>.e<JOBNUM>
#PBS -o job.out
#PBS -e job.err
#
#  Send mail at begin, end, abort, or never (b, e, a, n). Default is "a".
#PBS -m bea pdarancet@anl.gov
 
# change into the directory where qsub will be executed
cd \$PBS_O_WORKDIR
echo \$PBS_JOBID \$PWD >> \$HOME/jobdirs

# start MPI job over default interconnect; count allocated cores on the fly.
mpirun -machinefile  \$PBS_NODEFILE -np \$PBS_NP \
        vasp" > ./pbsscript

qsub pbsscript

