#!/bin/bash

MACHINELIST="machinename1 machinename2"

echo '# SCRIPT WAS CALLED BY: #'
echo "${0} $@"
echo 'MACHINE:      ' $1
echo '# ********************************************************#'
echo '# The script is checking its inputs...                    #'
#ADD SUPPORTED OPTION PRINTOUT IN ONESHOT.GENERAL
FOUND=".FALSE."
for MACHINECHECK in  $MACHINELIST
do
if [ "$1" = "$MACHINECHECK" ] ; then
FOUND=".TRUE."
fi
done

if [ "$FOUND" = ".TRUE." ] ; then
echo "Machine input is correct"
Machine="$1"
else
echo "ERROR: $1 is NOT supported... EXITING NOW"
exit 
fi


for Machine in $MACHINELIST
do 
ssh  $Machine "ssh-keygen -t rsa -f ~/.ssh/id_rsa && cat ~/.ssh/id_rsa.pub"
done

