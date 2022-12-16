#!/bin/bash
MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"
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

echo "Machine: $Machine ..."
echo "  ...script"
rsync -atuvz /home/darancet/scripts/ darancet@$Machine:/home/darancet/scripts
echo "  ...bin"
rsync -atuvz /home/darancet/bin/ darancet@$Machine:/home/darancet/bin/
echo "  ...Done"
echo "  ...INPUT FILES"
rsync -atuvz --max-size=1M /home/darancet/Desktop/DATA/   darancet@$Machine:/home/darancet/DATA
echo "  ...Done"
rsync -atuvz /home/darancet/Desktop/DATA/COMMON/   darancet@$Machine:/home/darancet/DATA/COMMON

