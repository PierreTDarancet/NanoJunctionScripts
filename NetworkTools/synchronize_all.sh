#!/bin/bash
MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"

for Machine in $MACHINELIST
do
echo "Machine: $Machine ..."
echo "  ...script"
rsync -atuvz /home/darancet/scripts/ darancet@$Machine\.apam:/home/darancet/scripts
echo "  ...bin"
rsync -atuvz /home/darancet/bin/ darancet@$Machine\.apam:/home/darancet/bin/
echo "  ...Done"
echo "  ...INPUT FILES"
rsync -atuvz --max-size=1M /home/darancet/Desktop/DATA/   darancet@$Machine\.apam:/home/darancet/DATA
echo "  ...Done"
rsync -atuvz /home/darancet/Desktop/DATA/COMMON/   darancet@$Machine\.apam:/home/darancet/DATA/COMMON
done
