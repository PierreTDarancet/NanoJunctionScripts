#!/bin/bash
MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"
for Machine in $MACHINELIST
do
echo "Machine: $Machine ..."
echo "  ...script"
rsync -auvz /home/darancet/scripts/ darancet@$Machine:/home/darancet/scripts
echo "  ...bin"
rsync -auvz /home/darancet/bin/ darancet@$Machine:/home/darancet/bin/
echo "  ...COMMON FILES"
rsync -auvz  /home/darancet/Desktop/DATA/COMMON   darancet@$Machine:/home/darancet/DATA/COMMON
echo "  ...Done"
echo "  ...INPUT FILES"
rsync -auvz --max-size=3M /home/darancet/Desktop/DATA   darancet@$Machine:/home/darancet/DATA
echo "  ...Done"
done
