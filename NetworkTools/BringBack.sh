#!/bin/bash

MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"


for Machine in $MACHINELIST
do 
ssh  $Machine "chmod a+x /home/darancet/synchronize_toHOME.sh && /home/darancet/synchronize_toHOME.sh"
done

