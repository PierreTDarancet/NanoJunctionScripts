#!/bin/bash

MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"




for Machine in $MACHINELIST
do 
ssh $Machine 'rm ~/report ; cd /home/ ; USRLIST=`ls` ; for usr in $USRLIST; do cd /home/$usr ; du -h  > ~/tmp.out ; echo $usr >> ~/report ; tail -n 1 ~/tmp.out >> ~/report ; done ; hostname; cat ~/report' 
done


for Machine in $MACHINELIST
do 
echo "$Machine"
ssh $Machine 'df -h | grep "sda"' 
done
