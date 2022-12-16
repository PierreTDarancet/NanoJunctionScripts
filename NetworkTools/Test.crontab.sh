#!/bin/bash

#MACHINELIST="soho tribeca gramercy chinatown les dumbo bronx harlem midtown downtown fidi"
MACHINELIST="bowery bronx chinatown civiccenter downtown dumbo fidi flatiron gramercy harlem hellskitchen les midtown moma noho nolita soho tribeca ues unionsquare"
for Machine in $MACHINELIST
do 
ssh  darancet@$Machine "echo "$Machine" >/tmp/Try.dat && /usr/bin/mail -s "Test" pierre.darancet@gmail.com < /tmp/Try.dat"
done

