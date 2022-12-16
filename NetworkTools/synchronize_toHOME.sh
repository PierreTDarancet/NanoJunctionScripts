#!/bin/bash

MACHINELIST="broadway"
for Machine in $MACHINELIST
do
echo "Machine: $Machine ..."
echo "  ...FILES"
rsync -auvz  /home/darancet/DATA/* darancet@$Machine:/home/darancet/Desktop/DATA/
echo "  ...Done"
done
