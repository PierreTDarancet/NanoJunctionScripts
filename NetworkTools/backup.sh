#!/bin/bash
BACKUPDIR="/home/darancet/BACKUP"
MACHINELIST="machine1 machine2"
CURRENTBACKUP=`date "+%Y-%m-%dT%H:%M:%S"`

for Machine in $MACHINELIST
do
echo "Machine: $Machine ..."
echo "  ...Backing up"
rsync -azP --link-dest=$BACKUPDIR/current /home/darancet/Desktop/DATA   $Machine:$BACKUPDIR/back-$CURRENTBACKUP 
echo "  ...Done"
echo "  ...Links"
ssh  $Machine "rm -f $BACKUPDIR/current && ln -s $BACKUPDIR/back-$CURRENTBACKUP $BACKUPDIR/current"
echo "  ...Done"


#echo "  ...Backing up"
#rsync -azP --max-size='500M' --exclude "*" --include "DATA*"   /home/darancet/Desktop/   $Machine:$BACKUPDIR/BACKUP
#echo "  ...Done"

done
