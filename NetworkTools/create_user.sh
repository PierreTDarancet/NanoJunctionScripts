echo creating user $1
if  [ -a /usr/local/share/public_keys/id_dsa.pub_$1 ] ; then 
  echo creating user $1
  sudo adduser --gecos emptyname --disabled-password --quiet $1
  # for now give everyone sudo... may revisit this later...
  sudo adduser $1 sudo
  sudo mkdir /home/$1/.ssh
  sudo chown $1 /home/$1/.ssh/
  sudo chgrp $1 /home/$1/.ssh/
  sudo cp /usr/local/share/public_keys/id_dsa.pub_$1 /home/$1/.ssh/authorized_keys2
  sudo chown $1 /home/$1/.ssh/authorized_keys2
  sudo chgrp $1 /home/$1/.ssh/authorized_keys2
  sudo chmod 600 /home/$1/.ssh/authorized_keys2
  sudo echo source /usr/local/share/group_bashrc >> /home/$1/.bashrc
else 
  echo key for user $1 not present in /usr/local/share/public_keys, user will not be created...
fi

