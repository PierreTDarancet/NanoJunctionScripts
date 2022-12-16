boomerang.py
====================
**SOURCE:** :download:`boomerang.py <//home/cam1/git_projects/software/utils/boomerang.py>`

Boomerang is a simple program to submit runs to a cluster of machines using SSH. 
The code is object oriented and therefore should be easily incorporated into other programs.
Also, it is easy to call from the command line. The program stores information in a series
of data files (described below). The basic essence of the program invloves a few simple steps: 

#. First, queue jobs on your local machine for a given cluster.
#. Second, determine the status of the cluster by finding out how many jobs are running on each. This is accomplished by a universally updated status, which works for all machines. 
#. Submit the jobs to machines that contain available CPU and memory for your jobs.
#. While the jobs are running boomerang has convenient functions to let you view
   the results. 
#. Finally you can retreive any runs that have been finished.

**Disk space**: Effective August 2014, each of the compute nodes will will go through each directory in ``/home/*/tmp_runs/*`` (where boomerang runs) and, if NO file has been touched in 31 days, delete WAVECAR, CHGCAR, and vasprun.xml files. If you need one of these files, just rename it, e.g. ``gzip CHGCAR``, and it will not be deleted. 

Files
---------------

All data files pertaining to this program are stored in a hidden directory in your home directory .boomerang. "Global" files are stored at ~/.boomerang/file. "Local" files are stored at ~/.boomerang/default/file. Some files can be global or local; this means that if there exists a file at ~/.boomerang/default/file then it will be used; else the file at ~/.boomerang/file will be used. ("default" here can be any queue name, e.g. "dualsocket".) We will discuss applications of this in a bit. 

Here are the files that boomerang uses:

  * *environment* (global) - This file contains the environment variable LD_LIBRARY_PATH (and any others) that need to be set on the
    remote machines. The format of this file is one line per machine, with the first word being the machine name and the
    following words being a bash command. The default value will be applied to all machines unless the machine is explicitly
    named. This file is actually not required; the default value is hard coded into the program.

  * *queue* (local) - Shows which directories on your local machine are queued to be submitted.

  * *running* (local) - Shows which files are currently running on the remote machines.

  * *finished* (local) - This contains a list of all the jobs that have finished.

  * *required* (global or local) - A list of files which must be present in order to queue a job.

  * *exclude* (global or local) - Do not transfer these files back from the remote machine. (See "Disk Space" above for how long they are retained on the remote machine.)

  * *machlist* (global or local) - contains a list of machines to submit to. These may be shortcuts defined in your ssh config file or complete IP names. Only one machine per line should be listed. Lines beginning with # will not be read. 

  * *status* (global) - The query method will look up how many jobs there are on each machine and write the result here.

  * *script* (global or local) - This is the default script that will be submitted unless 
    a file name boomerang_submit.sh exists in the
    directory. An example script is shown here. The first step is to copy the
    vasp executable to the local directory with a distinct name so you know
    which executable corresponds to which job. boomerang will replace any
    instance of XXXX in your script with a four-digit random integer. A time
    command is included for convenience, though it may not be necessary.  There
    is a an if statement which executes a postprocessing script if it is
    present in the running directory. Finally, the vasp executable is removed.
    Please note that this is
    only an example and the script can be whatever you want. This is not tied
    to vasp in any particular way.

.. code-block:: bash 

  ## SCRIPT DEFAULT
  #BOOMMEM=10 # Memory expected in GB
  #BOOMCPU=400 # should be 100*ncores
  ncores=4
  npar=2

  # Find executables
  vasploc=`which vasp5-intel-parallel`
  if [ "$vasploc" = '' ]; then vasploc="/usr/local/bin/vasp5-intel-parallel"; fi
  mpiloc=`which mpirun`
  if [ "$mpiloc"  = '' ]; then mpiloc="/opt/intel/impi/4.1.3.048/intel64/bin/mpirun"; fi
  cp $vasploc vaspexe.XXXX

  # Setup and run
  echo "NPAR=$npar" >> INCAR
  ulimit -s unlimited
  /usr/bin/time -a -o out -f "\n%PCPU \n%U user-sec \n%S system-sec \n%E elapsed-hr:min:sec \n%e elapsed-sec" $mpiloc -n $ncores ./vaspexe.XXXX

  if [ -f postp.sh ]; then ./postp.sh ; fi

  rm ./vaspexe.XXXX

If you're using this script, notice a few things:
  
  - The location has a default; boomerang does not source ~/.bashrc, so may have trouble finding the executable in the PATH
  - We change the executable to exe.XXXX, which lets us see e.g. whether the executable is running (boomerang.py -check) or which directory is hogging all the memory (top)
  - NPAR and number of cores (ncores) may depend on the queue (e.g. dual-socket machines), so this just puts it straight into the INCAR

The software repository has a "dotboomerang" directory that contains information for a few different clusters and options; see below how to get started with this. 

  - The default queue runs VASP 5 on single-core (3770 and 3930)
  - The dual queue runs VASP5 on dual-core (mostly 2670v2). This has a modified machlist
  - The lowmem queue runs VASP5 on single-core machines, but is willing to scatter with only 3GB of memory available. This has a modified script. 


Multiple clusters: If you want a different cluster (e.g. high memory), we can make a new directory ~/.boomerang/highmem, which contains all the local files and an adjusted "machlist" file. If you want a different script (e.g. abinit instead of vasp), we can make a directory ~/.boomerang/abinit with all the local files and an adjusted "script" file. We do this as follows:

.. code-block:: bash

  $ mkdir ~/.boomerang/highmem
  $ touch ~/.boomerang/highmem/{queue,finished,running}
  $ vi ~/.boomerang/highmem/machlist

Setup
---------------

To setup boomerang for the first time for VASP, just copy the software repository's dotboomerang/*  into ~/.boomerang run:

.. code-block:: bash
  
  $ mkdir ~/.boomerang
  $ cp dotboomerang/* ~/.boomerang/

Make sure you have boomerang.py in your PATH and inputs.py in your PYTHONPATH. You can adjust these in ~/.bashrc.

Universal status
---------------------

As of May 2014, boomerang can check for a status maintained across all machines in our cluster. This alleviates the need for everyone to query each machine. There is a script that automatically updates the universal status file every few minutes, and every time someone submits something, boomerang updates the status file accordingly. For more details, see the readme file in the boomerang_public directory (currently on columbuso; hopefully soon to be moved to grandcentral). 

If one person uses boomerang without the universal status, there is a greater chance of multiple people submitting to one machine within a short time span, which will make jobs run slower. That said, if you want to do it, here are some ways to do it:

  1. Command-line usage: When using a command that needs querying (e.g. scatter, force, status), add the command "-queryall" to manually check each machine in machlist (alternatively: "-noremote"). If you queried it recently and want to use the saved status in ~/.boomerang/status, use the command "-noquery" instead. Note that your status may be out of date, and there is no way to figure that out without querying. 

  2. Python usage: There is no automatic querying; just instantiate the class as bb=boomerang(noremote=True), and the query method won't look for a universal status.

  3. Change boomerang.py: during __init__, change the default value of noremote to True.


Command-line usage
-------------------------

We will now discuss how to do things from the command line. For the most part, each task is associated with a given
method. The methods themselves are documented in the code and these comments are pulled out into the docstring section
below. The most important methods are: query, queue, scatter, get.

* **Identifying the cluster** - The name keyword may be combined with all of the below commands. If name is not given, it
  will be assumed that the name is default. One may want to have different clusters for different purposes; see above.

  .. code-block:: bash 

    $ boomerang.py name=dual

* **Preparing remote machines** - Create the remote run directories on the of the remote hosts (i.e. currently called ~/tmp_runs). Each run will be placed at ~/tmp_runs/directorynameXXXX, where XXXX is a random 4-digit number. This may have been already accomplished when your profile was set up on the remote machine. 

  .. code-block:: bash 

    $ boomerang.py -prep

* **Querying availability** - Check all of the machines to see how many jobs are running on each. This updates the status file. There are a few ways to do this:

  - "boomerang.py -query" : Attempts to retrieve status from the master (universal) status file; any failure or machines that are missing are updated "by hand"
  - "boomerang.py -query -quick" : Attempts to retrieve status from the master (universal) status file; any machines that are missing are NOT updated "by hand"
  - "boomerang.py -queryall" : Queries all machines manually, without inquiring of universal status
  - "boomerang.py -noquery" : Reads the (local) status file that was last written, although it is very likely out of date. This is not recommended. 

* **Queueing** - Go to the directory which contains the files you want to submit and execute:

  .. code-block:: bash 

    $ boomerang.py -queue

* **Summary** - Show what jobs are in the queue, running.

  .. code-block:: bash 

    $ boomerang.py -sum

* **Status** - Show how many jobs are running on the cluster and by whom

  .. code-block:: bash 

    $ boomerang.py -status

* **Scatter queued jobs** - Send queued jobs to available slots on remote machines. Unless you write -noquery, it will ensure that the usage is current. 

  Other command-line flags exist for scatter: hogmax=XX tells boomerang to take no more than XX% (default 75) of the cluster, to leave space for other users; hogmin=XX tells boomerang to take XX% (default 30) even if nothing is left for others; the script will take some amount in between these, depending on how much is empty. 
  To use all default values, just run:

  .. code-block:: bash
  
    $ boomerang.py -query
    $ boomerang.py -scatter


* **Cat running jobs** - Cat a given file for each running directory on the remote machine. If no file is given, the default is OSZICAR:

  .. code-block:: bash 

    $ boomerang.py -cat INCAR

* **Go to remote machine** - Change to a running directory and execute boomerang with the go option and it will take you to the running directory on the remote machine:

  .. code-block:: bash 

    $ boomerang.py -go

* **Get finished job** - Retreive the results from all finished calculations:

  .. code-block:: bash 

    $ boomerang.py -get

* **Force** - Force a job to start on a given remote machine regardless of its status. There is no need to queue the job in this case, and you should be in the directory of interest upon executing the command:

  .. code-block:: bash 

    $ boomerang.py force=gramercy

* **Kill** - This will kill a job. You can provide the call with a filename, and this file should contain a list of all directories that should be killed. If no filename is given, it will try to kill the current directory:

  .. code-block:: bash 

    $ boomerang.py -kill [filename]

* **Force Get** - This will force boomerang to try and get files from a remote host for the current directory. It is not relevant if they have already been retrieved. Of course, it will fail if they have already been deleted.

  .. code-block:: bash 

    $ boomerang.py -forceget

* **Check** - This will check each directory that is currently running to ensure the executable is still running. It returns the directories for which no executable is found. (It assumes that the executable has the four-number XXXX in the name; see the above example script.)

* **remove** - This will remove the directory from any finished, queue, or running files. 


Docstrings from boomerang.py 
-----------------------------

.. automodule:: boomerang
   :members:

