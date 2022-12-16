#!/usr/bin/env python 
import os,sys,glob,random,subprocess,re,pwd,time,math
import signal
try:
  from inputs import inputs
except ImportError:
  from mysub import inputs


########################
# SSH INTERFACE
########################

# Issue: When a naive ssh -oConnectTimeout flag is used, it can hang if the node has little memory, because it's not timing out, but it's slow as molasses. 
# This is not OK for things like query, which need to be done relatively quickly. If I can't get through in 1 second, just forget it and try the next node.
# Resolution: Have Python track the actual time it takes
# Much of this is based on http://stackoverflow.com/questions/24921527/option-for-ssh-to-timeout-after-a-short-time-clientalive-connecttimeout-dont and http://stackoverflow.com/questions/1191374/subprocess-with-timeout/1191537#1191537
class TimeoutException(Exception):   # Custom exception class
  pass
def TimeoutHandler(signum, frame):   # Custom signal handler
  raise TimeoutException
def get_process_children(pid):
    p = subprocess.Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
              stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    return [int(p) for p in stdout.split()]

def ssh1(mach,remotecommands,sshargs=[],timeout=30):
  """ 
  SSH into mach, with ssh arguments sshargs, and execute remotecommands. If the entire process takes longer than timeout, the process and its children are killed.
  Timeout here is the python-goverened timeout, and has nothing to do with when SSH ConnectTimeout
  """
  # Change the behavior of SIGALRM
  OriginalHandler = signal.signal(signal.SIGALRM,TimeoutHandler)
  command=['ssh']+sshargs+[mach,remotecommands]
  # Start the timer. Once 30 seconds are over, a SIGALRM signal is sent.
  signal.alarm(timeout)
  # This try/except loop ensures that you'll catch TimeoutException when it's sent.
  p = subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  try:
    out,err= p.communicate()
    signal.alarm(0) # reset alarm
  except TimeoutException:
    err= "SSH command timed out."
    out=''
    pids=[p.pid]
    pids.extend(get_process_children(p.pid))
    for pid in pids:
      # process might have died before getting to this line so wrap to avoid OSError: no such process
      try: 
        os.kill(pid, signal.SIGKILL)
      except OSError:
        pass
  # Reset the alarm stuff.
  signal.signal(signal.SIGALRM,OriginalHandler)
  return out,err


############################
### Main Boomerang script
############################

class boomerang:
  """This class is a simple queueing system which relies only on SSH. To use:\n
  (1) Choose a script / machine list\n
  (2a) Queue directory(ies) and scatter, or\n
  (2b) Send (force) a directory to a chosen machine\n
  (3) Get the runs back\n
  """
  def __init__(self,name,noremote=False,debug=False):
    """Initialize the key variables of the class. Only the name of the cluster is required. """
    if not name:name='default'
    self.name=name
    self.debug=debug
    self.noremote=noremote
    self.localnodes=[]
    self.username = pwd.getpwuid(os.getuid())[0]
    if self.username == 'localname': self.username = 'remotename' # If user has a different username on the cluster, change the strings here appropraitely
    self.nodes={} # structure: {mach:[usedCPU,freeCPU,freeMEM,totCPU,userdic]}  userdic={username:usedCPU}
    self.users={} # this has to be updated each time nodes is updated
    path=os.environ["HOME"]; self.path=path
    # Information where universal status is stored
    self.ustatus={'mach':'columbuso',
        'dir':'/home/mckornbluth/boomerang_public'}
    #self.ustatus={'mach':'grandcentral',
    #    'dir':'/home/boomeranger/boomerang_public'}
    self.file={}
    # first define files which may be global or local...
    for ff in "script exclude required environment machlist".split():
      tfileg="%s/.boomerang/%s"%(path,ff)
      tfilel="%s/.boomerang/%s/%s"%(path,name,ff)
      if os.path.exists(tfileg): self.file[ff]=tfileg
      if os.path.exists(tfilel): self.file[ff]=tfilel
      if not self.file.keys().count(ff):self.file[ff]=''
    # now save specific files, and overwrite global if present...
    for ff in "running finished queue".split():
      self.file[ff]="%s/.boomerang/%s/%s"%(path,name,ff)
    #global files only
    for ff in 'status'.split():
      self.file[ff]="%s/.boomerang/%s"%(path,ff)

    # check for required files...
    for i in 'machlist'.split():
      if not os.path.exists(self.file[i]):sys.stderr.write("Cannot run without %s\n"%i);sys.exit()
          
    # create the machine list...
    self.machl=[s for s in open(self.file["machlist"]).read().split() if s.strip()[0]!='#']
    if not self.machl:print "No machines given in machine list %s"%self.file["machl"];sys.exit()

    # make a set from required files...
    if os.path.exists(self.file['required']):
      self.required=set(open(self.file["required"]).read().strip().replace("\n"," ").split())

    # get the environment variables...
    self.environ={}
    self.environ['default']="source /opt/intel/icsxe/2013.1.046/ictvars.sh > /dev/null"
    if os.path.exists(self.file['environment']):
      for i in open(self.file['environment']).readlines():
        self.environ[i.split()[0]]=" ".join(i.split()[1:])

    # parse script for requirements info
    self.parse_script()

  def prep(self):
    """Create temporary run directory called tmp_runs on each remote machine in the cluster. """
    print "Prepping remote boomerang tmp_runs directory..."
    for i in self.machl: os.system("ssh %s 'mkdir tmp_runs'"%i)

  def kill(self,inp=''):
    """Kill all or some selected set of jobs. If inp is empty it is assumed that the job in the 
       current directory will be killed. Otherwise, inp can be a filename which contains a list
       of directories to be killd.
       Note that a killed directory will not be retrieved with the "get" method, but "force_get" should still work"""
    # ssh br ps -A | grep 7419
    if inp:
      if os.path.exists(inp):inp=[s.strip() for s in open(inp).readlines() if s.split()]
      else:print "file %s does not exist"%inp;sys.exit()
    else:inp=['./']
    for i in inp:
      os.chdir(i.strip())
      temp=glob.glob("running_boomerang_*")
      if not temp:raise ValueError( "no run to kill...")
      else: temp=temp[0]
      mach=temp.split("_")[-2]
      ran=temp.split("_")[-1]
      pid=os.popen("ssh %s ps -A | grep %s"%(mach,ran)).readline().split()
      if not pid:sys.stderr.write('Job in directory %s not running...\n'%i);continue
      else:pid=pid[0]
      os.system("ssh %s 'kill -9 %s'"%(mach,pid))
      print "job %s was killed..."%i.strip()
      os.system("touch killed_job")
      # Remove from running and remove running file; add to finished
      with open(self.file['running']) as f:
        contents=[s for s in f.readlines() if not s.strip()==os.path.abspath(i)+'/']
      with open(self.file['running'],'w') as f:
        f.write(''.join(contents))
      os.system("rm %s/running_boomerang_*"%i)
      with open(self.file['finished'],'a') as f:
        f.write(i+'\n')
      
      #self.njobs[mach]-=1 # everything is simpler if it's just updated by querying 

  def parse_script(self,script='',defaultcpu=400,defaultmem=13):
    """
    cpu is given in 100*ncores; memory is given in GB
    Parse script for information about run requirements
    If no script argument is given, it uses self.file['script']
    The script should include the following lines:
    #BOOMMEM=8
    #BOOMCPU=200
    (spaces don't matter)
    """
    self.jobcpu=defaultcpu
    self.jobmem=defaultmem # in GB
    if not script:
      with open(self.file['script']) as f:
        script=f.read()
    script=script.strip().splitlines()
    for line in script:
      dat=line.replace(' ','')
      if dat[:9]=='#BOOMMEM=': self.jobmem=float(dat.split('#')[1][8:])
      if dat[:9]=='#BOOMCPU=': self.jobcpu=float(dat.split('#')[1][8:])
    # TODO: detect from mpirun / mpiexec call 
    if self.jobcpu<=0: raise ValueError("CPU must be positive")
    if self.jobmem<=0: raise ValueError("MEM must be positive")
      
  ###############################################################
  # Methods to return information about current runs
  # cat: cat file for each running job
  # go: ssh to the machine & directory
  # print_summary: prints queue, running
  # check: Ensure executable hasn't failed
  ###############################################################

  def cat(self,file=''):
    """Cat the results of the "file" for each running job. If no file is given, the file OSZICAR will be catted."""
    if not file:file='OSZICAR'
    INP=[ s.strip() for s in open(self.file["running"]) if s.strip()]
    for i in INP:
      temp=glob.glob("%s/running_boomerang_*"%i)
      if temp:
        t1,t2=temp[0].split("/")[-1].split("_")[-2:]
        t3=i.split("/")[-2]
        print i,t1
        print os.popen("ssh -o ConnectTimeout=10  %s 'cat tmp_runs/%s%s/%s'"%(t1,t3,t2,file)).read()
        print
      else: print "\nDirectory %s is listed as running but running_boomerang_* is not present...\n"%i

  def go(self,ddir=''):
    """Go to the machine running the job in directory dir or in the current directory if dir is empty.  This will work even if the job has finished and been transferred back as long as the remote directory has not been removed. """
    file='bsubmit.sh' if not ddir else ddir+'/'+'bsubmit.sh'
    if not os.path.exists(file):sys.stderr.write("File %s does not exist...\n"%file);sys.exit()
    t1,t2=open(file).readline().split()[1:]
    os.system("ssh -t %s 'cd tmp_runs/%s;/bin/bash -i'"%(t1,t2))

  def print_summary(self):
    """Print out a summary of the queued files, the running files, but not status of the cluster (use print_status() for that). """
    print "="*15+' QUEUE '+"="*15
    print open(self.file["queue"]).read().strip()
    print "="*15+' RUNNING '+"="*15
    print open(self.file["running"]).read().strip()

  def check(self,goto=False,timeout=6):
    """Checks all running directories to ensure that it is running. This works only if the script includes creation of an executable with XXXX in its name.
    """
    running = [s.strip() for s in open(self.file["running"]) if s.strip()]
    for i in running:
      if not os.path.exists(i):sys.stderr.write("Directory %s does not exist...\n"%i); continue
      temp=glob.glob("%s/running_boomerang_*"%i)
      if not temp:sys.stderr.write('running_boomerang file missing in %s\n'%i);continue
      t1,t2=temp[0].split("/")[-1].split("_")[-2:]
      t3=i.split("/")[-2]
      t4=temp[0].split("/")[-1]
      dcontent=os.popen("ssh -o ConnectTimeout=%s  %s 'ls tmp_runs/%s%s'"%(timeout,t1,t3,t2)).read()
      # t1=mach, t2=rand, t3=dir, t4 = running_boomerang_etc
      if dcontent.strip() and dcontent.count(t4): 
        out=os.popen("ssh -o ConnectTimeout=%s  %s 'ps -e | grep %s'"%(timeout,t1,t2)).read()
        if not out.strip():
          print 'no executable found for %s'%i
          if goto:
            self.go(i)


  #################################################
  # Methods to get and save availability information
  # To get information:
  #   query_one: Gets info from one node. Saves and returns info.
  #   query_all: Gets info from all nodes in machlist. Saves info. 
  #   parse_status: Transforms file contents into memory.
  #   fetch_remote: Gets info from remote. Saves info. 
  #   fetch_local: Gets info from remote. Saves info. 
  #   query: Try remote; if not, local; fill in missing nodes
  # To display information:
  #   print_status: Prints status currently saved
  # To save information:
  #   request_update: Requests an update on remote 
  #   (un)lock_remote: (Un)locks remote file
  #   format_status: Transforms memory into file contents. Return string.
  #   write_local: Writes status to local. 
  #   write_remote: Writes status to remote

  def query_one(self,mach,timeout=2,minjp=60,minmem=15,quiet=False):
    """Query one node. Returns either 1 (for error), or the self.nodes entry:
      [usedCPU,freeCPU,freeMEM,totCPU,userdic]
      where userdic={username:usedCPU}
      Also writes to self.nodes.
      Arguments: machine name (as in .ssh/config), timeout for connection, min cpu to count as a job, min memory to count as a job, and quiet mode
    """
    #command=["ssh","-o ConnectTimeout=%s"%timeout,mach,"nproc; free -m; ps -eouser:16,pcpu,pmem"]
    #out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    out,err = ssh1(mach,"nproc; free -m; ps -eouser:16,pcpu,pmem",[],timeout)
    if err:
      if not quiet: sys.stderr.write("*** %s has an error\n"%mach)
      if mach in self.nodes: del self.nodes[mach]
      return 1
    else:
      outdat=out.strip().splitlines()
      totcpu=100*int(outdat[0])
      dat = [ [f(t) for f,t in zip((str,float,float),s.split()) ] for s in outdat[6:]]
      # Lines 0-5 are nproc, mem info, and column titles
      pcpu= sum([s[1] for s in dat if s[0]!='root'])
      # For free memory, we want to ignore the cache, bc it can be freed up easily.
      # See here: http://askubuntu.com/questions/155768/how-do-i-clean-or-disable-the-memory-cache
      # This is as simple as changing outdat[2] -> outdat[3]
      usage=[ pcpu, totcpu-pcpu, float(outdat[3].split()[3])*1e-3, totcpu, {} ]
      users=set( [self.username]+[s[0] for s in dat if (s[1]>minjp or s[2]>minmem) and s[0]!='root' ])
      for us in users:
        usage[4][us]= sum([s[1] for s in dat if s[0]==us])
      self.nodes[mach]=usage
      if not quiet:
        sys.stderr.write("%16s: Free: %4i%% CPU, %4.1f GB. Users: "%(mach,int(self.nodes[mach][1]),self.nodes[mach][2]))
        for cpu,usname in sorted([(t,s) for s,t in self.nodes[mach][4].items()],reverse=True):
          if cpu<50.: continue
          sys.stderr.write("%s (%d%%) "%(usname,cpu))
        sys.stderr.write('\n')
      return usage

  def query_all(self,timeout=2,minjp=60,minmem=15,quiet=False):
    """Find out how many jobs are running on each node of the cluster and who is running them.  """
    self.nodes={}
    self.users={}
    # query all the machines to see how many jobs are running on each...
    for m in self.machl:
      usage=self.query_one(m,timeout,minjp,minmem,quiet)
    self.calc_users()
    if not quiet: self.print_status(nomachs=True) # users only

  def parse_status(self,stringorlist):
    """ Method to read a status file. Info is loaded into memory
    Input is either the filestring or the filestring.strip().splitlines()
    File structure is space-delimited. Each line is 
    machname usedCPU freeCPU freeMEM totCPU user1 user1CPU user2 user2CPU ...
    Returns list of machines successfully parsed
    """
    self.users={}
    machs=[]
    if type(stringorlist)==type(''): 
      alldat=stringorlist.strip().splitlines()
    else:
      alldat=stringorlist
    for line in alldat:
      if line.strip()[0]=='#': continue # don't even keep it; prob from old version of boomerang
      dat=line.strip().split()
      if len(dat)<=2: continue # prob from old version of boomerang
      if not dat[0] in self.machl: continue # not relevant for this machinelist
      # position 0 is machname. [1,4] is machine info. i,i+1 is user,CPU for i in range(5,len(dat),2)
      self.nodes[dat[0]]= [float(s) for s in dat[1:5]] + [ dict( (dat[i],float(dat[i+1])) for i in range(5,len(dat),2) )]
      machs.append(dat[0])
      if not self.username in self.nodes[dat[0]][4]: self.nodes[dat[0]][4][self.username]=0.
    return machs

  def calc_users(self):
    """ Uses self.nodes to get information about users, to be stored in self.users={username:totCPU}, summed only over nodes in the machlist. Returns the dictionary self.users
    """
    self.users={self.username:0.}
    for key,val in self.nodes.items():
      if not key in self.machl: continue
      for usname,cpu in val[4].items():
        if not usname in self.users: self.users[usname]=0.
        self.users[usname]+=cpu
    return self.users

  def fetch_remote(self,timeout=10,ntries=3,updatewait=10,lockwait=10):
    """ Fetch and load info from remote status, as defined in self.ustatus['mach'] and self.ustatus['dir']. Files in that directory are status, lock, and update. 
    Returns list of machines successfully retrieved
    Default wait time is 10 seconds... after 3 tries, this is still shorter than the about 45 seconds to query all
    """
    command=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"if [ -f %s/update ]; then echo 1; elif [ -f %s/lock ]; then echo 2; else cat %s/status; fi"%tuple([self.ustatus['dir']]*3)]
    fetched=False
    for itry in range(ntries):
      out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
      if err:  #failed to fetch
        sys.stderr.write("Attempt %s/%s: Error in connecting to remote.\n"%(itry+1,ntries))
        continue
      alldat=out.strip().splitlines()
      if alldat==['1']: # needs update
        sys.stderr.write("Attempt %s/%s: Remote update needed. Waiting %s seconds\n"%(itry+1,ntries, updatewait))
        # Wait updatewait seconds and try again, unless it's the last try
        if itry+1<ntries: time.sleep(updatewait)
        continue
      elif alldat==['2']: # locked
        sys.stderr.write("Attempt %s/%s: File locked. Waiting %s seconds\n"%(itry+1,ntries, lockwait))
        # Wait lockwait seconds and try again, unless it's the last try
        if itry+1<ntries: time.sleep(lockwait)
        continue
      # Iff we're here, it looks like we got it
      machs = self.parse_status(alldat[:]) # formerly this was [1:], which was a bug
      sys.stderr.write("Successfully fetched status from remote source.\n")
      fetched=True; self.remote=True;self.changed=False
      return machs
    if not fetched: raise IOError("Failed to fetch status from remote source.")

  def fetch_local(self):
    """Fetch and load file from local status. Does NOT update this status."""
    with open(self.file['status']) as f:
      self.parse_status(f.read())
    sys.stderr.write("Successfully fetched status from local source.\n")

  def query(self,remotetimeout=10,onetimeout=1,minjp=60,minmem=15,quiet=False,quick=False):
    """Tries to fetch status remotely; if cannot, queries machines. If only gets some machines and quick==True, queries the rest. """
    # Clear previous status
    self.nodes={} 
    self.users={}
    try:
      if self.noremote: raise IOError # go straight to local
      # fetch remotely
      machs = self.fetch_remote(remotetimeout)
      for m in self.machl:
        if m in machs or quick: continue # already loaded or don't want to load
        sys.stderr.write("Querying %s... "%m)
        usage=self.query_one(m,onetimeout,minjp,minmem,True)
        self.localnodes.append(m)
        sys.stderr.write("Done\n")
      if not quiet: self.print_status()
    except IOError: # couldn't get it remotely
      # query locally
      self.query_all(onetimeout,minjp,minmem,quiet)

  def print_status(self,allmachs=False,nomachs=False):
    """Prints information about usage of each node and each user. Argument can request all machines, even those not in self.machl"""
    if not nomachs:
      for mach in sorted(self.nodes):
        if not (allmachs or mach in self.machl): continue 
        sys.stderr.write("%16s: Free: %4i%% CPU, %4.1f GB. Users: "%(mach,int(self.nodes[mach][1]),self.nodes[mach][2]))
        for cpu,usname in sorted([(t,s) for s,t in self.nodes[mach][4].items()],reverse=True):
          if cpu<50.: continue
          sys.stderr.write("%s (%d%%) "%(usname,cpu))
        sys.stderr.write('\n')
    self.calc_users()
    for user in sorted(self.users):
      sys.stderr.write("%16s: Using %6.1f cores across cluster\n"%(user,self.users[user]/100.))
    totuse=sum([self.nodes[m][0] for m in self.nodes])
    tothav=sum([self.nodes[m][3] for m in self.nodes])
    sys.stderr.write('%16s: Using %4d out of %4d cores across cluster (%2d%%)\n'%tuple(['**total**']+[int(round(x)) for x in [totuse/100.,tothav/100.,100.*totuse/tothav if tothav!=0 else 0]]) )

  def request_update(self,timeout=10):
    command=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"touch %s/update"%self.ustatus['dir']]
    out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if err: raise IOError("Failed to request update")

  def lock_remote(self,timeout=10):
    command=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"touch %s/lock"%self.ustatus['dir']]
    out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if err: raise IOError("Failed to lock remote")

  def unlock_remote(self,timeout=10,ignoredone=False):
    command=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"rm -f %s/lock"%self.ustatus['dir']]
    out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if err: 
      if 'No such file or directory' in err and ignoredone: return
      command=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"touch %s/unlock"%self.ustatus['dir']]
      out,err=subprocess.Popen(command,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
      raise IOError("Failed to unlock remote. Error message: ",err)

  def format_status(self,prevfile='',remote=False):
    """ Method to format a status file for writing. Output is the string to go to file
    If prevfile exists (should be the file's contents, a string or splitted string), any machine that does not exist in self.machlist is retained in the output file. 
    File structure is space-delimited. Each line is 
    machname usedCPU freeCPU freeMEM totCPU user1 user1CPU user2 user2CPU ...
    If remote is enabled, it  does not include usages of zero
    (formerly, "remote" also prohibited writing nodes that were not read from remote.)
    """
    out=''
    for key,val in sorted(self.nodes.items()):
      #if key in self.localnodes and remote: continue # don't write it because it is not being updated automatically.  # On second thought, it's ok because that node will be deleted when the updater can't get in touch with it. 
      out+='%16s %5d %5d %3d %5d '%tuple([key]+val[:4])
      for us,cpu in sorted(val[4].items()):
        if remote and cpu<50: continue
        out+='%16s %5d '%(us,cpu)
      out+='\n'
    if type(prevfile)==type(''): prevfile=[s for s in prevfile.splitlines() if s.strip()]
    for line in prevfile:
      if line.strip()[0]=='#': continue # don't write; prob from old version of boomerang
      dat=line.split()
      if len(dat)<=2: continue # prob from old version of boomerang
      if not (dat[0] in self.nodes): 
        out += line+'\n'
    return out

  def write_local(self):
    """Writes status to local file (usually ~/.boomerang/status) """
    with open(self.file['status']) as f:
      prevfile=f.read()
    with open(self.file['status'],'w') as f:
      f.write(self.format_status(prevfile))
    sys.stderr.write("Successfully wrote local status.\n")

  def write_remote(self,timeout=10,ntries=3,overrideremote=False,selfunlock=True,lockwait=5):
    """ Writes status to remote location. 
    Checks if someone else has locked it
    overrideremote= write to remote even though nothing was read from remote.
    selfunlock = write even if locked, provided it was locked by same username. 
    """
    # Check it is sensible to write to remote
    if not dir(self).count('remote') or not self.remote:
      if not overrideremote: 
        raise ValueError("You are trying to write remote without having read from remote. Override option is necessary.")
    # Generate file and check
    readcommand=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],"cat %s/status"%self.ustatus['dir']]
    out,err=subprocess.Popen(readcommand,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    prevfile=out
    # Generate file
    out = self.format_status(prevfile,remote=True)
    if "'" in out: raise NotImplementedError("Single-quote in file; need to clean contents for ssh")
    cmdcheck=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],r"stat -c %U "+'%s/lock 2>&1'%(self.ustatus['dir']) ]
    cmdwrite=["ssh","-o ConnectTimeout=%s"%timeout,self.ustatus['mach'],'''echo '%s' > %s/status'''%(out,self.ustatus['dir'])]
    written=False
    for itry in range(ntries):
      # Ensure it's not already locked by someone else:
      fetched=False
      out,err=subprocess.Popen(cmdcheck,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
      if err:  #failed 
        sys.stderr.write("Attempt %s/%s: Error in checking for lock.\n"%(itry+1,ntries))
        continue
      if not 'stat: cannot stat' in out: #locked
        # If it is locked by someone else, he/she is writing to it. This could be a problem: Say Alice scattered 1 job to bronx, then Bob scattered 1 job to brooklyn. When Alice writes the queue, bronx will be correct and bkln will be wrong. When Bob writes the queue, he will overwrite bronx to the original. So it's easiest to just request a full update. 
        # Catch: If there are multiple instances of boomerang running by one user, they could interfere. 
        if (not selfunlock) or out.strip()!=self.username:
          sys.stderr.write("Attempt %s/%s: File locked by %s. Aborting write; requesting update instead.\n"%(itry+1,ntries,out.strip()))
          self.request_update()
          break
      # Iff we're here, it looks like we're safe
      try:
        self.lock_remote()
        out,err=subprocess.Popen(cmdwrite,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        if err: 
          sys.stderr.write("Failed to write remote")
          raise IOError("Failed to write remote")
        written=True
        sys.stderr.write("Successfully wrote remote status.\n")
        break # exit loop
      except:
        self.request_update()
      finally:
        self.unlock_remote()
    if not written: raise IOError("Failed to write to remote source.")

  ###############################################################
  # Methods to handle submission and retrieval
  # submit_run: Submits one run
  # force: submits a run to a machine
  # force_get: retrieves files. 
  # get: Get all finished runs
  # add_to_queue
  # find_resources: Figure out which machines are available 
  # scatter
  # remove: Removes directory from queue, running, and finished
  ###############################################################

  def submit_run(self,ddir,dest):
    """This function handles all the tasks with submitting a job. This includes creating a 
       finalized submission script, transfering files, and submit the script. Please note
       that a final line is added to the submission script which removes the file running_boomerang_X_X.
       This is the only way that boomerang will know that a job has finished. """
    ran=random.randint(1000,9999)
    # Choose libraries:
    # If there is a machine name, choose that. Elif there is a queue name, choose that. Else, default
    if dest in self.environ:
      libs=self.environ[dest]
    else:
      libs=self.environ['default']
    if not ddir[-1]=='/': ddir = ddir+'/'
    # make sure directory is not already running... redundancy in two conditions...
    if glob.glob("%s/running_boomerang_*"%ddir): 
      sys.stderr.write("Cannot submit %s; running_boomerang already exists...\n"%ddir);return
    if [re.sub("/*$","/",s.strip()) for s in open(self.file['running']).readlines()].count(ddir):
      sys.stderr.write("Cannot submit %s... Already in running file\n"%ddir);return
    print "transferring %s to %s"%(ddir,dest)
    if os.path.exists("%s/boomerang_submit.sh"%ddir):script="%s/boomerang_submit.sh"%ddir
    elif os.path.exists(self.file["script"]):script=self.file["script"]
    else:sys.stderr.write("Cannot find boomerang submission script...");sys.exit()
    self.changed=True
    self.nodes[dest][0]+=self.jobcpu
    self.nodes[dest][1]-=self.jobcpu
    self.nodes[dest][2]-=self.jobmem
    self.nodes[dest][4][self.username]+=self.jobcpu
    os.system("touch %s/running_boomerang_%s_%s"%(ddir,dest,ran))
    OUT=open("%s/bsubmit.sh"%ddir,'w')
    OUT.write("# %s %s%s\n"%(dest,ddir.split("/")[-2],ran))
    OUT.write(open(script).read().replace("XXXX",str(ran)))
    OUT.write("rm ~/tmp_runs/%s%s/running_boomerang_%s_%s\n"%(ddir.split("/")[-2],ran,dest,ran))
    OUT.close()
    flag=0
    try:
      os.system("chmod u+x %s/bsubmit.sh"%ddir)
      os.system("scp -qpr %s %s:tmp_runs/%s%s"%(ddir,dest,ddir.split("/")[-2],ran))
      os.system("ssh %s 'cd tmp_runs/%s%s;%s;pwd>>~/.running;nohup /bin/bash ./bsubmit.sh 2>err.bsubmit > out.bsubmit &'"%(dest,ddir.split("/")[-2],ran,libs))
      OUT=open(self.file["running"],'a');OUT.write(ddir+'\n');OUT.close()
    except: flag=1
    return flag

  def force(self,inp):
    """Force the job in the current directory to be submitted to machine inp. """
    if not inp:print "must give destination machine to force";sys.exit() 
    ddir=os.getcwd()
    self.submit_run(ddir,inp)

  def force_get(self,ddir=''):
    """Force boomerang to get the files from the current run directory. It is not relevant
    when this run was performed. An attempt will be made to transfer the files from 
    the remote node. Of course, this will fail if they have already been deleted."""
    file='bsubmit.sh' if not ddir else ddir+'/'+'bsubmit.sh'
    if not os.path.exists(file):sys.stderr.write("File %s does not exist... cannot force_get\n"%file);sys.exit()
    t1,t2=open(file).readline().split()[1:]
    #os.system("rsync --exclude-from=%s -a %s 'cd tmp_runs/%s;/bin/bash -i'"%(self.file['exclude'],t1,t2))
    sys.stdout.write("rsync --exclude-from=%s -a %s:tmp_runs/%s/ %s/ \n"%(self.file['exclude'],t1,t2,'./' if not ddir else ddir))
    try: os.system("rsync --exclude-from=%s -a %s:tmp_runs/%s/ %s/ \n"%(self.file['exclude'],t1,t2,'./' if not ddir else ddir))
    except:sys.stderr.write("Problem forcing transfer of files for %s"%ddir)

  def get(self,timeout=6):
    """Get all running files that have finished. This function checks to see that running_boomerang_X_X
       no longer exists on the remote and if not it will copy all files back."""
    nodone=[]
    done=[]
    INP=[ s.strip() for s in open(self.file["running"]) if s.strip()]
    for i in INP:
      if not os.path.exists(i):
        sys.stderr.write("Directory %s does not exist...\n"%i)
        nodone.append(i)
        continue
      temp=glob.glob("%s/running_boomerang_*"%i)
      if not temp:
        sys.stderr.write('running_boomerang file not present in %s\nCannot get it...\n'%i)
        nodone.append(i)
        continue
      t1,t2=temp[0].split("/")[-1].split("_")[-2:]
      t3=i.split("/")[-2]
      t4=temp[0].split("/")[-1]
      dcontent=os.popen("ssh -o ConnectTimeout=%s  %s 'ls tmp_runs/%s%s'"%(timeout,t1,t3,t2)).read()
      if dcontent.strip() and not dcontent.count(t4): 
        print "from %s to   %s"%(t1.ljust(16),i)
        os.system("rsync --exclude-from=%s -a %s:tmp_runs/%s%s/ %s"%(self.file['exclude'],t1,t3,t2,i))
        os.system("rm %s/running_boomerang_*"%i)
        done.append(i)
        #if self.njobs.keys().count(t1):self.njobs[t1]-=1 #if it completed already, we'll be subtracting too much!
      else:nodone.append(i)
    OUT=open(self.file["running"],'w');OUT.write("\n".join(nodone)+'\n');OUT.close()
    if done: OUT=open(self.file["finished"],'a');OUT.write("\n".join(done)+"\n");OUT.close()
  
  def add_to_queue(self,direc=''):
    """Add a directory to the queue. """
    if not direc:direc=os.getcwd()
    direc=re.sub("/*$","/",direc)
    for d in sorted(glob.glob(direc)):
      # do not queue for many reasons... many are redundant...
      if os.path.exists("%s/queue_boomerang"%d):sys.stderr.write("Already queued %s...\n"%d);continue
      if os.path.exists("%s/bsubmit.sh"%d):sys.stderr.write("Cannot queue %s... delete bsubmit.sh\n"%d);continue
      if glob.glob("%s/running_boomerang_*"%d):sys.stderr.write("Cannot queue %s... run in progress\n"%d);continue
      if [re.sub("/*$","/",s.strip()) for s in open(self.file['queue']).readlines()].count(d):
        sys.stderr.write("Cannot queue %s... Already queued\n"%d);continue
      if [re.sub("/*$","/",s.strip()) for s in open(self.file['running']).readlines()].count(d):
        sys.stderr.write("Cannot queue %s... Already running\n"%d);continue
      # do not queue if required files are not present...
      if self.required-set([ s.split("/")[-1] for s in  glob.glob("%s/*"%d)]):
        sys.stderr.write("Required input files not present for %s\n"%d)
        continue
      sys.stderr.write("Adding %s to queue.\n"%d)
      open("%s/queue_boomerang"%d,'w').close()
      OUT=open(self.file["queue"],'a');OUT.write(d+'\n');OUT.close()
# check queue file before adding...

  def find_resources(self,hogmin=25,hogmax=75,verbose=False):
    """ Figure out which machines are available for jobs. Returns list of machines to scatter to (including some repeated elements if appropriate to scatter >1 job to it)
    Hog specifies the amount to hog the cluster in percent: If there is hogmin% (25 default) or less available, I will hog it all up. On the other extreme, if the cluster is entirely empty, I will only take hogmax%% (75 default). The middle is interpolated.\n
    This won't work right if your login username isn't the username on the cluster. So you should hardcode self.username above, e.g. add the line \n
      if self.username == 'localname': self.username = 'remotename'
    \n
    Algorithm for equitable distribution: 
    If myfrac=.75 (hogmax), I leave .25 (1-hogmax) free. If myfrac<=.3 (hogmin), I leave 0 free. Interpolate: If I have myfrac, I leave freefrac=(1-hogmax)/(hogmax-hogmin) * (myfrac-hogmin) [ but if freefrac<0, set freefrac=0]
    Let mi, mf = myfrac initial (before) and final (after scattering). Let hp,hm equal hogmax,hogmin. Let a = frac available, l = frac I leave, t = frac I take. .
    Then the three equations we need are:
    1) t = mf - mi
    2) a = t + l 
    3) l = (mf-hm)*(1-hp)/(hp-hm)
    We want to eliminate l and mf, and solve for t, so we get
    t = hm - mi - (hp-hm)*(1-a-mi)/(1-hm)
    or
    takefrac = hogmax - myfrac - (hogmax-hogmin)*(1-availfrac-myfrac)/(1-hogmin)
    With the provision that 0 <= t <= a
    """
    if not self.nodes.keys():
      sys.stderr.write("Cluster status file does not exist; running query\n")
      self.query()
    if self.users=={}:
      self.calc_users()
    # only submit to machines that query has found online and are part of the machl... 
    machl=[s for s in self.nodes if s in self.machl]
    if not machl: raise ValueError("No machines found")
    # determine available resources...
    # How much have I used?
    total=sum([self.nodes[m][3] for m in machl]) 
    myfrac=self.users[self.username]/float(total)
    # How much CPU is available? It is available iff there is both necessary memory and necessary CPU to run the job. Allowed to overload CPU by up to 50. 
    availables=[]
    # OLD: # Sort by CPU usage; on average this will be random with a tendency to the unused machines. (Will fill one machine fully before starting an empty machine, but hopefully that is no big deal.) 
    for cpuused,mach in sorted([ (self.nodes[m][0],m) for m in machl]):
      maxjobsmem=int(self.nodes[mach][2]/self.jobmem)
      maxjobscpu=int((self.nodes[mach][1]+50)/self.jobcpu) 
      availables+=[mach]*min(maxjobsmem,maxjobscpu)
    availfrac=len(availables)*self.jobcpu/total
    sys.stderr.write("%s processes are available (%d%%; job: %s%% CPU, %s GB). You are taking %s%% of the cluster. "%(len(availables),100*availfrac,self.jobcpu,self.jobmem,int(round(myfrac*100))))
    if verbose: sys.stderr.write("Machines: "+', '.join(availables)+'. ')
    # See algorithm above for why these are the numbers
    hogmin = .01 * hogmin; hogmax = .01 * hogmax
    assert 0<=hogmin<=1 and 0<=hogmax<=1
    takefrac = hogmax - myfrac - (hogmax-hogmin)*(1-availfrac-myfrac)/(1-hogmin)
    if takefrac<=0:
      takers=[]
    elif 0<takefrac<availfrac:
      takers = availables[ : int(takefrac * total / self.jobcpu) ]
    elif takefrac>=availfrac:
      takers=availables
    sys.stderr.write('Scatter takes up to %s. \n'%(len(takers)))
    # Shuffle it so that it's not alphabetical, which could be an issue in terms of (A) disk space, and (B) if the remote isn't updated, bigger chance of submitting to the same machine.
    #   In principle, we could've programmed it to take the least-free node first. But that will have more chance of two users, scattering at the same time, submitting to the same machine. 
    #   People frequently launch with a script that runs every X minutes, so it's plausible two people will be launching at the same time. 
    #   If it's random, we get to avoid this--it's less likely two people will submit to the same machine. 
    random.shuffle(takers)
    return takers

  def scatter(self,hogmin=25,hogmax=75):
    """Scatter the queued jobs out to all the avaiable nodes on the cluster.\n
    Hog specifies the amount to hog the cluster in percent: See find_resources for details\n
    This won't work right if your login username isn't the username on the cluster. So you should hardcode self.username above, e.g. add the line \n
      if self.username == 'localname': self.username = 'remotename'
    """
    libs=self.environ.copy()
    # determine available resources...
    machines=self.find_resources(hogmin,hogmax)

    with open(self.file['queue']) as f:
      queue=f.readlines()
    done=[];running=[]
    with open(self.file['queue'],'w') as f:
      pass # clear queue since we read it. 
      # This way anything written while scattering won't be lost
    try: # catch any errors so we don't lose anything from the queue 
      # We need to keep everything in one list (queue) so that try/finally can write without losing anything
      for iq,q in enumerate(queue):
        q=q.strip()
        if not q:
          done.append(iq)
          continue # it's just a blank line
        elif q[0]=='#': 
          # q is already broken
          sys.stderr.write(" %s is broken\n"%(q[1:]))
          continue
        # check that job is a valid path, if not store in broken list...
        elif not os.path.exists(q):
          # q is not a valid path
          sys.stderr.write(" %s does not exist\n"%(q))
          queue[iq]='#'+q
          continue
        elif not os.path.isfile(q+'/queue_boomerang'):
          # check that queue_boomerang file still exists
          # it may have been already submitted or just unqueued
          sys.stderr.write("====> %s does not have queue_boomerang file\n"%q)
          queue[iq]='#'+q
          continue
        elif not machines: 
          # nothing to submit to
          break
        # if we're here, let's submit
        # machines says which mach to submit to
        os.remove("%s/queue_boomerang"%q) # this way the file is deleted before rsync; also, we don't resubmit if it throws an error and still actually goes through 
        if not self.submit_run(q,machines.pop(0)):
          done.append(iq)
        else: # an error came up in submission
          sys.stderr.write("====> Error submitting %s\n"%q)
          queue[iq]='#'+q
    finally: # no matter what happened, rewrite the queue
      with open(self.file["queue"],'a') as OUT: # append in case something was added to the queue in the meantime
        OUT.write('\n'.join([q.strip() for iq,q in enumerate(queue) if not iq in done])+'\n')

  def remove(self,inp=''):
    """Removes a directory from queue, running, and finished files. If it is running on the remote computer, it IS NOT killed. If it is queued, running, or finished, the associated local files (queue_boomerang, running_boomerang, bsubmit.sh) are NOT deleted. Input can be a filename with a list of directories, or blank for the current directory"""
    if inp:
      if os.path.exists(inp):inp=[s.strip() for s in open(inp).readlines() if s.split()]
      else:print "file %s does not exist"%inp;sys.exit()
    else:inp=[os.getcwd()]
    for direc in inp:
      if direc[-1]!='/': direc+='/'
      with open(self.file['queue']) as f:
        queueds = [s.strip() for s in f.readlines()]
      with open(self.file['running']) as f:
        runnings = [s.strip() for s in f.readlines()]
      with open(self.file['finished']) as f:
        finisheds = [s.strip() for s in f.readlines()]
      noaction=True
      if direc in queueds:
        noaction=False
        sys.stderr.write("Removing %s from queue.\n"%direc)
        newqueue = [s for s in queueds if not s==direc]
        with open(self.file['queue'],'w') as out:
          out.write('\n'.join(newqueue)+'\n')
      if direc in runnings:
        noaction=False
        sys.stderr.write("Removing %s from running.\n"%direc)
        newrunning = [s for s in runnings if not s==direc]
        with open(self.file['running'],'w') as out:
          out.write('\n'.join(newrunning)+'\n')
      if direc in finisheds:
        noaction=False
        sys.stderr.write("Removing %s from finished.\n"%direc)
        newfinished = [s for s in finisheds if not s==direc]
        with open(self.file['finished'],'w') as out:
          out.write('\n'.join(newfinished)+'\n')
      if noaction: 
        sys.stderr.write("Did not find %s in queue, running, or finished.\n"%direc)

if __name__=="__main__":
  from inspect import getargspec
  av=inputs(""" name=str noremote=str
      query=str status=str queryall=str noquery=str quick=str
      prep=str kill=str 
      cat=str go=str sum=str check=str goto=str
      force=str forceget=str get=str queue=str scatter=str remove=str
      timeout=int hogmin=int hogmax=int 
      """)
  av.update(sys.argv[1:])
  def gsp(func):return av.get_dict(getargspec(func)[0][1:])
  qq=boomerang(av.name if dir(av).count('name') else '',noremote=(dir(av).count('noremote')!=0))

  # Get: force_get, get
  if dir(av).count('forceget'):   qq.force_get()
  elif dir(av).count('get'):       qq.get()

  # Manipulate: go, remove, kill, prep, queue
  elif dir(av).count('go'):        qq.go(av.go)
  elif dir(av).count('remove'):    qq.remove(av.remove)
  elif dir(av).count('kill'):      qq.kill(av.kill)
  elif dir(av).count('prep'):      qq.prep()
  elif dir(av).count('queue'):     qq.add_to_queue(av.queue)

  # Read information without query: cat, check, sum
  elif dir(av).count("sum"):       qq.print_summary()
  elif dir(av).count("cat"):       qq.cat(av.cat)
  elif dir(av).count("check"):     qq.check(True if dir(av).count('goto') else False)

  else: # for these we need info about cluster status
    # Get info: queryall, query. We need to query for the other functions anyway
    timeout=av.timeout if dir(av).count('timeout') else 2
    if dir(av).count("queryall"):    qq.query_all(timeout)
    elif dir(av).count('noquery'):   qq.fetch_local()
    #elif dir(av).count('query'):     qq.query(10,timeout)
    else:                            qq.query(10,timeout,quick=(dir(av).count('quick')!=0) )

    # Print info: status
    if dir(av).count("status"):       qq.print_status()
    # Give: force, scatter, queue
    elif dir(av).count("force"):      qq.force(av.force)
    elif dir(av).count("scatter"):
      qq.scatter( av.hogmin if dir(av).count('hogmin') else 25, \
          av.hogmax if dir(av).count('hogmax') else 75)

    # write the final state of the queue...
    qq.write_local()
    if hasattr(qq,'remote') and qq.remote and qq.changed: qq.write_remote()
