#!/usr/bin/env python 
import os,sys  
from mysub import pmatrix,inputs,gs_orthog,lstr
import numpy as np
from periodica import structure
import glob
import scipy.linalg as la
from densm_vasp import get_dm
import cmath


mass={
'H':1.0079,
'He':4.0026,
'Li':6.941,
'Be':9.0122,
'B':10.811,
'C':12.0107,
'N':14.0067,
'O':15.9994,
'F':18.9984,
'Ne':20.1797,
'Na':22.9897,
'Mg':24.305,
'Al':26.9815,
'Si':28.0855,
'P':30.9738,
'S':32.065,
'Cl':35.453,
'K':39.0983,
'Ar':39.948,
'Ca':40.078,
'Sc':44.9559,
'Ti':47.867,
'V':50.9415,
'Cr':51.9961,
'Mn':54.938,
'Fe':55.845,
'Ni':58.6934,
'Co':58.9332,
'Cu':63.546,
'Zn':65.39,
'Ga':69.723,
'Ge':72.64,
'As':74.9216,
'Se':78.96,
'Br':79.904,
'Kr':83.8,
'Rb':85.4678,
'Sr':87.62,
'Y':88.9059,
'Zr':91.224,
'Nb':92.9064,
'Mo':95.94,
'Tc':98.0,
'Ru':101.07,
'Rh':102.9055,
'Pd':106.42,
'Ag':107.8682,
'Cd':112.411,
'In':114.818,
'Sn':118.71,
'Sb':121.76,
'I':126.9045,
'Te':127.6,
'Xe':131.293,
'Cs':132.9055,
'Ba':137.327,
'La':138.9055,
'Ce':140.116,
'Pr':140.9077,
'Nd':144.24,
'Pm':145.0,
'Sm':150.36,
'Eu':151.964,
'Gd':157.25,
'Tb':158.9253,
'Dy':162.5,
'Ho':164.9303,
'Er':167.259,
'Tm':168.9342,
'Yb':173.04,
'Lu':174.967,
'Hf':178.49,
'Ta':180.9479,
'W':183.84,
'Re':186.207,
'Os':190.23,
'Ir':192.217,
'Pt':195.078,
'Au':196.9665,
'Hg':200.59,
'Tl':204.3833,
'Pb':207.2,
'Bi':208.9804,
'Po':209.0,
'At':210.0,
'Rn':222.0,
'Fr':223.0,
'Ra':226.0,
'Ac':227.0,
'Pa':231.0359,
'Th':232.0381,
'Np':237.0,
'U':238.0289,
'Am':243.0,
'Pu':244.0,
'Cm':247.0,
'Bk':247.0,
'Cf':251.0,
'Es':252.0,
'Fm':257.0,
'Md':258.0,
'No':259.0,
'Rf':261.0,
'Lr':262.0,
'Db':262.0,
'Bh':264.0,
'Sg':266.0,
'Mt':268.0,
'Rg':272.0,
'Hs':277.0,
}


#'H':1.008,'C':12.011,'O':15.9994,'Pu':244.,'Pb':207.2,'Te':127.6,'Mo':95.96,'S':32.06,'Co':58.933,'Li':6.94}
#names=['Pu','O','O']
#names=['Pu','Pu','O','O','O','O']


class phonons_k:
  def __init__(self,pos,kvec,supa=''):
    if type(kvec)==str: kvec=[eval(s) for s in kvec.split(',')]
    if len(kvec)!=3:print "Wrong k-vector format..."; sys.exit()
    self.kvec=[float(s) for s in kvec]
    self.pos=structure(pos)
    # create supercell corresponding to the given k-vector...
    if not supa: supa=self.pos.get_supa_fromk(self.kvec)
    print "Creating Supercell %s"%supa
    self.pos.supa(supa)
    self.pos.title+=" "+supa
    self.norbs=len(self.pos.primcell.orbitals)
    self.get_basis()
    self.remove_acoustic()
    self.get_displaced(vec=0,delta=0.04)


  def get_basis(self):
    # start with the abstract vector...
    self.ab_basis=[{} for s in range(self.norbs)]
    #print self.pos.primcell.orbitals
    #print self.pos.supa_orbitals
    for i,ii in enumerate(self.pos.primcell.orbitals): 
      temp=[s for s in self.pos.supa_orbitals if s[0]==ii[0] and s[-1]==ii[-1] ]
      tempc=[np.exp(1j*2*np.pi*np.dot(s[1],self.kvec)) for s in temp]
      const=np.sqrt(np.dot(np.array(tempc).conj(),tempc))
      tempc=[s/const for s in tempc]
      for j,jj in enumerate(temp):self.ab_basis[i][jj]=tempc[j]
    # now make a vector...
    self.basis=[]
    for i,ii in enumerate(self.ab_basis):
      self.basis.append([ ii[s] if ii.keys().count(s) else 0.  for s in 
                         [(m[0],m[1],k) for m in self.pos.supa_patom for k in ['p_x','p_y','p_z']]   ])
    self.basis=np.array(self.basis)

  def remove_acoustic(self):
    if not tuple(self.kvec)==(0,0,0):return
    natoms=self.pos.natoms
    N=1./np.sqrt(float(natoms))
    acoustic=[]
    acoustic.append([N,0,0]*natoms)
    acoustic.append([0,N,0]*natoms)
    acoustic.append([0,0,N]*natoms)
    # now we need to orthogonalize out the acoustic modes...
    for x,xx in enumerate(acoustic):
      for v,vv in enumerate(self.basis):
        if np.dot(xx,vv)==1: break
        elif 1*10**-6<np.dot(xx,vv)<1:self.basis[v]=[s for s in xx];break
    self.basis=gs_orthog(np.array(self.basis))
    # now remove the acoustic modes... and tack them on at end...
    phonons=[]
    for v,vv in enumerate(self.basis):
      temp=0
      for x,xx in enumerate(acoustic): temp+=np.dot(xx,vv)
      if temp<=10**-6:phonons.append(vv)
    #pmatrix(phonons,6)
    self.basis=np.array([s for s in phonons+acoustic])


  def get_displaced(self,vec,delta,output=''):
    delta=delta/max(np.array(self.basis[vec]).real)
    disp=np.reshape(self.basis[vec],(self.pos.natoms,-1)).real
    bb=self.pos.copy()
    bb.atc+=disp*delta
    bb.at=bb.atc*1.
    bb.coord='c'
    return bb

  def construct_dispk(self,delta,wdir='',inp_path=''):
    """This method will make a directory for each 3N phonon displacements for wavevector k.
       If inp_path is given, the code looks in that directory for standard input files (ie. INCAR, KPOINTS, etc).
    """
    if not wdir:wdir="k_%.3f_%.3f_%.3f"%tuple(self.kvec)
    if os.path.exists(wdir):print "directory %s already exists... stopping."%wdir;sys.exit()
    if not inp_path or inp_path[0:2]=='..':inp_path='../'+inp_path
    os.mkdir(wdir); os.chdir(wdir)
    delta=float(delta)
    os.mkdir("unperturbed")
    self.pos.pposcar("unperturbed/POSCAR")
    for i in range(self.norbs): 
      output='dist_%.2i'%i
      if not os.path.exists(output):
        os.mkdir(output)
        self.get_displaced(i,delta).pposcar(output+'/POSCAR')
        for f in ['INCAR','KPOINTS','POTCAR','dc.input','density_matrix_up.inp','density_matrix_dn.inp']:
          if os.path.exists(inp_path+'/'+f):
            os.system("cp %s/%s %s"%(inp_path,f,output))
            os.system("cp %s/%s unperturbed/"%(inp_path,f))
      else:print "%s already exists..."%output

  def units(self,inp=''):
    if inp=="meV":     return 1000 # meV/eV
    elif inp=='THz':   return 241.8 # THz/eV
    elif inp=='cm-1':  return 8065.6 # cm^-1/eV
    elif inp=='':      return 1. # eV
    else:print "units not recognized...";sys.exit()
    
  def compute_phonons(self,dirlist='dist_*',units=''):
    """Computes the phonons using frozen phonon displacement runs. """
    # compute the mass matrix...
    mass_matrix=[[ (kk if k==n else 0) for k,kk in enumerate(sum([[mass[s]]*3 for s in self.pos.names  ],[])) ] for n in range(self.pos.natoms*3)]
    mass_matrix=np.dot(self.basis,np.dot(mass_matrix,np.transpose(self.basis).conj())).real
    #pmatrix(mass_matrix)

    # compute the force constants...
    dynm=np.zeros((self.norbs,self.norbs),complex)
    for i in sorted(glob.glob(dirlist)): 
      pf=int(i.split("_")[-1])
      if os.path.exists('%s/OUTCAR'%i):
          aa=structure("%s/POSCAR"%i)
          delta=np.dot(self.basis,np.reshape(aa.atc,(1,-1))[0]-np.reshape(self.pos.atc,(1,-1))[0])
          ind=[s for s,ss in enumerate(delta) if abs(ss)>10**-6]
          if not ind:print "Error: no distortion found for %s..."%i;sys.exit()
          else: 
            if len(ind)>1:print "More than one displacement detected for %s... stopping."%i;sys.exit()
            else:ind=ind[0]
            if ind!=pf:print "Calculated mode does not match directory %s label..."%i;sys.exit()
            else:delta=delta[ind]
          print "Distortion %s delta=%.4f+i*%.4f ind=%s"%(i,delta.real,delta.imag,ind)
          force=[[float(s) for s in kk.split()[3:]] for k,kk in enumerate(os.popen("grep -A %s TOTAL-FORCE %s/OUTCAR"%(aa.natoms+1,i)).readlines()[2:])]
          force=np.reshape(force,(1,-1))[0]
          dynm[pf]=[((ss+dynm[pf,s])/2 if dynm[pf,s]>10**-3 else ss) for s,ss in enumerate(-1*np.dot(self.basis.conj(),force)/delta) ]
          #pmatrix([np.dot(np.array(self.basis).conj(),force).real])
          #pmatrix([np.dot(np.array(self.basis).conj(),force).imag])
      else:print "directory %s not started..."%i 
     

    print "Final Dynamical Matrix real..."
    pmatrix(dynm.real)
    print "Final Dynamical Matrix imag..."
    pmatrix(dynm.imag)
    

    dynm=(dynm+dynm.transpose().conj())/2.
    #eig=[np.sqrt(s/(1.66054*10**-27*6.2415*10**-2))*6.582*10**-16*self.units(units) if s>=0 else s for s in la.eigh(dynm,mass_matrix)[0]]
    eig=[cmath.sqrt(s/(1.66054*10**-27*6.2415*10**-2))*6.582*10**-16*self.units(units) for s in la.eigh(dynm,mass_matrix)[0]]
    eigf=la.eigh(dynm,mass_matrix)[1]

    # normalize the eigenfunctions with the mass matrix...
    #temp=eigf[:,0]
    #tempmass=np.array([np.sqrt(mass[s]) for s in self.pos.names])
    #eigf=np.array([s*tempmass for s in eigf.transpose()]).transpose()

    eigf=np.transpose(eigf)
    print "Eigenfunctions are rows"
    print "Real"
    pmatrix(eigf.real,4)
    print 'Imag'
    pmatrix(eigf.imag,4)
    print "Eigenvalues"
    #print " ".join(["%.4f"%s for s in eig])
    print lstr([s.real if abs(s.real)>abs(s.imag) else -abs(s.imag) for s in eig])
    
     
      


if __name__=='__main__':
  av=inputs("""source='vasp' pos='poscar' delta=float k=str path=str inc=str kdens=5 
               so=str nosym=str ldapu=str kmesh='auto' dc=float units='cm-1'
               usedm=str nsp=str supa='' verb=str inp_path=''  """)
  av.update(sys.argv[1:])

  
  # build the displacement directories...
  if dir(av).count("delta"):
    kv=phonons_k(av.pos,av.k,supa=av.supa)
    kv.construct_dispk(av.delta,inp_path=av.inp_path)
    sys.exit()

  # process the displacement directories...
  if dir(av).count("k"): 
    kv=phonons_k(av.pos,av.k,supa=av.supa)
  else: # automatically find kpoint from directory name...
    temp=os.getcwd().split("/")[-1]
    if temp.count("k_"):
      kk=tuple([float(s) for s in temp[2:].split("_")])
      kv=phonons_k(av.pos,kk,supa=av.supa)
    else:print "No k-vector supplied or found...";sys.exit()

  # print things depending on verbosity setting...
  if dir(av).count('verb'):
    print 'real'
    pmatrix(np.array(kv.basis).real)
    print 'imag'
    pmatrix(np.array(kv.basis).imag)

  # compute the phonons...
  kv.compute_phonons(units=av.units)

