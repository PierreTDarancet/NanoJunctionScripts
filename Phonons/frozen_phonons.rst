frozen_phonon.py
==============


Intro
^^^^


frozen_phonon.py is a python script working with vasp that:

	(1) Generates the geometries for the calculation of dynamical matrix (second derivative of the total energy wrt the displacement field)
	(2) Calculates the dynamical matrix from the change in forces (first derivative wrt displacement) as a function of displacement. 


Before you start:
^^^^^^^^^^^^^^

	(1) Check the mass dictionary present at the beginning of ``frozen_phonon.py``. It reads, ``mass={'H':1.008,'Li':6.94,'C':12.011,'S':32.06,'Co':58.933,'Mo':95.96,'Te':127.6,'O':15.9994,'Pb':207.2,'Pu':244.,'Ta':180.9479}``. Since it obviously doesn't contain all the needed atomic masses, be sure to add the one you need. 
	(2) Be aware that ``frozen_phonon.py`` makes use of several modules (mysub.py, periodica.py), meaning a ``export PYTHONPATH=~/pythondir`` should be called in the bashrc or in the session. 
	(3) Prepare a POSCAR file where a fourth column is present and reads "AtomicSymbol:p" (for example Ta:p). In that case, the orbital doesn't refer to the valence/conduction bands but to the symmetry of the deformation potential. Name this file "poscar" (lower case). For all displacements, use the orbital code "p".  

Step (1) will be called by 

.. code-block:: bash 

	$cd workdir/
	$if [ -f ./poscar] frozen_phonon.py k=1/3.,1/3.,0 delta=0.04


``frozen_phonon.py`` will create 3*Natoms + 1 subdirectories in a k_kx_ky_kz directory. In every single one of these subdirectories besides the 'undistorted', one atom is moved from the distance delta (in Angstroms) along one of the axis in the supercell defined by the q-vector of the phonon. Ideally, delta should be as small as numerically possible (linear response). At the very least the convergence should be tested wrt to this parameter. 


The postprocessing Step (2)
^^^^^^^^^^^^^^^^^^^^^^^^^

It will by called by

.. code-block:: bash 

	$cd workdir/k_0.333_0.333_0.000/
	$frozen_phonon.py pos=../poscar k=1/3.,1/3.,0 

With this syntax, the code is going to parse through all the dist_* directories and grep the forces on each ion, and build the dynamical matrix :math:`\frac{d^2 E}{d\mathbf{R}_A d\mathbf{R}_B}`. Diagonalizing this matrix gives the eigenmode energies as well as their wavefunctions projected on the ions. 

The code will output a lot of information:

First, the parsing:: 


	Creating Supercell 1 0 0   0 1 0   0 0 1
	Distortion dist_00 delta=0.0405+i*0.0000 ind=0
	Distortion dist_01 delta=0.0405+i*0.0000 ind=1
	...
	Eigenfunctions are rows
	Real
	... some long line of numbers
	Imag::
	...  some long line of numbers, a lot of zeros at G
	 Eigenvalues
	  -1.3499  -0.4228  -0.2559  42.1015 ...

Beware of the units! They can be changed with the flag ``units=``. The default is ``units='cm-1'`` . 


Troubleshooting:
^^^^^^^^^^^^^^^


supa error: AttributeError: structure instance has no attribute 'get_supa_fromk'
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

I have been getting the following bug::

	$~/bin/frozen_phonon.py pos=../poscar k=0,0,0 delta=0.04
	Traceback (most recent call last):
	  File "/home/darancet/bin/frozen_phonon.py", line 184, in <module>
	    kv=phonons_k(av.pos,av.k,supa=av.supa)
	  File "/home/darancet/bin/frozen_phonon.py", line 25, in __init__
	    if not supa: supa=self.pos.get_supa_fromk(self.kvec)
	AttributeError: structure instance has no attribute 'get_supa_fromk'

Supa doesn't get which q-point I am trying to calculate. (Even for the Gamma point). This can be solved by adding a supa flag to the initial call, e.g. ``/bin/frozen_phonon.py pos=../poscar k=0,0,0 supa="1 0 0 0 1 0 0 0 1" delta=0.04``.

Code doesn't understand it needs to postprocess
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This one hurts but withdrawing the flag ``delta=``  is what switches the code to postprocessing mode 

supa error: did not find all translations within supercell... algorithm problem
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Use the latest version of periodica from the software repository. The supercell atom search needs the more thorough version, not the fast version. If you don't see "Trying enhanced search for translated atoms", you have the old code. 


