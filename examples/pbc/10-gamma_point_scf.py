#!/usr/bin/env python

'''
Gamma point Hartree-Fock/DFT

The 2-electron integrals are computed using Poisson solver with FFT by default.
In most scenario, it should be used with pseudo potential.
'''

# Note import path which is different to molecule code
from pyscf.pbc import gto, scf, dft
import numpy

cell = gto.M(
    h = numpy.eye(3)*3.5668,
    atom = '''C     0.      0.      0.    
              C     0.8917  0.8917  0.8917
              C     1.7834  1.7834  0.    
              C     2.6751  2.6751  0.8917
              C     1.7834  0.      1.7834
              C     2.6751  0.8917  2.6751
              C     0.      1.7834  1.7834
              C     0.8917  2.6751  2.6751''',
    basis = 'gth-szv',
    pseudo = 'gth-pade',
    gs = [10]*3,
    verbose = 4,
)

mf = scf.RHF(cell)
ehf = mf.kernel()
print("HF energy (per unit cell) = %.17g" % ehf)

mf = dft.RKS(cell)
mf.xc = 'm06'
edft = mf.kernel()
print("DFT energy (per unit cell) = %.17g" % edft)

#
# By default, DFT use uniform cubic grids.  It can be replaced by atomic grids.
#
mf = dft.RKS(cell)
mf.grids = dft.gen_grid.BeckeGrids(cell)
mf.xc = 'bp86'
mf.kernel()
