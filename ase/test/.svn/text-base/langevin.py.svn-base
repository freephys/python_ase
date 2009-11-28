"""Test the Langevin dynamics.

Tests Langevin dynamics using the EMT Copper potential.

For time reasons, long term temperature average is only calculated with asap.
"""

import time
import numpy as np
from ase import Atoms
from ase.calculators import EMT, ASAP
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin

try:
    from Scientific.Functions.LeastSquares import leastSquaresFit
except ImportError:
    usescipy = False
else:
    usescipy = True

try:
    import Asap
except ImportError:
    useasap = False
else:
    useasap = True

if useasap:
    nequil = 1000
    nequilprint = 25
    nsteps = 20000
    nprint = 250
    reltol = 0.025
else:
    nsteps = 2000
    nequil = 500
    nequilprint = 25
    nprint = 25
    reltol = 0.05
nminor = 25
timestep = 0.5

# Set up atoms in a regular simple-cubic lattice.
a = 3.5
atoms = Atoms('Cu4', a * np.array([[0,0,0],
                                   [0.5, 0.5, 0],
                                   [0.5, 0, 0.5],
                                   [0, 0.5, 0.5]]),
              cell=(a, a, a))
atoms *= (1, 5, 5)

if useasap:
    atoms.set_calculator(ASAP())
else:
    atoms.set_calculator(EMT())
    print """
WARNING: You have chosen to use the Python-implemented PairPotential.
It is exceedingly slow, and will not accumulate enough statistics for
a good test of the dynamics.  For orders of magnitude better
perfomance, make sure Asap is installed!"""
    
# Least-squares fit model during equilibration
def targetfunc(params, t):
    return params[0] * np.exp(-params[1] * t) + params[2]

def test(temp, frict):
    output = file('Langevin.dat', 'w')
    
    # Make a small perturbation of the momenta
    atoms.set_momenta(1e-6 * np.random.random([len(atoms), 3]))
    print 'Initializing ...'
    predyn = VelocityVerlet(atoms, 0.5)
    predyn.run(2500)

    dyn = Langevin(atoms, timestep, temp, frict)
    print ''
    print ('Testing Langevin dynamics with T = %f eV and lambda = %f' %
           (temp, frict))
    ekin = atoms.get_kinetic_energy()/len(atoms)
    print ekin
    output.write('%.8f\n' % ekin)

    print 'Equilibrating ...'

    # Initial guesses for least-squares fit
    a = 0.04
    b = 2*frict
    c = temp
    params = (a,b,c)
    fitdata = [(0, 2.0 / 3.0 * ekin)]

    tstart = time.time()
    for i in xrange(1,nequil+1):
        dyn.run(nminor)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        fitdata.append((i*nminor*timestep, 2.0/3.0 * ekin))
        if usescipy and i % nequilprint == 0:
            (params, chisq) = leastSquaresFit(targetfunc, params, fitdata)
            print '%.6f  T_inf = %.6f (goal: %f), tau = %.2f,  k = %.6f' % \
                  (ekin, params[2], temp, 1.0/params[1], params[0])
        output.write('%.8f\n' % ekin)
    tequil = time.time() - tstart
    print 'This took %s minutes.' % (tequil / 60)
    output.write('&\n')
    assert abs(temp-params[2]) < 0.25*temp, 'Least-squares fit is way off'
    assert nequil*nminor*timestep > 3.0/params[1], 'Equiliberation was too short'
    fitdata = np.array(fitdata)

    print 'Recording statistical data - this takes ten times longer!'
    temperatures = []
    tstart = time.time()
    for i in xrange(1,nsteps+1):
        dyn.run(nminor)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        temperatures.append(2.0/3.0 * ekin)
        if i % nprint == 0:
            tnow = time.time() - tstart
            tleft = (nsteps-i) * tnow / i
            print '%.6f    (time left: %.1f minutes)' % (ekin, tleft/60)
        output.write('%.8f\n' % ekin)
    output.write('&\n')
    output.close()

    temperatures = np.array(temperatures)
    mean = sum(temperatures) / len(temperatures)
    print 'Mean temperature:', mean, 'eV'
    print
    print 'This test is statistical, and may in rare cases fail due to a'
    print 'statistical fluctuation.'
    print
    assert abs(mean - temp) <= reltol*temp, 'Deviation is too large.'
    print 'Mean temperature:', mean, ' in ', temp, ' +/- ', reltol*temp

    return fitdata, params, temperatures

if __name__ in ['__main__', '__builtin__']:
    temp = 0.01
    frict = 0.001
    if useasap:
        fitdata, params, temperatures = test(temp, frict)
