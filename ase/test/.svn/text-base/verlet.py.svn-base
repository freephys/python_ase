from ase import *
from ase.calculators import TestPotential
np.seterr(all='raise')
a = Atoms('4N', 
          positions=[(0, 0, 0),
                     (1, 0, 0),
                     (0, 1, 0),
                     (0.1, 0.2, 0.7)],
          calculator=TestPotential())
print a.get_forces()
md = VelocityVerlet(a, dt=0.005)
def f():
    print a.get_potential_energy(), a.get_total_energy()
md.attach(f, 500)
traj = PickleTrajectory('4N.traj', 'w', a)
md.attach(traj.write, 100)
print md.observers
md.run(steps=10000)
qn = QuasiNewton(a)
qn.attach(traj.write)
qn.run()
