import numpy as np

from ase.md import MolecularDynamics


class VelocityVerlet(MolecularDynamics):
    def __init__(self, atoms, dt, trajectory=None):
        MolecularDynamics.__init__(self, atoms, dt, trajectory)
            
    def step(self, f):
        atoms = self.atoms
        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * f
        self.atoms.set_positions(self.atoms.get_positions() +
            self.dt * p / self.atoms.get_masses()[:,np.newaxis])
        # We need to store the momenta on the atoms before calculating
        # the forces, as in a parallel Asap calculation atoms may
        # migrate during force calculations, and the momenta need to
        # migrate along with the atoms.  For the same reason, we
        # cannot use self.masses in the line above.
        self.atoms.set_momenta(p)
        f = self.atoms.get_forces()
        atoms.set_momenta(self.atoms.get_momenta() + 0.5 * self.dt * f)
        return f
