import os
import pickle

from ase.calculators import SinglePointCalculator
from ase.atoms import Atoms
from ase.parallel import rank
from ase.utils import devnull
from ase.neb import NEB


class PickleTrajectory:
    "Reads/writes Atoms objects into a .traj file."
    # Per default, write these quantities
    write_energy=True
    write_forces=True
    write_stress=True
    write_momenta=True
    def __init__(self, filename, mode='r', atoms=None, master=None,
                 write_first_image=True):
        """A PickleTrajectory can be created in read, write or append mode.

        Parameters:

        filename:
            The name of the parameter file.  Should end in .traj.

        mode='r':
            The mode.

            'r' is read mode, the file should already exist, and
            no atoms argument should be specified.

            'w' is write mode.  If the file already exists, is it
            renamed by appending .bak to the file name.  The atoms
            argument specifies the Atoms object to be written to the
            file, if not given it must instead be given as an argument
            to the write() method.

            'a' is append mode.  It acts a write mode, except that
            data is appended to a preexisting file.

        atoms=None:
            The Atoms object to be written in write or append mode.

        master=None:
            Controls which process does the actual writing. The
            default is that process number 0 does this.  If this
            argument is given, processes where it is True will write.

        write_first_image=True:
        
            If this argument is True and atoms is specified, the atoms
            are written to the file immediately.  This is intended to
            write the initial configuration in a simulation.  Note
            that in append mode you probably want to set this to
            False, as the first configuration is likely to already be
            in the file.
        
        """
        self.offsets = []
        if master is None:
            master = (rank == 0)
        self.master = master
        self.set_atoms(atoms)
        self.open(filename, mode)

        if write_first_image and atoms is not None:
            self.write()
        
    def open(self, filename, mode):
        """Opens the file.

        For internal use only.
        """
        self.fd = filename
        if mode == 'r':
            if isinstance(filename, str):
                self.fd = open(filename, mode + 'b')
            self.read_header()
        elif mode == 'a':
            exists = True
            if isinstance(filename, str):
                exists = os.path.isfile(filename)
                self.fd = open(filename, mode + 'b+')
            if exists:
                self.read_header()
        elif mode == 'w':
            if self.master:
                if isinstance(filename, str):
                    if os.path.isfile(filename):
                        os.rename(filename, filename + '.bak')
                    self.fd = open(filename, 'wb')
            else:
                self.fd = devnull
        else:
            raise ValueError('mode must be "r", "w" or "a".')

    def set_atoms(self, atoms=None):
        """Associate an Atoms object with the trajectory.

        Mostly for internal use.
        """
        if atoms is not None and not hasattr(atoms, 'get_positions'):
            raise TypeError('"atoms" argument is not an Atoms object.')
        self.atoms = atoms

    def read_header(self):
        try:
            if self.fd.read(len('PickleTrajectory')) != 'PickleTrajectory':
                raise IOError('This is not a trajectory file!')
            d = pickle.load(self.fd)
        except EOFError:
            raise EOFError('Bad trajectory file.')
        self.pbc = d['pbc']
        self.numbers = d['numbers']
        self.tags = d.get('tags')
        self.constraints = d['constraints']
        self.offsets.append(self.fd.tell())

    def write(self, atoms=None):
        """Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        if atoms is None:
            atoms = self.atoms

        if isinstance(atoms, NEB):
            neb = atoms
            for image in neb.images:
                self.write(image)
            return

        if len(self.offsets) == 0:
            self.write_header(atoms)

        if atoms.has('momenta'):
            momenta = atoms.get_momenta()
        else:
            momenta = None

        d = {'positions': atoms.get_positions(),
             'cell': atoms.get_cell(),
             'momenta': momenta}


        if atoms.get_calculator() is not None:
            if self.write_energy:
                d['energy'] = atoms.get_potential_energy()
            if self.write_forces:
                assert(self.write_energy)
                d['forces'] = atoms.get_forces(apply_constraint=False)
            if self.write_stress:
                assert(self.write_energy)
                try:
                    d['stress'] = atoms.get_stress()
                except NotImplementedError:
                    pass

            try:
                if atoms.calc.get_spin_polarized():
                    d['magmoms'] = atoms.get_magnetic_moments()
            except (NotImplementedError, AttributeError):
                pass

        if 'magmoms' not in d and atoms.has('magmoms'):
            d['magmoms'] = atoms.get_initial_magnetic_moments()
            
        if self.master:
            pickle.dump(d, self.fd, protocol=-1)
        self.fd.flush()
        self.offsets.append(self.fd.tell())

    def write_header(self, atoms):
        self.fd.write('PickleTrajectory')
        if atoms.has('tags'):
            tags = atoms.get_tags()
        else:
            tags = None
        d = {'pbc': atoms.get_pbc(),
             'numbers': atoms.get_atomic_numbers(),
             'tags': tags,
             'constraints': atoms.constraints}
        pickle.dump(d, self.fd, protocol=-1)
        self.header_written = True
        self.offsets.append(self.fd.tell())
        
    def close(self):
        """Close the trajectory file."""
        self.fd.close()

    def __getitem__(self, i=-1):
        N = len(self.offsets)
        if 0 <= i < N:
            self.fd.seek(self.offsets[i])
            try:
                d = pickle.load(self.fd)
            except EOFError:
                raise IndexError
            if i == N - 1:
                self.offsets.append(self.fd.tell())
            try:
                magmoms = d['magmoms']
            except KeyError:
                magmoms = None    
            atoms = Atoms(positions=d['positions'],
                          numbers=self.numbers,
                          cell=d['cell'],
                          momenta=d['momenta'],
                          magmoms=magmoms,
                          tags=self.tags,
                          pbc=self.pbc,
                          constraint=[c.copy() for c in self.constraints])
            if 'energy' in d:
                calc = SinglePointCalculator(
                    d.get('energy', None), d.get('forces', None),
                    d.get('stress', None), magmoms, atoms)
                atoms.set_calculator(calc)
            return atoms

        if i >= N:
            for j in range(N - 1, i + 1):
                atoms = self[j]
            return atoms

        i = len(self) + i
        if i < 0:
            raise IndexError('Trajectory index out of range.')
        return self[i]

    def __len__(self):
        N = len(self.offsets) - 1
        while True:
            self.fd.seek(self.offsets[N])
            try:
                pickle.load(self.fd)
            except EOFError:
                return N
            self.offsets.append(self.fd.tell())
            N += 1

    def __iter__(self):
        del self.offsets[1:]
        return self

    def next(self):
        try:
            return self[len(self.offsets) - 1]
        except IndexError:
            raise StopIteration


def read_trajectory(filename, index=-1):
    traj = PickleTrajectory(filename, mode='r')

    if isinstance(index, int):
        return traj[index]
    else:
        return [traj[i] for i in range(len(traj))[index]]

def write_trajectory(filename, images):
    """Write image(s) to trajectory.

    Write also energy, forces, and stress if they are already
    calculated."""

    traj = PickleTrajectory(filename, mode='w')

    if not isinstance(images, (list, tuple)):
        images = [images]
        
    for atoms in images:
        # Avoid potentially expensive calculations:
        calc = atoms.get_calculator()
        if (calc is not None and
            (not hasattr(calc, 'calculation_required') or
             calc.calculation_required(atoms,
                                       ['energy', 'forces', 'stress']))):
            traj.write_energy=False
            traj.write_forces=False
            traj.write_stress=False
        
        traj.write(atoms)
    traj.close()
