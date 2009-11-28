"""This module defines an ASE interface to SIESTA.

http://www.uam.es/departamentos/ciencias/fismateriac/siesta
"""

import os
from os.path import join, isfile, islink

import numpy as np

from ase.data import chemical_symbols


class Siesta:
    """Class for doing SIESTA calculations.

    The default parameters are very close to those that the SIESTA
    Fortran code would use.  These are the exceptions::
    
      calc = Siesta(label='siesta', xc='LDA', pulay=5, mix=0.1)

    Use the set_fdf method to set extra FDF parameters::
    
      calc.set_fdf('PAO.EnergyShift', 0.01 * Rydberg)

    """
    def __init__(self, label='siesta', xc='LDA', kpts=None, nbands=None,
                 width=None, meshcutoff=None, charge=None,
                 pulay=5, mix=0.1, maxiter=120,
                 basis=None, ghosts=[],
                 write_fdf=True):
        """Construct SIESTA-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.fdf, label.txt, ...).
            Default is 'siesta'.
        xc: str
            Exchange-correlation functional.  Must be one of LDA, PBE,
            revPBE, RPBE.
        kpts: list of three int
            Monkhost-Pack sampling.
        nbands: int
            Number of bands.
        width: float
            Fermi-distribution width in eV.
        meshcutoff: float
            Cutoff energy in eV for grid.
        charge: float
            Total charge of the system.
        pulay: int
            Number of old densities to use for Pulay mixing.
        mix: float
            Mixing parameter between zero and one for density mixing.
        write_fdf: bool
            Use write_fdf=False to use your own fdf-file.

        Examples
        ========
        Use default values:
        
        >>> h = Atoms('H', calculator=Siesta())
        >>> h.center(vacuum=3.0)
        >>> e = h.get_potential_energy()
        
        """
        
        self.label = label#################### != out
        self.xc = xc
        self.kpts = kpts
        self.nbands = nbands
        self.width = width
        self.meshcutoff = meshcutoff
        self.charge = charge
        self.pulay = pulay
        self.mix = mix
        self.maxiter = maxiter
        self.basis = basis
        self.ghosts = ghosts
        self.write_fdf_file = write_fdf
    
        self.converged = False
        self.fdf = {}

    def update(self, atoms):
        if (not self.converged or
            len(self.numbers) != len(atoms) or
            (self.numbers != atoms.get_atomic_numbers()).any()):
            self.initialize(atoms)
            self.calculate(atoms)
        elif ((self.positions != atoms.get_positions()).any() or
              (self.pbc != atoms.get_pbc()).any() or
              (self.cell != atoms.get_cell()).any()):
            self.calculate(atoms)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(self.numbers):
            if a in self.ghosts:
                Z = -Z
            if Z not in self.species:
                self.species.append(Z)

        if 'SIESTA_PP_PATH' in os.environ:
            pppaths = os.environ['SIESTA_PP_PATH'].split(':')
        else:
            pppaths = []

        for Z in self.species:
            symbol = chemical_symbols[abs(Z)]
            name = symbol + '.vps'
            name1 = symbol + '.psf'
            found = False
            for path in pppaths:
                filename = join(path, name)
                filename1 = join(path,name1)
                if isfile(filename) or islink(filename):
                    found = True
                    if path != '.':
                        if islink(name) or isfile(name):
                            os.remove(name)
                        os.symlink(filename, name)

                elif isfile(filename1) or islink(filename1):
                    found = True
                    if path != '.':
                        if islink(name1) or isfile(name1):
                            os.remove(name1)
                        os.symlink(filename1, name1)
            if not found:
                raise RuntimeError('No pseudopotential for %s!' % symbol)

        self.converged = False
        
    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)

        if force_consistent:
            return self.efree
        else:
            # Energy extrapolated to zero Kelvin:
            return  (self.etotal + self.efree) / 2

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces.copy()
    
    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress.copy()

    def calculate(self, atoms):
        self.positions = atoms.get_positions().copy()
        self.cell = atoms.get_cell().copy()
        self.pbc = atoms.get_pbc().copy()

        if self.write_fdf_file:
            self.write_fdf(atoms)

        siesta = os.environ['SIESTA_SCRIPT']
        locals = {'label': self.label}
        execfile(siesta, {}, locals)
        exitcode = locals['exitcode']
        if exitcode != 0:
            raise RuntimeError(('Siesta exited with exit code: %d.  ' +
                                'Check %s.txt for more information.') %
                               (exitcode, self.label))
        
        self.read()

        self.converged = True
        
    def set_fdf(self, key, value):
        """Set FDF parameter."""
        self.fdf[key] = value
                
    def write_fdf(self, atoms):
        """Write input parameters to fdf-file."""
        fh = open(self.label + '.fdf', 'w')

        fdf = {
            'SystemLabel': self.label,
            'AtomicCoordinatesFormat': 'Ang',
            'LatticeConstant': 1.0,
            'NumberOfAtoms': len(atoms),
            'MeshCutoff': self.meshcutoff,
            'NetCharge': self.charge,
            'ElectronicTemperature': self.width,
            'NumberOfEigenStates': self.nbands,
            'DM.UseSaveDM': self.converged,
            'PAO.BasisSize': self.basis,
            'SolutionMethod': 'diagon',
            'DM.NumberPulay': self.pulay,
            'DM.MixingWeight': self.mix,
            'MaxSCFIterations' : self.maxiter
            }
        
        if self.xc != 'LDA':
            fdf['xc.functional'] = 'GGA'
            fdf['xc.authors'] = self.xc

        magmoms = atoms.get_initial_magnetic_moments()
        if magmoms.any():
            fdf['SpinPolarized'] = True
            fh.write('%block InitSpin\n')
            for n, M in enumerate(magmoms):
                if M != 0:
                    fh.write('%d %.14f\n' % (n + 1, M))
            fh.write('%endblock InitSpin\n')
        
        fdf['Number_of_species'] = len(self.species)

        fdf.update(self.fdf)

        for key, value in fdf.items():
            if value is None:
                continue

            if isinstance(value, list):
                fh.write('%%block %s\n' % key)
                for line in value:
                    fh.write(line + '\n')
                fh.write('%%endblock %s\n' % key)
            else:
                unit = keys_with_units.get(fdfify(key))
                if unit is None:
                    fh.write('%s %s\n' % (key, value))
                else:
                    if 'fs**2' in unit:
                        value /= fs**2
                    elif 'fs' in unit:
                        value /= fs
                    fh.write('%s %f %s\n' % (key, value, unit))

        fh.write('%block LatticeVectors\n')
        for v in self.cell:
            fh.write('%.14f %.14f %.14f\n' %  tuple(v))
        fh.write('%endblock LatticeVectors\n')

        fh.write('%block Chemical_Species_label\n')
        for n, Z in enumerate(self.species):
            fh.write('%d %s %s\n' % (n + 1, Z, chemical_symbols[abs(Z)]))
        fh.write('%endblock Chemical_Species_label\n')

        fh.write('%block AtomicCoordinatesAndAtomicSpecies\n')
        a = 0
        for pos, Z in zip(self.positions, self.numbers):
            if a in self.ghosts:
                Z = -Z
            a += 1
            fh.write('%.14f %.14f %.14f' %  tuple(pos))
            fh.write(' %d\n' % (self.species.index(Z) + 1))
        fh.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')

        if self.kpts is not None:
            fh.write('%block kgrid_Monkhorst_Pack\n')
            for i in range(3):
                for j in range(3):
                    if i == j:
                        fh.write('%d ' % self.kpts[i])
                    else:
                        fh.write('0 ')
                fh.write('%.1f\n' % (((self.kpts[i] + 1) % 2) * 0.5))
            fh.write('%endblock kgrid_Monkhorst_Pack\n')
        
        fh.close()
        
    def read(self):
        """Read results from SIESTA's text-output file."""
        text = open(self.label + '.txt', 'r').read().lower()
        assert 'error' not in text
        lines = iter(text.split('\n'))

        # Get the number of grid points used:
        for line in lines:
            if line.startswith('initmesh: mesh ='):
                self.grid = [int(word) for word in line.split()[3:8:2]]
                break

        # Stress:
        for line in lines:
            if line.startswith('siesta: stress tensor (total) (ev/ang**3):'):
                self.stress = np.empty((3, 3))
                for i in range(3):
                    self.stress[i] = [float(word)
                                      for word in lines.next().split()]
                break
        else:
            raise RuntimeError

        # Energy:
        for line in lines:
            if line.startswith('siesta: etot    ='):
                self.etotal = float(line.split()[-1])
                self.efree = float(lines.next().split()[-1])
                break
        else:
            raise RuntimeError

        # Forces:
        lines = open(self.label + '.FA', 'r').readlines()
        assert int(lines[0]) == len(self.numbers)
        assert len(lines) == len(self.numbers) + 1
        self.forces = np.array([[float(word)
                                 for word in line.split()[1:4]]
                                for line in lines[1:]])
        

def fdfify(key):
    return key.lower().replace('_', '').replace('.', '').replace('-', '')


keys_with_units = {
    'paoenergyshift': 'eV',
    'zmunitslength': 'Bohr',
    'zmunitsangle': 'rad',
    'zmforcetollength': 'eV/Ang',
    'zmforcetolangle': 'eV/rad',
    'zmmaxdispllength': 'Ang',
    'zmmaxdisplangle': 'rad',
    'meshcutoff': 'eV',
    'dmenergytolerance': 'eV',
    'electronictemperature': 'eV',
    'oneta': 'eV',
    'onetaalpha': 'eV',
    'onetabeta': 'eV',
    'onrclwf': 'Ang',
    'onchemicalpotentialrc': 'Ang',
    'onchemicalpotentialtemperature': 'eV',
    'mdmaxcgdispl': 'Ang',
    'mdmaxforcetol': 'eV/Ang',
    'mdmaxstresstol': 'eV/Ang**3',
    'mdlengthtimestep': 'fs',
    'mdinitialtemperature': 'eV',
    'mdtargettemperature': 'eV',
    'mdtargetpressure': 'eV/Ang**3',
    'mdnosemass': 'eV*fs**2',
    'mdparrinellorahmanmass': 'eV*fs**2',
    'mdtaurelax': 'fs',
    'mdbulkmodulus': 'eV/Ang**3',
    'mdfcdispl': 'Ang',
    'warningminimumatomicdistance': 'Ang',
    'rcspatial': 'Ang',
    'kgridcutoff': 'Ang',
    'latticeconstant': 'Ang'}

