import numpy as np
from ase.units import Bohr

def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4):
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    sep = '---------------'
    while sep not in fileobj.readline(): # Skip to after first seperator line
        pass

    for line in fileobj:
        if sep in line: # Stop at last seperator line
            break

        words = line.split()
        if len(words) != 6:
            raise IOError('Number of columns in ACF file incorrect!\n'
                          'Check that Bader program version >= 0.25')

        atom = atoms[int(words[0]) - 1]
        atom.charge = atom.number - float(words[4])

        if displacement is not None: # check if the atom positions match
            xyz = np.array([float(w) for w in words[1:4]]) * Bohr
            assert np.linalg.norm(atom.position - xyz) < displacement
