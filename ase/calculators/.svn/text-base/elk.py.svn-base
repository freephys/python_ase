import os

from ase.units import Bohr, Hartree

elk_parameters = {
    'swidth': Hartree,
    }

class ELK:
    def __init__(self, dir='.', xc=None, tasks=[0], **kwargs):
        for key, value in kwargs.items():
            if key in elk_parameters:
                kwargs[key] /= elk_parameters[key]

        if xc is not None:
            if 'xctype' in kwargs:
                raise ValueError("You can't use both 'xc' and 'xctype'!")
            else:
                kwargs['xctype'] = {'LDA': 3, # PW92
                                    'PBE': 20,
                                    'REVPBE': 21,
                                    'PBESOL': 22,
                                    'WC06': 26,
                                    'AM05': 30}[xc.upper()]

        kwargs['tasks'] = tasks

        self.parameters = kwargs

        self.dir = dir
        self.energy = None

    def get_potential_energy(self, atoms):
        self.write(atoms)
        assert os.system('cd %s; elk' % self.dir) == 0
        self.read()
        return self.energy
    
    def write(self, atoms):
        if not os.path.isdir(self.dir):
            os.mkdir(self.dir)
        fd = open('%s/elk.in' % self.dir, 'w')
        for key, value in self.parameters.items():
            fd.write('%s\n' % key)
            if isinstance(value, bool):
                fd.write('.%s.\n\n' % ('false', 'true')[value])
            elif isinstance(value, (int, float)):
                fd.write('%s\n\n' % value)
            else:
                fd.write('%s\n\n' % ' '.join([str(x) for x in value]))

        fd.write('avec\n')
        for vec in atoms.cell:
            fd.write('%.14f %.14f %.14f\n' % tuple(vec / Bohr))
        fd.write('\n')

        fd.write("sppath\n'%s'\n\n" % os.environ['ELK_SPECIES_PATH'])

        species = {}
        symbols = []
        for a, symbol in enumerate(atoms.get_chemical_symbols()):
            if symbol in species:
                species[symbol].append(a)
            else:
                species[symbol] = [a]
                symbols.append(symbol)
        fd.write('atoms\n%d\n' % len(species))
        scaled = atoms.get_scaled_positions()
        for symbol in symbols:
            fd.write("'%s.in'\n" % symbol)
            fd.write('%d\n' % len(species[symbol]))
            for a in species[symbol]:
                fd.write('%.14f %.14f %.14f 0.0 0.0 0.0\n' % tuple(scaled[a]))

    def read(self):
        fd = open('%s/TOTENERGY.OUT' % self.dir, 'r')
        self.energy = float(fd.readlines()[-1]) * Hartree
