from ase import *

a = 5.0
d = 1.9
c = a / 2
atoms = Atoms('AuH',
              positions=[(c, c, 0), (c, c, d)],
              cell=(a, a, 2 * d),
              pbc=(0, 0, 1))
atoms *= (1, 1, 2)
images = [atoms.copy(), atoms.copy()]

r = ['xyz', 'traj', 'cube', 'pdb']
w = r + ['xsf']
try:
    import matplotlib
except ImportError:
    pass
else:
    w += ['png', 'eps']
    
for format in w:
    print format, 'O',
    fname1 = 'io-test.1.' + format
    fname2 = 'io-test.2.' + format
    write(fname1, atoms, format=format)
    if format not in ['cube', 'png', 'eps']:
        write(fname2, images, format=format)

    if format in r:
        print 'I'
        a1 = read(fname1)
        if format != 'cube':
            a2 = read(fname2)
            a3 = read(fname2, index=0)
            a4 = read(fname2, index=slice(None))
    else:
        print
