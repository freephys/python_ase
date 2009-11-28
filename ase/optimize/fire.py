import numpy as np

from ase.optimize import Optimizer


class FIRE(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 dt=0.1,
                 maxmove=0.2,dtmax=1.0,Nmin=5,finc=1.1,fdec=0.5,
                 astart=0.1,fa=0.99,a=0.1):
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

      #  if dt is not None:
        self.dt = dt

        self.Nsteps = 0
        self.maxmove = maxmove
        self.dtmax = dtmax
        self.Nmin = Nmin
        self.finc = finc
        self.fdec = fdec
        self.astart = astart
        self.fa = fa
        self.a = a
	
    def initialize(self):
        self.v = None
#        self.dt = 0.1

    def read(self):
        self.v, self.dt = self.load()
       
    def step(self,f):
#	print self.dt
        atoms = self.atoms
        if self.v is None:
            self.v = np.zeros((len(atoms), 3))
        else:
            vf = np.vdot(f,self.v)
#            self.v = (1.0-self.a)*self.v+self.a*f/np.sqrt(np.vdot(f,f)*np.vdot(self.v,self.v))
            if vf > 0.0:
                self.v = (1.0-self.a)*self.v + self.a*f / np.sqrt(np.vdot(f,f)) * np.sqrt(np.vdot(self.v,self.v))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt*self.finc,self.dtmax)
                    self.a *= self.fa
                self.Nsteps += 1
            else:
                self.v[:] *= 0.0
                self.a = self.astart
                self.dt *= self.fdec
                self.Nsteps = 0
#            if vf < 0.0:
#                self.v[:] = 0.0
#                self.a = self.astart
#                self.dt *= self.fdec
#                self.Nsteps = 0
#            else:
#                self.v = (1.0-self.a)*self.v+self.a*f*np.sqrt(np.vdot(f,f)/np.vdot(self.v,self.v))
#                if self.Nsteps > self.Nmin:
#                    dt = min(dt*self.finc,dtmax)
#                    self.a *= self.fa
#                    self.Nsteps += 1

            self.v += self.dt * f
#            self.v += 0.5*self.dt * f
            dr = self.dt*self.v
            if np.sqrt(np.vdot(dr,dr))>self.maxmove:
                dr = self.maxmove*dr/np.sqrt(np.vdot(dr,dr))
            r = atoms.get_positions()
            atoms.set_positions(r + dr)
            self.dump((self.v, self.dt))
