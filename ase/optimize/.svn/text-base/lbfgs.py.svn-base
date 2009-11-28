import numpy as np
import copy 
from ase.optimize import Optimizer
from ase.neb import *
import random

class _LBFGS(Optimizer):
    """Limited memory bfgs algorithm. Unlike the bfgs algorithm used in qn.py,
       the inverse of hessian matrix is updated. 

    Parameters:

    restart: string
        Pickle file used to store vectors for updating the inverse of
        hessian matrix. If set, file with such a name will be searched
        and information stored will be used, if the file exists.

    memory: int
        Number of steps to be stored. Default value is 25.

    method: string
        Two methods for determing atomic movement are available. If
        method = 'line', a line search will be performed to determine
        the atomic movement. An extra scf loop for the trail step in
        each atomic step. And if method = 'hess', the atomic step will
        be determined by hessian matrix, which means each atomic step
        include only one scf loop.  """

    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=0.2, dR=0.1,
                 memory=25, alpha=0.05, method=None, damping=1.):
        Optimizer.__init__(self, atoms, restart, logfile, trajectory)

        if maxstep > 1.0:
            raise ValueError(
                             'Wanna fly? I know the calculation is too slow. ' +
                             'But you have to follow the rules.\n'+
                             '            The maximum step size %.1f' % maxstep 
                             +' is too big! \n'+
                             '            Try to set the maximum step size'+
                             ' below 0.2.')
        self.maxstep = maxstep
        self.dR = dR
        self.memory = memory
        self.alpha = alpha
        self.method = method
        self.damping = damping

    def initialize(self):
        self.ITR = None
        self.f_old = None
        self.r_old = None

    def sign(self,w):
        if(w < 0.0):
            return -1.0
        return 1.0

    def read(self):
        (self.ITR, self.s, self.y, self.rho, self.r_old, 
         self.f_old) = self.load()

    def step(self, f):
        r = self.atoms.get_positions()
        self.update(r, f, self.r_old, self.f_old)

        du = self.d / np.sqrt(np.vdot(self.d, self.d))
        if(self.method == 'hess'): 
        # use the Hessian Matrix to predict the min
            dr = self.d
            if(abs(np.sqrt(np.vdot(dr, dr).sum())) > self.maxstep):
                dr = du * self.maxstep

        elif(self.method == 'line'):
        # Finite difference step using temporary point
            tmp_r = r.copy()
            tmp_r += (du * self.dR)
            self.tmp.set_positions(tmp_r)
        # Decide how much to move along the line du
            Fp1 = np.vdot(f, du)
            Fp2 = np.vdot(self.tmp.get_forces(), du)
            CR = (Fp1 - Fp2) / self.dR
            #RdR = Fp1*0.1
            if(CR < 0.0):
                print "negcurve"
                RdR = self.maxstep
                #if(abs(RdR) > self.maxstep):
                #    RdR = self.sign(RdR) * self.maxstep
            else:
                Fp = (Fp1 + Fp2) * 0.5
                RdR = Fp / CR 
                if(abs(RdR) > self.maxstep):
                    RdR = self.sign(RdR) * self.maxstep
                else:
                    RdR += self.dR * 0.5
            dr = du * RdR
        self.r_old = r.copy()
        self.f_old = f.copy()
        r += dr * self.damping
        self.atoms.set_positions(r)

    def update(self, r, f, r_old, f_old):
        a = np.zeros(self.memory + 1, 'd')
        self.tmp = self.atoms
        self.Ho = np.ones((np.shape(r)[0], 3), 'd')
        if (self.method == 'hess'):
            self.Ho = self.Ho * self.alpha
        if(not self.ITR):
            self.ITR = 1
            self.s = [1.] # The 0'th element is not actually used
            # The point is to use 1-indexation
            self.y = [1.]
            self.rho = [1.]
        else:
            a1 = abs (np.vdot(f, f_old))
            a2 = np.vdot(f_old, f_old)
            if(self.method == 'line'):
                if(a1 <= 0.5 * a2 and a2 != 0):
                    reset_flag = 0
                else:
                    reset_flag = 1
            else:
                reset_flag = 0
            if(reset_flag == 0):
                ITR = self.ITR
                if(ITR > self.memory):
                    self.s.pop(1)
                    self.y.pop(1)
                    self.rho.pop(1)
                    ITR = self.memory
                self.s.append(r - r_old)
                self.y.append(-(f - f_old))
                self.rho.append(1 / np.vdot(self.y[ITR],self.s[ITR]))
                self.ITR += 1
            else:
                self.ITR = 1
                self.s = [1.]
                self.y = [1.]
                self.rho = [1.]
        self.dump((self.ITR, self.s, self.y, self.rho, r_old, f_old))

        r_old = r.copy()
        f_old = f.copy()
        if(self.ITR <= self.memory):
            BOUND = self.ITR
        else:
            BOUND = self.memory
        q = -1.0 * f
        for j in range(1,BOUND):
            k = (BOUND - j)
            a[k] = self.rho[k] * np.vdot(self.s[k], q)
            q -= a[k] * self.y[k]
        d = self.Ho * q 
        for j in range(1,BOUND):
            B = self.rho[j] * np.vdot(self.y[j], d)
            d= d + self.s[j] * (a[j] - B)
        self.d = -1.0 * d

    def replay_trajectory(self, traj):
        """Initialize hessian from old trajectory."""
        if isinstance(traj, str):
            from ase.io.trajectory import PickleTrajectory
            traj = PickleTrajectory(traj, 'r')
        atoms = traj[0]
        r_old = atoms.get_positions()
        f_old = atoms.get_forces()
        for i in range(0,len(traj)-1):
            r = traj[i].get_positions()
            f = traj[i].get_forces()
            self.update(r, f, r_old, f_old)
            r_old = r
            f_old = f
        self.r_old = traj[-2].get_positions()
        self.f_old = traj[-2].get_forces()

class LineLBFGS(_LBFGS):
    def __init__(self, *args, **kwargs):
        kwargs['method'] = 'line'
        _LBFGS.__init__(self, *args, **kwargs)

class HessLBFGS(_LBFGS):
    def __init__(self, *args, **kwargs):
        kwargs['method'] = 'hess'
        _LBFGS.__init__(self, *args, **kwargs)

