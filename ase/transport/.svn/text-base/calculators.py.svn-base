import numpy as np

from numpy import linalg
from ase.transport.selfenergy import LeadSelfEnergy, BoxProbe
from ase.transport.greenfunction import GreenFunction
from ase.transport.tools import subdiagonalize, cutcoupling, tri2full, dagger


class TransportCalculator:
    """Determine transport properties of device sandwiched between
    semi-infinite leads using nonequillibrium Green function methods.
    """

    def __init__(self, **kwargs):
        """Bla Bla XXX
        
        energies is the energy grid on which the transport properties
        should be determined.
        
        h1 (h2) is a matrix representation of the Hamiltonian of two
        principal layers of the left (right) lead, and the coupling
        between such layers.
        
        h is a matrix representation of the Hamiltonian of the
        scattering region. This must include at least one lead
        principal layer on each side. The coupling in (out) of the
        scattering region is by default assumed to be identical to the
        coupling between left (right) principal layers.  However,
        these couplings can also be specified explicitly through hc1
        and hc2.
        
        s, s1, and s2 are the overlap matrices corresponding to h, h1,
        and h2. Default is the identity operator. sc1 and sc2 are the
        overlap matrices corresponding to the optional couplings hc1
        and hc2.
        
        align_bf specifies the principal layer basis index used to
        align the fermi levels of the lead and scattering regions.
        """
        
        # The default values for all extra keywords
        self.input_parameters = {'energies': None,
                                 'h': None,
                                 'h1': None,
                                 'h2': None,
                                 's': None,
                                 's1': None,
                                 's2': None,
                                 'hc1': None,
                                 'hc2': None,
                                 'sc1': None,
                                 'sc2': None,
                                 'box': None,
                                 'align_bf': None,
                                 'eta1': 1e-3,
                                 'eta2': 1e-3,
                                 'eta': 1e-3,
                                 'logfile': None, # '-',
                                 'eigenchannels': 0,
                                 'dos': False,
                                 'pdos': [],
                                 }
        self.initialized = False # Changed Hamiltonians?
        self.uptodate = False # Changed energy grid?
        self.set(**kwargs)

    def set(self, **kwargs):
        for key in kwargs:
            if key in ['h', 'h1', 'h2', 'hc1', 'hc2',
                       's', 's1', 's2', 'sc1', 'sc2',
                       'eta', 'eta1', 'eta2', 'align_bf', 'box']:
                self.initialized = False
                self.uptodate = False
                break
            elif key in ['energies', 'eigenchannels', 'dos', 'pdos']:
                self.uptodate = False
            elif key not in self.input_parameters:
                raise KeyError, '\'%s\' not a vaild keyword' % key

        self.input_parameters.update(kwargs)
        log = self.input_parameters['logfile']
        if log is None:
            class Trash:
                def write(self, s):
                    pass
                def flush(self):
                    pass
            self.log = Trash()
        elif log == '-':
            from sys import stdout
            self.log = stdout
        elif 'logfile' in kwargs:
            self.log = open(log, 'w')

    def initialize(self):
        if self.initialized:
            return

        print >> self.log, '# Initializing calculator...'

        p = self.input_parameters
        if p['s1'] == None:
            p['s1'] = np.identity(len(p['h1']))
        if p['s2'] == None:
            p['s2'] = np.identity(len(p['h2']))
        if p['s'] == None:
            p['s'] = np.identity(len(p['h']))
            
        h_mm = p['h']
        s_mm = p['s']
        pl1 = len(p['h1']) / 2
        pl2 = len(p['h2']) / 2
        h1_ii = p['h1'][:pl1, :pl1]
        h1_ij = p['h1'][:pl1, pl1:2 * pl1]
        s1_ii = p['s1'][:pl1, :pl1]
        s1_ij = p['s1'][:pl1, pl1:2 * pl1]
        h2_ii = p['h2'][:pl2, :pl2]
        h2_ij = p['h2'][pl2: 2 * pl2, :pl2]
        s2_ii = p['s2'][:pl2, :pl2]
        s2_ij = p['s2'][pl2: 2 * pl2, :pl2]
        
        if p['hc1'] is None:
            nbf = len(h_mm)
            h1_im = np.zeros((pl1, nbf), complex)
            s1_im = np.zeros((pl1, nbf), complex)
            h1_im[:pl1, :pl1] = h1_ij
            s1_im[:pl1, :pl1] = s1_ij
        else:
            h1_im = p['hc1']
            if p['sc1'] is not None:
                s1_im = p['sc1']
            else:
                s1_im = np.zeros(h1_im.shape, complex)

        if p['hc2'] is None:
            h2_im = np.zeros((pl2, nbf), complex)
            s2_im = np.zeros((pl2, nbf), complex)
            h2_im[-pl2:, -pl2:] = h2_ij
            s2_im[-pl2:, -pl2:] = s2_ij
        else:
            h2_im = p['hc2']
            if p['sc2'] is not None:
                s2_im[:] = p['sc2']
            else:
                s2_im = np.zeros(h2_im.shape, complex)

        align_bf = p['align_bf']
        if align_bf != None:
            diff = (h_mm[align_bf, align_bf] - h1_ii[align_bf, align_bf]) \
                   / s_mm[align_bf, align_bf]
            print >> self.log, '# Aligning scat. H to left lead H. diff=', diff
            h_mm -= diff * s_mm

        #setup lead self-energies
        self.selfenergies = [LeadSelfEnergy((h1_ii, s1_ii), 
                                            (h1_ij, s1_ij),
                                            (h1_im, s1_im),
                                            p['eta1']),
                             LeadSelfEnergy((h2_ii, s2_ii), 
                                            (h2_ij, s2_ij),
                                            (h2_im, s2_im),
                                            p['eta2'])]
        box = p['box']
        if box is not None:
            print 'Using box probe!'
            self.selfenergies.append(
                BoxProbe(eta=box[0], a=box[1], b=box[2], energies=box[3],
                         S=s_mm, T=0.3))
        
        #setup scattering green function
        self.greenfunction = GreenFunction(selfenergies=self.selfenergies,
                                           H=h_mm,
                                           S=s_mm,
                                           eta=p['eta'])

        self.initialized = True
    
    def update(self):
        if self.uptodate:
            return
        
        p = self.input_parameters
        self.energies = p['energies']
        nepts = len(self.energies)
        nchan = p['eigenchannels']
        pdos = p['pdos']
        self.T_e = np.empty(nepts)
        if p['dos']:
            self.dos_e = np.empty(nepts)
        if pdos != []:
            self.pdos_ne = np.empty((len(pdos), nepts))
        if nchan > 0:
            self.eigenchannels_ne = np.empty((nchan, nepts))

        for e, energy in enumerate(self.energies):
            Ginv_mm = self.greenfunction.retarded(energy, inverse=True)
            lambda1_mm = self.selfenergies[0].get_lambda(energy)
            lambda2_mm = self.selfenergies[1].get_lambda(energy)
            a_mm = linalg.solve(Ginv_mm, lambda1_mm)
            b_mm = linalg.solve(dagger(Ginv_mm), lambda2_mm)
            T_mm = np.dot(a_mm, b_mm)
            if nchan > 0:
                t_n = linalg.eigvals(T_mm).real
                self.eigenchannels_ne[:, e] = np.sort(t_n)[-nchan:]
                self.T_e[e] = np.sum(t_n)
            else:
                self.T_e[e] = np.trace(T_mm).real

            print >> self.log, energy, self.T_e[e]
            self.log.flush()

            if p['dos']:
                self.dos_e[e] = self.greenfunction.dos(energy)

            if pdos != []:
                self.pdos_ne[:, e] = np.take(self.greenfunction.pdos(energy),
                                             pdos)
        
        self.uptodate = True

    def print_pl_convergence(self):
        self.initialize()
        pl1 = len(self.input_parameters['h1']) / 2
        
        h_ii = self.selfenergies[0].h_ii
        s_ii = self.selfenergies[0].s_ii
        ha_ii = self.greenfunction.H[:pl1, :pl1]
        sa_ii = self.greenfunction.S[:pl1, :pl1]
        c1 = np.abs(h_ii - ha_ii).max()
        c2 = np.abs(s_ii - sa_ii).max()
        print 'Conv (h,s)=%.2e, %2.e' % (c1, c2)

    def plot_pl_convergence(self):
        self.initialize()
        pl1 = len(self.input_parameters['h1']) / 2       
        hlead = self.selfenergies[0].h_ii.real.diagonal()
        hprincipal = self.greenfunction.H.real.diagonal[:pl1]

        import pylab as pl
        pl.plot(hlead, label='lead')
        pl.plot(hprincipal, label='principal layer')
        pl.axis('tight')
        pl.show()

    def get_transmission(self):
        self.initialize()
        self.update()
        return self.T_e

    def get_dos(self):
        self.initialize()
        self.update()
        return self.dos_e

    def get_eigenchannels(self, n=None):
        """Get ``n`` first eigenchannels."""
        self.initialize()
        self.update()
        if n is None:
            n = self.input_parameters['eigenchannels']
        return self.eigenchannels_ne[:n]

    def get_pdos(self):
        self.initialize()
        self.update()
        return self.pdos_ne

    def subdiagonalize_bfs(self, bfs):
        self.initialize()
        bfs = np.array(bfs)
        p = self.input_parameters
        h_pp = p['h']
        s_pp = p['s']
        ht_pp, st_pp, c_pp, e_p = subdiagonalize(h_pp, s_pp, bfs)
        c_pp = np.take(c_pp, bfs, axis=0)
        c_pp = np.take(c_pp, bfs, axis=1)
        return ht_pp, st_pp, e_p, c_pp

    def cutcoupling_bfs(self, bfs):
        self.initialize()
        bfs = np.array(bfs)
        p = self.input_parameters
        h_pp = p['h'].copy()
        s_pp = p['s'].copy()
        cutcoupling(h_pp, s_pp, bfs)
        return h_pp, s_pp
        
    def get_left_channels(self, energy, nchan=1):
        self.initialize()
        g_s_ii = self.greenfunction.retarded(energy)
        lambda_l_ii = self.selfenergies[0].get_lambda(energy)
        lambda_r_ii = self.selfenergies[1].get_lambda(energy)

        if self.greenfunction.S is None:
            s_s_qsrt_ii = s_s_isqrt = np.identity(len(g_s_ii))
        else:
            s_mm = self.greenfunction.S
            s_s_i, s_s_ii = linalg.eig(s_mm)
            s_s_i = np.abs(s_s_i)
            s_s_sqrt_i = np.sqrt(s_s_i) # sqrt of eigenvalues  
            s_s_sqrt_ii = np.dot(s_s_ii * s_s_sqrt_i, dagger(s_s_ii))
            s_s_isqrt_ii = np.dot(s_s_ii / s_s_sqrt_i, dagger(s_s_ii))

        lambdab_r_ii = np.dot(np.dot(s_s_isqrt_ii, lambda_r_ii),s_s_isqrt_ii)
        a_l_ii = np.dot(np.dot(g_s_ii, lambda_l_ii), dagger(g_s_ii))
        ab_l_ii = np.dot(np.dot(s_s_sqrt_ii, a_l_ii), s_s_sqrt_ii)
        lambda_i, u_ii = linalg.eig(ab_l_ii)
        ut_ii = np.sqrt(lambda_i / (2.0 * np.pi)) * u_ii
        m_ii = 2 * np.pi * np.dot(np.dot(dagger(ut_ii), lambdab_r_ii),ut_ii)
        T_i,c_in = linalg.eig(m_ii)
        T_i = np.abs(T_i)
        
        channels = np.argsort(-T_i)[:nchan]
        c_in = np.take(c_in, channels, axis=1)
        T_n = np.take(T_i, channels)
        v_in = np.dot(np.dot(s_s_isqrt_ii, ut_ii), c_in)

        return T_n, v_in
