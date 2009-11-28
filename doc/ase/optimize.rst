.. _structure_optimizations:

======================
Structure optimization
======================
.. module:: optimize
   :synopsis: Structure Optimization

There are currently 4 different optimization algorithms available:
``QuasiNewton``, ``MDMin``, ``FIRE``, and ``LBFGS``.

``MDMin`` and ``FIRE`` both use Newtonian dynamics with added
friction, to converge to an energy minimum, whereas ``QuasiNewton``
uses the forces of consecutive steps to dynamically update a Hessian
describing the curvature of the potential energy landscape. ``GLBFGS``
is an experimental optimizer designed for simultaneous update of the
images along a nudged elastic band trajectory.

All optimizer classes have the following structure::

  class Optimizer:
      def __init__(self, atoms, restart=None, logfile=None):
      def run(self, fmax=0.05, steps=100000000):
      def get_number_of_steps():

The convergence criterion is that the force on all individual atoms
should be less than *fmax*:

.. math:: \max_a |\vec{F_a}| < fmax


QuasiNewton
-----------
.. module:: optimize.qn
   :synopsis: Quasi-Newton

The ``QuasiNewton`` object is one of the minimizers in the ASE
package.  Let's try to use it to optimize the structure of a water
molecule.  We start with the experimental geometry::

  from ase import *
  d = 0.9575
  t = pi / 180 * 104.51
  water = Atoms('H2O',
                positions=[(d, 0, 0),
                           (d * cos(t), d * sin(t), 0),
                           (0, 0, 0)],
                calculator=EMT())
  dyn = QuasiNewton(water)
  dyn.run(fmax=0.05)
  QuasiNewton:   0        6.445801      51.6847
  QuasiNewton:   1        2.418583      27.2946
  QuasiNewton:   2        0.551767      12.1607
  QuasiNewton:   3       -0.039301       4.0520
  QuasiNewton:   4       -0.128045       0.8479
  QuasiNewton:   5       -0.132312       0.0397

When doing structure optimization, it is useful to write the
trajectory to a file, so that the progress of the optimization run can
be followed during or after the run::

  dyn = QuasiNewton(water, trajectory='H2O.traj')
  dyn.run(fmax=0.05)
  
Use the command ``ag H2O.traj`` to see what is going on (more here:
:mod:`gui`).  The trajectory file can also be accessed using the
module :mod:`ase.io.trajectory`.

The ``attach`` method takes an optional argument ``interval=n`` that can
be used to tell the structure optimizer object to write the
configuration to the trajectory file only every ``n`` steps.

During a structure optimization, the :class:`QuasiNewton` and
:class:`LBFGS` optimizers use two quantities to decide where to move
the atoms on each step:

 * the forces on each atom, as returned by the associated :class:`Calculator`
   object
 * the Hessian matrix, i.e. the matrix of second derivatives
   :math:`\frac{\partial^2 E}{\partial x_i \partial x_j}` of the
   total energy with respect to nuclear coordinates.

If the atoms are close to the minimum, such that the potential energy
surface is locally quadratic, the Hessian and forces accurately
determine the required step to reach the optimal structure.  The
Hessian is very expensive to calculate *a priori*, so instead the
algorithm estimates it by means of an initial guess which is adjusted
along the way depending on the information obtained on each step of
the structure optimization.

It is frequently practical to restart or continue a structure
optimization with a geometry obtained from a previous relaxation.
Aside from the geometry, the Hessian of the previous run can and
should be retained for the second run.  Use the ``restart`` keyword to
specify a file in which to save the Hessian::

  dyn = QuasiNewton(system, trajectory='qn.traj', restart='qn.pckl')

This will create an optimizer which saves the Hessian to
:file:`qn.pckl` (using the Python :mod:`pickle` module) on each
step.  If the file already exists, the Hessian will also be
*initialized* from that file.

The trajectory file can also be used to restart a structure
optimization, since it contains the history of all forces and
positions, and thus whichever information about the Hessian was
assembled so far::

  dyn = QuasiNewton(system, trajectory='qn.traj')
  dyn.replay_trajectory('history.traj')

This will read through each iteration stored in :file:`history.traj`,
performing adjustments to the Hessian as appropriate.  Note that these
steps will not be written to :file:`qn.traj`.  If restarting with more than
one previous trajectory file, use :command:`ag` to concatenate them
into a single trajectory file first::

  $ ag part1.traj part2.traj -o history.traj

The file :file:`history.traj` will then contain all necessary
information.

When switching between different types of optimizers, e.g. between
``QuasiNewton`` and ``LBFGS``, the pickle-files specified by the
``restart`` keyword are not compatible, but the Hessian can still be
retained by replaying the trajectory as above.

LBFGS
-----
.. module:: optimize.lbfgs

The LBFGS is the limited memory version of BFGS algorithm, where 
the inverse of Hessian matrix is updated instead of the Hessian
itself. Tow ways exist in the LBFGS class for determining atomic
step, which are called ``HessLBFGS`` and ``LineLBFGS``. For the 
first one, both the directions and lengths of the atomic steps 
are determined by the approximated Hessian matrix. While for the 
latter one, the approximated Hessian matrix is only used to find 
out the directions of the line searches and atomic steps, the 
step lengths are determined by the forces. 

To start a structure optimization with LBFGS algorithm is similar to
QuasiNewton. A typical optimization should look like::

  dyn = HessLBFGS(system, trajectory='lbfgs.traj', restart='lbfgs.pckl')

where the trajectory and the restart save the trajectory of the 
optimization and the vectors needed to generate the Hessian Matrix.

FIRE
----
.. module:: optimize.fire

...

MDMin
-----
.. module:: optimize.mdmin

The MDmin algorithm is a modification of the usual velocity-Verlet
molecular dynamics algorithm.  Newtons second law is solved
numerically, but after each time step the dot product between the
forces and the momenta is checked.  If it is zero, the system has just
passed through a (local) minimum in the potential energy, the kinetic
energy is large and about to decrease again.  At this point, the
momentum is set to zero.  Unlike a "real" molecular dynamics, the
masses of the atoms are not used, instead all masses are set to one.

The MDmin algorithm exists in two flavors, one where each atom is
tested and stopped individually, and one where all coordinates are
treated as one long vector, and all momenta are set to zero if the
dotproduct between the momentum vector and force vector (both of
length 3N) is zero.  This module implements the latter version.

Although the algorithm is primitive, it performs very well because it
takes advantage of the physics of the problem.  Once the system is so
near the minimum that the potential energy surface is approximately
quadratic it becomes advantageous to switch to a minimization method
with quadratic convergence, such as `Conjugate Gradient` or `Quasi
Newton`.
