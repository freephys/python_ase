.. _overview:

========
Overview
========

ASE is an Atomistic Simulation Environment written in the
Python_ programming language with the aim of setting up, stearing, and
analyzing atomistic simulations. The ASE has been constructed with a
number of "design goals" as:


:Simplicity in use:
  Setting up an atomistic total energy calculation or molecular
  dynamics simulation with ASE is simple and straightforward. The Python
  scripts are almost self-explanatory
  (see :ref:`python_info` for a short introduction)
  and look like well
  commented input files to an atomistic simulation program.

:Flexibility in use:
  Since ASE is based on the Python scripting language it is possible
  without any code modifications to perform very complicated simulation
  tasks. For example a sequence of calculations may be performed with
  the use of simple "for-loop" constructions or simulations of different
  types (:term:`DFT` and classical molecular mechanics potentials) may
  be coupled together.

:Simplicity in development:
  The ASE defines a set of interfaces for different objects, i.e. an
  :class:`~ase.atoms.Atoms` object is required to posses a method with the name
  :meth:`~ase.atoms.Atoms.get_positions` which returns the coordinates of
  the atom. By following a few such standard interfaces it is easy for
  new users to get access to all of the functionality of ASE.

:Flexibility in development:
  The Python code in ASE is structured in different modules intended for
  different purposes. There are :mod:`calculators` for calculating
  energies, forces and stresses, :mod:`md` and :mod:`optimize` modules
  for controlling the motion of atoms, :mod:`constraint <constraints>`
  objects and filters for performing :mod:`nudged-elastic-band <neb>`
  calculations etc. The modularity of the code and the documented
  interfaces make it simple to contribute new functionality to ASE.

:Pythonic:
  It fits nicely into the rest of the Python world with heavy
  use of the popular :term:`NumPy` package for numerical work
  (see :ref:`numpy` for a short introduction). The
  use of the Python language allows ASE to be used both interactively
  as well as in scripts.

:Open to participation:
  The CAMPOS Atomic Simulation Environment is released under the GNU
  Lesser General Public License.  See the files :trac:`COPYING` and
  :trac:`COPYING.LESSER` which accompany the downloaded files, or see
  the license at GNU's web server at
  http://www.gnu.org/copyleft/lgpl.html.  Everybody is invited to
  participate in using and :ref:`developing the code <devel>`.

.. _Python: http://www.python.org
