.. _faq:

==========================
Frequently Asked Questions
==========================

General
=======

Citation: how should I cite ASE?
--------------------------------

If you find ASE useful in your research please cite:

   | S. R. Bahn and K. W. Jacobsen
   | `An object-oriented scripting interface to a legacy electronic structure code`__
   | Comput. Sci. Eng., Vol. **4**, 56-66, 2002

   __ http://dx.doi.org/10.1109/5992.998641

BibTex (:svn:`doc/ASE.bib`):

.. literalinclude:: ASE.bib

Download
========

Trying to checkout the code via SVN resulted::

 [~]$ svn checkout "https://svn.fysik.dtu.dk/projects/ase/trunk"
 svn: Unrecognized URL scheme 'https://svn.fysik.dtu.dk/projects/ase/trunk'

This error is diplayed in case the library 'libsvn_ra_dav' is missing on your system. The library is used by SVN, but is not installed by default. 
