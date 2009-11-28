.. _python_codingstandard:

==================
Coding Conventions
==================

Importing modules
=================

We distinguish between scripts and code.  In scripts, like the ones
found in the examples and tutorials, we use::

  from ase import *

which gives us the most used symbols.

I code, like the implementation of ASE, we must *not* use the
``import *`` syntax.  Import everthing explicitly from exactly the
place where it's defined::

  from ase.io import read, write

Python Coding Conventions
=========================

The rules for the Python part are almost identical
to those used by the `Docutils project`_:

Contributed code will not be refused merely because it does not
strictly adhere to these conditions; as long as it's internally
consistent, clean, and correct, it probably will be accepted.  But
don't be surprised if the "offending" code gets fiddled over time to
conform to these conventions.

The project shall follow the generic coding conventions as
specified in the `Style Guide for Python Code`_ and `Docstring
Conventions`_ PEPs, summarized, clarified, and extended as follows:

* 4 spaces per indentation level.  No hard tabs.

* Avoid introducing `trailing whitespaces <http://www.gnu.org/software/emacs/manual/html_node/emacs/Useless-Whitespace.html>`_.

* Try to use only 7-bit ASCII, no 8-bit strings.

* No one-liner compound statements (i.e., no ``if x: return``: use two
  lines & indentation), except for degenerate class or method
  definitions (i.e., ``class X: pass`` is OK.).

* Lines should be no more than 78 characters long.

* Use "StudlyCaps" for class names.

* Use "lowercase" or "lowercase_with_underscores" for function,
  method, and variable names.  For short names, maximum two words,
  joined lowercase may be used (e.g. "tagname").  For long names with
  three or more words, or where it's hard to parse the split between
  two words, use lowercase_with_underscores (e.g.,
  "note_explicit_target", "explicit_target").  If in doubt, use
  underscores.

* Avoid lambda expressions, which are inherently difficult to
  understand.  Named functions are preferable and superior: they're
  faster (no run-time compilation), and well-chosen names serve to
  document and aid understanding.

* Avoid functional constructs (filter, map, etc.).  Use list
  comprehensions instead.

* Avoid ``from __future__ import`` constructs.  They are inappropriate
  for production code.

* Use 'single quotes' for string literals, and """triple double
  quotes""" for :term:`docstring`\ s.  Double quotes are OK for
  something like ``"don't"``.

.. _Style Guide for Python Code: http://www.python.org/peps/pep-0008.html
.. _Docstring Conventions: http://www.python.org/peps/pep-0257.html
.. _Docutils project: http://docutils.sourceforge.net/docs/dev/policies.html#python-coding-conventions

.. attention::

   Thus spake the Lord: Thou shalt indent with four spaces. No more, no less.
   Four shall be the number of spaces thou shalt indent, and the number of thy
   indenting shall be four. Eight shalt thou not indent, nor either indent thou
   two, excepting that thou then proceed to four. Tabs are right out.

                                          Georg Brandl

General advice
==============

 * Get rid of as many ``break`` and ``continue`` statements as possible.

Writing documentation in the code
=================================

Here is an example of how to write good docstrings:

  http://projects.scipy.org/numpy/browser/trunk/doc/example.py
