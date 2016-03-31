The PyDy source code is actually spread across both the SymPy and PyDy
repositories. The two projects are very intimately tied together in terms of
development. Much of our development practices come from SymPy's methods so
reading their documentation for contributing is the best place to start:

https://github.com/sympy/sympy/wiki/introduction-to-contributing

Pull Request Check Lists
========================

We do have some differences though. In particular, we use a checklists to make
sure each pull request is ready to be merged into the master branch (which
should always be stable). The following gives some example checklists for
various types of pull requests:

New Feature
-----------

::

   - [ ] There are no merge conflicts.
   - [ ] If there is a related issue, a reference to that issue is in the
     commit message.
   - [ ] Unit tests have been added for the new feature.
   - [ ] The PR passes tests both locally (run `nosetests`) and on Travis CI.
   - [ ] All public methods and classes have docstrings. (We use the [numpydoc
     format](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).)
   - [ ] An explanation has been added to the online documentation. (`docs`
     directory)
   - [ ] The code follows PEP8 guidelines. (use a linter, e.g.
     [pylint](http://www.pylint.org), to check your code)
   - [ ] The new feature is documented in the [Release
     Notes](https://github.com/pydy/pydy#release-notes).
   - [ ] The code is backwards compatible. (All public methods/classes must
     follow deprecation cycles.)
   - [ ] All reviewer comments have been addressed.

Bug Fix
-------

::

   - [ ] There are no merge conflicts.
   - [ ] If there is a related issue, a reference to that issue is in the
     commit message.
   - [ ] Unit tests have been added for the bug. (Please reference the issue #
     in the unit test.)
   - [ ] The tests pass both locally (run `nosetests`) and on Travis CI.
   - [ ] The code follows PEP8 guidelines. (use a linter, e.g.
     [pylint](http://www.pylint.org), to check your code)
   - [ ] The bug fix is documented in the [Release
     Notes](https://github.com/pydy/pydy#release-notes).
   - [ ] The code is backwards compatible. (All public methods/classes must
     follow deprecation cycles.)
   - [ ] All reviewer comments have been addressed.

Example
-------

::

   - [ ] There are no merge conflicts.
   - [ ] A figure (preferrably SVG) has been added for the example.
   - [ ] The required versions of all dependencie have been noted. (For
     notebooks, the [version information
     extension](https://github.com/jrjohansson/version_information) is helpful
     for this.)
   - [ ] A `run.py` file that executes the example is included.
   - [ ] All reviewer comments have been addressed.

The appropriate checklist should be pasted in the top most comment block of the
pull request. It is then the duty of reviewers to check off the items as they
are completed. If any questions are not relevant to your particular pull
request, a reviewer will simply check it off as done.

Style
=====

Code written in python should roughly adhere to PEP8:

http://www.python.org/dev/peps/pep-0008/.

There are some exceptions, such as naming of variables. Use your best judgment.

Commit Guidelines
=================

When submitting, please follow the commit guidelines used for Git:

http://git-scm.com/book/ch5-2.html#Commit-Guidelines

The first line should have a soft limit of 50 characters, and subsequent lines
limitted to 72 characters. The commit messages should be written in the
imperative present tense.
