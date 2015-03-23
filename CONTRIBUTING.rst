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

  - [ ] Are there merge conflicts?
  - [ ] If there is an issue for this feature is it referenced in the commit
    message?
  - [ ] Have you added unit tests?
  - [ ] Did it pass the tests? Locally (run `nosetests`)? Travis CI?
  - [ ] Do all public methods/classes have docstrings? (We use the [numpydoc
    format](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).)
  - [ ] Did you add an explanation to the online documentation? (`docs`
    directory)
  - [ ] Does it follow PEP8 guidelines? (use a linter, e.g.
    [pylint](http://www.pylint.org), to check your code)
  - [ ] Is this new feature documented in the [Release
    Notes](https://github.com/pydy/pydy/tree/contributing#release-notes)?
  - [ ] Is it backwards compatible? All public methods/classes must follow
    deprecation cycles.
  - [ ] Have all reviewer comments been addressed?

Bug Fix
-------

::

  - [ ] Are there merge conflicts?
  - [ ] If there is an issue for this bug is it referenced in the commit
    message?
  - [ ] Have you created a unit test for the bug? (please reference the issue #
    in the unit test)
  - [ ] Did it pass the tests? Locally (run `nosetests`)? Travis CI?
  - [ ] Does it follow PEP8 guidelines? (use a linter, e.g.
    [pylint](http://www.pylint.org), to check your code)
  - [ ] Is this change documented in the [Release
    Notes](https://github.com/pydy/pydy/tree/contributing#release-notes)?
  - [ ] Have all reviewer comments been addressed?

Example: IPython Notebook
-------------------------

::

  - [ ] Are there merge conflicts?
  - [ ] Have you included a figure (preferrably SVG) for the example?
  - [ ] Have you noted the required versions of all dependencies? (The [version
    information extension](https://github.com/jrjohansson/version_information)
    is helpful for this.)
  - [ ] Have you included a `run.py` file that executes the example?
  - [ ] Have all reviewer comments been addressed?

Example: Script
---------------

::

  - [ ] Are there merge conflicts?
  - [ ] Have you included a figure (preferrably SVG) for the example?
  - [ ] Have you included the required versions for all of the dependencies?
  - [ ] Does the example execute with a single `python run.py` command?
  - [ ] Have all reviewer comments been addressed?

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
