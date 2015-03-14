The PyDy source code is actually spread across both the SymPy and PyDy
repositories. The two projects are very intimately tied together in terms of
development. Much of our development practices come from SymPy's methods so
reading their documentation for contributing is the best place to start:

https://github.com/sympy/sympy/wiki/introduction-to-contributing

We do have some differences though. In particular, we use a checklist to make
sure each pull request is ready to be merged into the master branch (which
should always be stable). In particular, every pull request into this
repository must satisfy a checklist:

- [ ] Are there no merge conflicts?
- [ ] Did it pass the tests? Locally? Travis CI? (run ``nosetests``)
- If this introduces new functionality:
  - [ ] Have you added unit tests?
  - [ ] Do all public methods/classes have docstrings? (We use the [numpydoc
    format](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).)
  - [ ] Did you add an explanation to the online documentation?
- [ ] Does it follow PEP8 guidelines? (use a linter, e.g.
  [pylint](http://www.pylint.org), to check your code)
- [ ] Is this change documented in the [Release
  Notes](https://github.com/pydy/pydy/tree/contributing#release-notes)?
- [ ] Have all reviewer comments been addressed?

This checklist should be pasted in the top most comment block of the pull
request. It is then the duty of reviewers to check off the items as they are
completed. If any questions are not relevant to your particular pull request, a
reviewer will simply check it off as done.
