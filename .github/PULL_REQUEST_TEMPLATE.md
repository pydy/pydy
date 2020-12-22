Please edit below based on the type of PR. It is then the duty of reviewers to check off the items as they are completed. If any questions are not relevant to your particular pull request, a reviewer will simply check it off as done.

New Feature
-----------

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

   - [ ] There are no merge conflicts with the `master` branch.
   - [ ] Examples rst file is placed in `docs/examples/`.
   - [ ] An entry with a 200x200 pixel image has been added to
     `docs/index.rst`.
   - [ ] RsT and image filenames are kebab-case (lower case hyphen-separated).
   - [ ] jupyter_sphinx RsT directives are used to ensure example is executed
     on build.
   - [ ] jupyter_sphinx `ipynb` and `py` download links are included.
   - [ ] The example has adequate text and figures (preferably SVG) explaining
     the problem.
   - [ ] All reviewer comments have been addressed.
