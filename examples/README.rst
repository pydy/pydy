These are examples of using PyDy to derive, simulate, and study the motion of
classical dynamic systems. The equations of motion for the systems are
typically derived with SymPy Mechanics in symbolic form and then numerical
analyses is done with PyDy and various other tools in the SciPy Stack. Although
some examples also show cross language support for numerical analyses.

Each folder contains the files for one example. To contribute an example, make
a pull request with a new directory. The new folder should include, at the
minimum, a README explaining the problem, a figure (preferably SVG), and the
source code for the example either in script form or as an IPython notebook.
There should also be a file named `run.py` that executes the example.

Script Format
=============

There should be a REAMDE that explains the example, what each file is, and how
to run the example. Use as many files as needed to organize the example. Python
files should be PEP8 compliant. The `run.py` file should execute the entire
example.

IPython Notebook Format
=======================

Title
  Please provide some informative title for the example.
Problem Setup
  This should include a vector drawing of the system as well as a basic
  description of said drawing and the premise of the example.
Equations of Motion
  This should lead the reader through the derivation of the symbolic equations
  of motion (non-linear and/or linear).
Numerical Analyses
  If the example includes numerical analyses it should be explained here.
Visualization
  If the example includes a visualization of the motion include it here.

Style Guidelines:

  - The code cells should be PEP8 compliant.
  - Continuous sections of prose should be merged into the same cell (don't
    have two adjacent text cells).
  - The title should be heading 1 (single #) and in its own cell.
  - Use the "heading" cell type of proper level for organizing the sections.
