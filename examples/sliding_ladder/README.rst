Sliding Ladder in a Frictionless World
======================================

We determine the equations of motion of a sliding ladder falling under the
action of gravity, supported by a frictionless wall and a frictionless floor.

The python file uses Kane's method to derive the equations of motion of the
ladder. These equations of motion are then integrated in matlab and the results
are used to animate the ladder.

Run the example with::

   $ python sliding_ladder.py
   $ matlab -nodisplay -nosplash -nodesktop -r "run('Root_SlidingLadder.m');exit;"
