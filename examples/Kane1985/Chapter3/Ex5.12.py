#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Exercise 5.12 from Kane 1985."""

from __future__ import division
import numpy as np

# n_a = 3/5*n_1 - 4/5*n_3. Substituting n_i for e_i results in
# n_a = 4/5*e_1 + 3/5*n_2.
a = np.matrix([4/5, 3/5, 0])
Iij = np.matrix([[169, 144, -96],
                 [144, 260, 72],
                 [-96, 72, 325]])

print("Moment of inertia of B with respect to a line that is parallel to")
print("line PQ and passes through point O.")
print("{0} kg m**2".format((a * Iij * a.T).item(0)))
