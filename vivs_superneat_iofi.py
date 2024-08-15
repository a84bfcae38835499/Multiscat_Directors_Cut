import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import datetime
import unicodedata
import re
import matplotlib.patheffects as pe
import matplotlib.patches as patches
import math

import matplotlib.cm as cm

n1n2OfInterest = [[1,0],
                  [0,0],
                  [-1,0],
                  [-2,0],
                  [-3,0],
                  [0,1],
                  [1,-1],
                  [1,-2],
                  [-1,2]]
n1n2Colours = [[1, 0.82, 0.149],
                  [1, 1, 1],
                  [0.067, 0.769, 0.451],
                  [0.149, 0.792, 0.7],
                  [0.296, 0.369, 1],
                  [1, 0.463, 0.722],
                  [0.518, 0.51, 0.941],
                  [0.596, 0, 1],
                  [1, 0.18, 0.408]]

Ninterest = sum(1 for _ in n1n2OfInterest)
Ninplane = 5
Noutplane = Ninterest - Ninplane
n1n2Inplane = n1n2OfInterest[:Ninplane]
n1n2Outplane = n1n2OfInterest[Ninplane:]
print("===")
print(n1n2Inplane)
print("---")
print(n1n2Outplane)
print("===")

Thetas = np.array([
    0/36,
    1/36,
    2/36,
    3/36,
    4/36,
    5/36,
    6/36,
    7/36
])
NIs = sum(1 for _ in n1n2OfInterest)

Is = np.arrat([
    
])

fig = plt.figure()
gs = fig.add_gridspec(2, hspace=0)
(axIn,axOut) = gs.subplots(sharex=True, sharey=True)
axI.plot()