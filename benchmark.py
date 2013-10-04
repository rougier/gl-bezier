#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ----------------------------------------------------------------------------
# Copyright (C) 2013 Nicolas P. Rougier. All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY NICOLAS P. ROUGIER ''AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
# EVENT SHALL NICOLAS P. ROUGIER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation are
# those of the authors and should not be interpreted as representing official
# policies, either expressed or implied, of Nicolas P. Rougier.
# ----------------------------------------------------------------------------
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from cubic_bezier import CubicBezier
from curves import curve4_bezier


def point_line_distance(x1,y1, x2,y2, x3,y3): # x3,y3 is the point
    """ """

    px,py = x2-x1, y2-y1
    d = px*px + py*py
    u =  ((x3 - x1) * px + (y3 - y1) * py) / float(d)
    if u > 1:   u = 1
    elif u < 0: u = 0
    x,y = x1 + u * px, y1 + u * py
    dx, dy = x - x3, y - y3
    return math.sqrt(dx*dx + dy*dy)

def measure_error_polylines(samples=10000, n=100, flatness=0.125):
    """ """

    E = np.zeros((samples,3+8))
    for k in range(samples):
        points = np.random.randint(100,700,8)
        bezier = CubicBezier(*points)
        segments = bezier.flatten_iterative(flatness=flatness)
        D = np.zeros(n)
        for i,t in enumerate(np.linspace(0,1,n,endpoint=True)):
            x,y = bezier(t)
            dmin = 1e9
            for j in range(len(segments)-1):
                x1,y1 = segments[j]
                x2,y2 = segments[j+1]
                dmin = min(dmin,point_line_distance(x1,y1,x2,y2,x,y))
            D[i] = dmin
        E[k][0] = D.mean()
        E[k][1] = D.std()
        E[k][2] = len(segments)
        E[k][3:] = points

        if E[k][0] > 1.0:
            print E[k][0], points

    return E


def measure_error_polylines(samples=10000, n=100, flatness=0.125):
    """ """

    E = np.zeros((samples,3+8))
    for k in range(samples):
        points = np.random.randint(100,700,8)
        bezier = CubicBezier(*points)

        # segments = bezier.flatten_iterative(flatness=flatness)
        # segments = curve4_bezier(bezier.p0, bezier.p1, bezier.p2, bezier.p3)
        segments = bezier.flatten_forward_iterative(n=25)

        D = np.zeros(n)
        for i,t in enumerate(np.linspace(0,1,n,endpoint=True)):
            x,y = bezier(t)
            dmin = 1e9
            for j in range(len(segments)-1):
                x1,y1 = segments[j]
                x2,y2 = segments[j+1]
                dmin = min(dmin,point_line_distance(x1,y1,x2,y2,x,y))
            D[i] = dmin
        E[k][0] = D.mean()
        E[k][1] = D.std()
        E[k][2] = len(segments)
        E[k][3:] = points

        if E[k][0] > 1.0:
            print E[k][0], points

    return E


# E = measure_error_polylines(samples=10000, n=100, flatness=.125)
# np.save('flatten_iterative_error_s10000_n100_f0125.npy', E)
# E = np.load('flatten_iterative_error_s10000_n100_f0125.npy')

# E = measure_error_polylines(samples=10000, n=100, flatness=.125)
# np.save('flatten_recursive_error_s10000_n100_f0125.npy', E)
# E = np.load('flatten_recursive_error_s10000_n100_f0125.npy')

# E = measure_error_polylines(samples=10000, n=100, flatness=.125)
# np.save('flatten_forward_error_s10000_n100_n25.npy', E)
# E = np.load('flatten_forward_error_s10000_n100_n25.npy')
E = np.load('flatten_forward_error_s10000_n100_n50.npy')


fg = 0.0,0.0,0.0
bg = 1.0,1.0,1.0
matplotlib.rcParams['xtick.major.width'] = .5
matplotlib.rcParams['ytick.major.width'] = .5
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['font.size'] = 12.0
matplotlib.rc('axes', facecolor = bg)
matplotlib.rc('axes', edgecolor = fg)
matplotlib.rc('xtick', color = fg)
matplotlib.rc('ytick', color = fg)
matplotlib.rc('figure', facecolor = bg)
matplotlib.rc('savefig', facecolor = bg)


fig = plt.figure(figsize=(12,6), dpi=72)
ax = plt.subplot(111, axisbelow=True)
ax.spines['bottom'].set_position(('data',-5))

plt.hist(E[:,0],range=(0.0,0.16), bins=50, edgecolor='w', facecolor='.5')


ax.spines['right'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_yticks([])

plt.xlim(0.0, 0.16)
plt.ylim(0, 2600)

ax.yaxis.set_major_locator(MultipleLocator(500))
ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
[t.set_color('0.5') for t in ax.xaxis.get_ticklabels()]
[t.set_color('0.5') for t in ax.yaxis.get_ticklabels()]
[t.set_alpha(0.0) for t in ax.yaxis.get_ticklines()]
plt.text(0.16,2520, 'Brute iterative (n=50)',
                    va='bottom',ha='right', color='0.0', fontsize=16)
plt.text(0.16,2475, 'Mean error over 10,000 curves',
                    va='top',ha='right', color='.5', fontsize=12)

#M = E[:,0].mean()
#plt.axvline(x=M,ymin=0,ymax=2600, color='0.5',lw=.75,zorder=-1,ls='--')

plt.text(0.0,2510, '# Curves', va='bottom',ha='left', color='.5', fontsize=12)
fig.savefig("forward-iterative-50.pdf", dpi=72)
plt.show()
