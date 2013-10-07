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
import sys
import os.path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import distance
from cubic_bezier import CubicBezier

NTESTS = 10000

# Measure errors on 10,000 curves
np.random.seed(1)
curves = np.random.randint(100,700,(NTESTS,4,2))
flatness = 0.125
angle = 15
n1, n2 = 25, 50


# -------------------------------------
def update_progress(progress):
    """
    Displays or updates a console progress bar

    Accepts a float between 0 and 1. Any int will be converted to a float.
    A value under 0 represents a 'halt'.
    A value at 1 or bigger represents 100%
    """

    barLength = 50
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "Error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rProgress: [%s] %.2f%% %s" % ("="*block + " "*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


# Forward method, n=25
# -------------------------------------
filename = 'forward-iterative-25.npy'
if not os.path.exists(filename):
    print "Computing", filename
    E1 = []
    for i in range(NTESTS):
        update_progress(i/float(NTESTS))
        p0,p1,p2,p3 = curves[i]
        C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        P = C.flatten_forward_iterative(n=n1)
        d = distance.polyline_to_cubic(P, p0, p1, p2, p3, n=100)
        E1.append(d)
    update_progress(1)
    np.save(filename, E1)
else:
    print "Loading", filename
    E1 = np.load(filename)


# Forward method, n=50
# -------------------------------------
filename = 'forward-iterative-50.npy'
if not os.path.exists(filename):
    print "Computing", filename
    E2 = []
    for i in range(NTESTS):
        update_progress(i/float(NTESTS))
        p0,p1,p2,p3 = curves[i]
        C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        P = C.flatten_forward_iterative(n=n2)
        d = distance.polyline_to_cubic(P, p0, p1, p2, p3, n=100)
        E2.append(d)
    update_progress(1)
    np.save(filename, E2)
else:
    print "Loading", filename
    E2 = np.load(filename)

# Smart iterative
# -------------------------------------
filename = 'smart-iterative.npy'
if not os.path.exists(filename):
    print "Computing", filename
    E3 = []
    for i in range(NTESTS):
        update_progress(i/float(NTESTS))
        p0,p1,p2,p3 = curves[i]
        C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        P = C.flatten_iterative(flatness=flatness, angle=angle)
        d = distance.polyline_to_cubic(P, p0, p1, p2, p3, n=100)
        E3.append(d)
    update_progress(1)
    np.save(filename, E3)
else:
    print "Loading", filename
    E3 = np.load(filename)

# Recursive
# -------------------------------------
filename = 'recursive.npy'
if not os.path.exists(filename):
    print "Computing", filename
    E4 = []
    for i in range(NTESTS):
        update_progress(i/float(NTESTS))
        p0,p1,p2,p3 = curves[i]
        C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        P = C.flatten_recursive(flatness=flatness, angle=angle)
        d = distance.polyline_to_cubic(P, p0, p1, p2, p3, n=100)
        E4.append(d)
    update_progress(1)
    np.save(filename, E4)
else:
    print "Loading", filename
    E4 = np.load(filename)


# Arcs
# -------------------------------------
filename = 'arc-iterative.npy'
if not os.path.exists(filename):
    print "Computing", filename
    E5 = []
    for i in range(NTESTS):
        update_progress(i/float(NTESTS))
        p0,p1,p2,p3 = curves[i]
        C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])
        P = C.flatten_behdad_arc(flatness=flatness)
        d = distance.polyarc_to_cubic(P, p0, p1, p2, p3, n=100)
        E5.append(d)
    update_progress(1)
    np.save(filename, E5)
else:
    print "Loading", filename
    E5 = np.load(filename)


def histogram_error(E, title):
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
    ax.spines['bottom'].set_position(('data',-NTESTS / 2000.))

    plt.hist(E[:,0],range=(0.0,0.16), bins=50, edgecolor='w', facecolor='.5')

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_yticks([])

    plt.xlim(0.0, 0.16)
    plt.ylim(0, NTESTS / 4)

    ax.yaxis.set_major_locator(MultipleLocator(NTESTS / 20))
    ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    [t.set_color('0.5') for t in ax.xaxis.get_ticklabels()]
    [t.set_color('0.5') for t in ax.yaxis.get_ticklabels()]
    [t.set_alpha(0.0) for t in ax.yaxis.get_ticklines()]

    plt.text(0.16,NTESTS * .2520, title,
             va='bottom',ha='right', color='0.0', fontsize=16)
    plt.text(0.16,NTESTS * .2475, 'Mean error over 10,000 curves',
             va='top',ha='right', color='.5', fontsize=12)

    M = E[:,0].mean()

    plt.axvline(x=M,ymin=0,ymax=NTESTS / 4, color='0.5',lw=.75,zorder=-1,ls='--')
    plt.text(0.0,NTESTS * .2510, '# Curves', va='bottom',ha='left', color='.5', fontsize=12)
    return fig


def histogram_length(E, title):
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
    ax.spines['bottom'].set_position(('data',-NTESTS / 2000.))

    plt.hist(E[:,2], bins=15, edgecolor='w', facecolor='.5')

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_yticks([])

    plt.xlim(0, 60)
    plt.ylim(0, NTESTS / 4)

    ax.yaxis.set_major_locator(MultipleLocator(NTESTS / 50))
    ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
    [t.set_color('0.5') for t in ax.xaxis.get_ticklabels()]
    [t.set_color('0.5') for t in ax.yaxis.get_ticklabels()]
    [t.set_alpha(0.0) for t in ax.yaxis.get_ticklines()]

    plt.text(60,NTESTS * .2420, title,
             va='bottom',ha='right', color='0.0', fontsize=16)
    plt.text(60,NTESTS * .2390, 'Mean size over 10,000 curves',
             va='top',ha='right', color='.5', fontsize=12)
    M = E[:,2].mean()
    plt.axvline(x=M,ymin=0,ymax=NTESTS / 4, color='0.5',lw=.75,zorder=-1,ls='--')
    plt.text(0.0,NTESTS * .2390, '# Curves', va='top',ha='left', color='.5', fontsize=12)
    return fig


# fig = histogram_length(E3, "Smart iterative")
# fig.savefig("smart-iterative-size.pdf", dpi=72)
# plt.show()

# fig = histogram_length(E4, "Recursive (agg)")
# fig.savefig("recursive-size.pdf", dpi=72)
# plt.show()

# fig = histogram_length(E5, "Arc iterative (glyphy)")
# fig.savefig("arc-iterative-size.pdf", dpi=72)
# plt.show()

# fig = histogram_error(E1, "Forward iterative (n=25)")
# fig.savefig("forwar-iterative-25-error.pdf", dpi=72)
# plt.show()

# fig = histogram_error(E2, "Forward iterative (n=50)")
# fig.savefig("forwar-iterative-50-error.pdf", dpi=72)
# plt.show()

# fig = histogram_error(E3, "Smart iterative")
# fig.savefig("smart-iterative-error.pdf", dpi=72)
# plt.show()

# fig = histogram_error(E4, "Recursive (agg)")
# fig.savefig("recursive-error.pdf", dpi=72)
# plt.show()

fig = histogram_error(E5, "Arc iterative (glyphy)")
fig.savefig("arc-iterative-error.pdf", dpi=72)
plt.show()
