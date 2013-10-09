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
import matplotlib
import numpy as np
matplotlib.rcParams['toolbar'] = 'None'
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import Arc
import matplotlib.patches as patches
import matplotlib.patheffects as PathEffects



# -----------------------------------------------------------------------------
def cubic_bezier(p0, p1, p2, p3, color = 'k', linewidth=1, alpha=1,
                 capstyle='butt', joinstyle='round'):
    """ """

    verts = np.array([p0,p1,p2,p3]).reshape(4,2)
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4 ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path,
                              edgecolor=color,
                              facecolor = 'None',
                              linewidth=linewidth,
                              alpha=alpha)
    patch.set_path_effects([PathEffects.Stroke(capstyle=capstyle)])
    plt.gca().add_patch(patch)

    
    plt.plot(verts[:2,0],verts[:2,1], color='k', lw=.5)
    plt.plot(verts[2:,0],verts[2:,1], color='k', lw=.5)
    plt.scatter(verts[:,0], verts[:,1], zorder=10,
                s=20, edgecolor='k', facecolor='w')
    plt.xticks([]), plt.yticks([])


# -----------------------------------------------------------------------------
def polyline(verts, color = 'b', linewidth=1, alpha=1,
             capstyle='butt', joinstyle='miter'):
    """ """

    verts = np.array(verts).reshape(len(verts),2)
    codes = [Path.MOVETO] + [Path.LINETO,]*(len(verts)-1)
    path = Path(verts, codes)
    patch = patches.PathPatch(path,
                              linewidth = linewidth,
                              edgecolor = color,
                              facecolor = 'None',
                              alpha     = alpha)
    patch.set_path_effects([PathEffects.Stroke(capstyle=capstyle,
                                               joinstyle=joinstyle )])
    plt.gca().add_patch(patch)
    plt.xticks([]), plt.yticks([])


# -----------------------------------------------------------------------------
def polyarc(arcs, color = 'b', linewidth=1, alpha=1,
            capstyle='butt', joinstyle='miter'):
    """ """

    for arc in arcs:
        center,radius, angle0, angle1, negative, a = arc
        angle0 = 180*angle0/math.pi
        angle1 = 180*angle1/math.pi
        if negative:
            angle0, angle1 = angle1, angle0
        arc = Arc(center, 2*radius, 2*radius, 0, angle0, angle1,
                  color=color, linewidth = linewidth, alpha=alpha)
        arc.set_path_effects([PathEffects.Stroke(capstyle=capstyle,
                                               joinstyle=joinstyle )])
        plt.gca().add_artist(arc)
    plt.xticks([])
    plt.yticks([])


# -----------------------------------------------------------------------------
def figure(width=800, height=800, on_key=None):
    """ """

    dpi = 72.0
    figsize= width/float(dpi), height/float(dpi)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
    if on_key:
        fig.canvas.mpl_connect('key_press_event', on_key)
    axes = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    axes.set_xlim(0, width)
    axes.set_ylim(0, height)
    plt.xticks([])
    plt.yticks([])
    plt.ion()
    plt.show()


# -----------------------------------------------------------------------------
if __name__ == '__main__':
    from distance import interpolate_path, polyline_to_cubic
    from cubic_bezier import CubicBezier
    from vec2 import vec2
    from bezier_type import cubic_bezier_type

    def on_key(*args):

        P = np.random.randint(100,700,(4,2))
        C = CubicBezier(*P)

        plt.cla()
        plt.ion()

        cubic_bezier(C.p0,C.p1,C.p2,C.p3)

        # P = C.flatten_forward_iterative(n=50)
        # P = np.array(P)
        # polyline(P, linewidth=100, alpha=.25)
        # plt.scatter(P[:,0], P[:,1], s=25,
        #             edgecolor='k', facecolor='w', zorder=10, lw=.1)

        P = C.flatten_iterative(flatness=.125, angle=15)
        P = np.array(P)
        polyline(P, linewidth=100, alpha=.1)
        plt.scatter(P[:,0], P[:,1], s=25, edgecolor='k', facecolor='w', zorder=10, lw=.5)

        # P = C.flatten_recursive(flatness=.125, angle=10) 
        # P = np.array(P)
        # polyline(P, linewidth=100, alpha=.1)
        # plt.scatter(P[:,0], P[:,1], s=25, edgecolor='k', facecolor='w', zorder=10, lw=.5)

        #A = C.flatten_behdad_arc(0.125)
        #polyarc(A, linewidth=100, alpha=.1)

        print "%d,%d %d,%d %d,%d, %d,%d : %s" % (
            C.p0[0],C.p0[1],
            C.p1[0],C.p1[1],
            C.p2[0],C.p2[1],
            C.p3[0],C.p3[1],
            cubic_bezier_type( vec2(C.p0[0],C.p0[1]),
                                         vec2(C.p1[0],C.p1[1]),
                                         vec2(C.p2[0],C.p2[1]),
                                         vec2(C.p3[0],C.p3[1]) ))


        plt.ioff()

    fig = figure(on_key=on_key)
    plt.ioff()
    plt.show()

