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


# ----------------------------------------------------------- point_to_line ---
def point_to_line(p, start, end):
    """
    Point to line distance

    Parameters
    ----------

    p : tuple of 2 floats
        Point to measure distance from

    start : tuple of 2 floats
        Line start

    end : tuple of 2 floats
        Line end
    """

    px, py = end[0]-start[0], end[1]-start[1]
    d = px*px + py*py
    u = ((p[0]-start[0])*px + (p[1]-start[1])*py)/float(d)
    if u > 1:
        u = 1
    elif u < 0:
        u = 0
    dx, dy = start[0]+u*px-p[0], start[1]+u*py-p[1]

    return math.sqrt(dx*dx + dy*dy)


# ------------------------------------------------------------ point_to_arc ---
def point_to_arc(p, center, radius, angles):
    """
    Point to arc distance

    Parameters
    ----------

    p : tuple of 2 floats
        Point to measure distance from

    center : tuple of 2 floats
        Arc center

    radius : float
        Arc radius

    angles : tuple of 2 floats
        Arc angles (radians)
    """

    dx, dy = p[0]-center[0], p[1]-center[1]
    angle0, angle1 = angles

    # angle0 = math.fmod(angle0+2*math.pi,2*math.pi)
    # angle1 = math.fmod(angle1+2*math.pi,2*math.pi)
    # if not angle0 <= angle1:
    #     angle0,angle1 = angle1,angle0
    angle = math.atan2(dy,dx)
    angle = math.fmod(angle+2*math.pi,2*math.pi)

    # Distance to implicit circle
    #if angle0 <= angle <= angle1:
    d = abs(math.sqrt((p[0]-center[0])**2 + (p[1]-center[1])**2) - radius)
    return d
    
    # Distance to arc limit (points)
    #p0 = center[0]+radius*math.cos(angle0), center[1]+radius*math.sin(angle0)
    #d0 = (p0[0]-p[0])**2 + (p0[1]-p[1])**2
    #p1 = center[0]+radius*math.cos(angle1), center[1]+radius*math.sin(angle1)
    #d1 = (p1[0]-p[0])**2 + (p1[1]-p[1])**2

    return math.sqrt(min(d0,d1))


# -------------------------------------------------------- interpolate_path ---
def interpolate_path(points, n=100):
    """
    Regular interpolation of a path

    Paramaters
    ----------

    point : list of points
        Point describing the path

    n : int
        Numnber of equidistant points to interpolate
    """

    points = np.array(points)
    X,Y = points[:,0], points[:,1]

    # Compute tangents
    Tx = (X[1:] - X[:-1]).reshape(len(X)-1,1)
    Ty = (Y[1:] - Y[:-1]).reshape(len(Y)-1,1)

    # Compute segment lengths
    L = np.sqrt((Tx*Tx+Ty*Ty).sum(axis=1))

    # Compute curvilinear coordinate
    C = np.zeros(len(X))
    C[1:] = L.cumsum()
    length = C[-1]

    # Interpolate
    T = np.linspace(0,length,n,endpoint=True)

    P_ = np.zeros((n,2))
    P_[:,0] = np.interp(T,C,X)
    P_[:,1] = np.interp(T,C,Y)

    return P_


# ------------------------------------------------------- interpolate_cubic ---
def interpolate_cubic(p0, p1, p2, p3, n=100):
    """
    Non-regular interpolatation of a cubic Bézier curve

    p0,p1,p2,p3: tuple of 2 floats
        cubic Bézier points

    n : int
        number of point to generate
    """

    t = np.linspace(0,1,n,endpoint=True).reshape(n,1)
    m_t = 1 - t
    b = m_t * m_t
    c = t * t
    d = c * t
    a = b * m_t
    b *= 3 * t
    c *= 3 * m_t
    return a*p0 + b*p1 + c*p2 + d*p3


# ------------------------------------------------------- polyline_to_cubic ---
def polyline_to_cubic(points, p0, p1, p2, p3, n=100):
    """
    Polyline to cubic bézier mean distance

    Parameters
    ----------

    points: list of points
        points describing the polyline

    p0,p1,p2,p3: tuple of 2 floats
        cubic Bézier points

    n : int
        number of point to generate to measure distance
    """

    # Sample n points on cubic
    P = interpolate_cubic(p0,p1,p2,p3,n)

    # Make them regularly spaced
    # P = interpolate_path(P,n)

    D = np.zeros(n)
    for i,p in enumerate(P):
        dmin = 1e9
        for j in range(len(points)-1):
            d = point_to_line(p, points[j], points[j+1])
            dmin = min(d,dmin)
        D[i] = dmin
    return D.mean(), D.std(), len(points)



# -------------------------------------------------------- polyarc_to_cubic ---
def polyarc_to_cubic(arcs, p0, p1, p2, p3, n=100):
    """
    Polyarc to cubic bézier mean distance

    arc: list of arcs
        arcs describing the polyarcs

    p0,p1,p2,p3: tuple of 2 floats
        cubic Bézier points
    """

    import behdad

    # Sample n points on cubic
    P = interpolate_cubic(p0,p1,p2,p3,n)

    # Make them regularly spaced
    # P = interpolate_path(P,n)

    D = np.zeros(n)
    for i,p in enumerate(P):
        dmin = 1e9
        for arc in arcs:
            center,radius,angle0,angle1,negative,a = arc
            if negative:
               angle0, angle1 = angle1, angle0
            d = a.distance_to_point (behdad.Point(*p))
            dmin = min(d,dmin)
        D[i] = dmin
    return D.mean(), D.std(), len(arcs)



# ------------------------------------------------------------------------------
if __name__ == '__main__':
    import matplotlib
    import matplotlib.pyplot as plt
    from cubic_bezier import CubicBezier

    p0,p1,p2,p3 = np.random.randint(100,700,(4,2))
    C = CubicBezier(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1],p3[0],p3[1])

    P = C.flatten_forward_iterative(n=25)
    print polyline_to_cubic(P, p0, p1, p2, p3, n=100)

    P = C.flatten_forward_iterative(n=50)
    print polyline_to_cubic(P, p0, p1, p2, p3, n=100)

    P = C.flatten_iterative(0.125)
    print polyline_to_cubic(P, p0, p1, p2, p3, n=100)

    P = C.flatten_recursive(0.125)
    print polyline_to_cubic(P, p0, p1, p2, p3, n=100)

    A = C.flatten_behdad_arc(0.125)
    print polyarc_to_cubic(A, p0, p1, p2, p3, n=100)

    # # Check distance to arc
    # T = np.linspace(0,2*math.pi,100)
    # for t in T:
    #     x,y = math.cos(t), math.sin(t)
    #     angle = math.atan2(y,x)
    #     angle = math.fmod(angle+2*math.pi,2*math.pi)
    #     print 180*t/math.pi, point_to_arc((x,y), (0,0), 1, (0, math.pi/2))
