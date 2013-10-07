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
from curves import curve4_bezier


class CubicBezier(object):
    """
    P(t) = (1-t)^3*P0 + 3*(1-t)^2*t*P1 + 3(1-t)*t*t*P2 + t^3*P3
    """

    def __init__(self, x0=0, y0=0, x1=0, y1=0, x2=0, y2=0, x3=0, y3=0):
        """
        """

        self.x0 = float(x0)
        self.y0 = float(y0)
        self.x1 = float(x1)
        self.y1 = float(y1)
        self.x2 = float(x2)
        self.y2 = float(y2)
        self.x3 = float(x3)
        self.y3 = float(y3)
 
    @property
    def p0(self):
        return np.array([self.x0,self.y0])

    @property
    def p1(self):
        return np.array([self.x1,self.y1])

    @property
    def p2(self):
        return np.array([self.x2,self.y2])

    @property
    def p3(self):
        return np.array([self.x3,self.y3])


    def __call__(self, t):
        """
        Evaluate bezier curve at t
        """
        m_t = 1.0 - t
        b = m_t * m_t
        c = t * t
        d = c * t
        a = b * m_t
        b *= 3. * t
        c *= 3. * m_t
        return (a*self.x0 + b*self.x1 + c*self.x2 + d*self.x3,
                a*self.y0 + b*self.y1 + c*self.y2 + d*self.y3)
        
    def split(self, t):
        """
        Split curve at t into left and right cubic bezier curves
        """

        left, right = CubicBezier(), CubicBezier()
        left.x0 = self.x0
        left.y0 = self.y0
        left.x1 = self.x0 + t*(self.x1-self.x0)
        left.y1 = self.y0 + t*(self.y1-self.y0)
        left.x2 = self.x1 + t*(self.x2-self.x1)
        left.y2 = self.y1 + t*(self.y2-self.y1)
        right.x2 = self.x2 + t*(self.x3-self.x2)
        right.y2 = self.y2 + t*(self.y3-self.y2)
        right.x1 = left.x2 + t*(right.x2-left.x2)
        right.y1 = left.y2 + t*(right.y2-left.y2)
        left.x2 = left.x1 + t*(left.x2-left.x1)
        left.y2 = left.y1 + t*(left.y2-left.y1)
        left.x3 = right.x0 = left.x2 + t*(right.x1-left.x2)
        left.y3 = right.y0 = left.y2 + t*(right.y1-left.y2)
        right.x3 = self.x3
        right.y3 = self.y3        
        return left, right


    def extract(self, t1, t2):
        """
        Extract a cubic bezier curve between t1 and t2
        """
        _,right = self.split(t1)
        left,_= right.split((t2-t1)/(1-t1))
        return left


    def inflection_points(self):
        """
        Find inflection points
        """

        A = -  self.p0 + 3*self.p1 - 3*self.p2 + self.p3
        B =  3*self.p0 - 6*self.p1 + 3*self.p2
        C = -3*self.p0 + 3*self.p1
        a = 3*np.cross(A,B)
        b = 3*np.cross(A,C)
        c =   np.cross(B,C)
        r = b*b - 4*a*c
        if r >= 0 and a:
            r = np.sqrt(r)
            ip1 = (-b + r) / (2*a)
            ip2 = (-b - r) / (2*a)
            return min(ip1,ip2), max(ip1,ip2)
        return None, None


    def inflection_domain(self, t, flatness=0.25):
        """
        Determine the domain around an inflection point where the curve is flat.
        """

        _, right = self.split(t)
        ax = -right.x0 + 3*right.x1 - 3*right.x2 + right.x3
        ay = -right.y0 + 3*right.y1 - 3*right.y2 + right.y3
        ex = 3 * (right.x1 - right.x2)
        ey = 3 * (right.y1 - right.y2)
        ex2ey2 = ex * ex + ey * ey
        if not ex2ey2:
            return t,t

        s4 = abs(6. * (ey * ax - ex * ay) / math.sqrt(ex2ey2))
        tf = math.pow(9. * flatness / s4, 1./3.)
        return t-tf*(1-t), t+tf*(1-t)

    def diamond_angle(self,y,x):
        """ See http://jsperf.com/diamond-angle-vs-atan2/2 """
        if y >= 0:
            if x >=0:
                return y/(x+y)
            else:
                return 1-x/(-x+y)
        else:
            if x < 0:
                return 2-y/(-x-y)
            else:
                return 3+x/(x-y)

    def radian_to_diamond_angle(self,rad):
        return self.diamond_angle(math.cos(rad), math.sin(rad))


    def angle(self):
        a23 = math.atan2(self.y2 - self.y1, self.x2 - self.x1)
        da1 = abs(a23 - math.atan2(self.y1 - self.y0, self.x1 - self.x0))
        da2 = abs(math.atan2(self.y3 - self.y2, self.x3 - self.x2) - a23)
        if da1 >= math.pi:
            da1 = 2*math.pi - da1
        if da2 >= math.pi:
            da2 = 2*math.pi - da2
        return da1 + da2

    def flatten(self, flatness=0.25, angle=15):
        angle *= math.pi/180.0

        P = []
        while 1:
            dx = self.x1 - self.x0
            dy = self.y1 - self.y0
            norm = math.sqrt(dx * dx + dy * dy)
            if not norm:
                break
            s3 = (self.x2-self.x0)*(self.y1-self.y0)-(self.y2-self.y0)*(self.x1-self.x0)
            s3 = abs(s3)/norm
            t = 2*math.sqrt(flatness /(3*s3))
            if t > 1:
                break

            # Check angle is below tolerance
            for i in range(10):
                left, right = self.split(t)
                if left.angle() > angle:
                   t /= 2.0
                else:
                    break

            self.x0, self.y0 = right.x0, right.y0
            self.x1, self.y1 = right.x1, right.y1
            self.x2, self.y2 = right.x2, right.y2
            self.x3, self.y3 = right.x3, right.y3

            P.append((self.x0, self.y0))
        return P


    def flatten_brute_iterative(self, n=50):
        """
        """
        P = [(self.x0,self.y0)]
        for t in np.linspace(0,1,n-1,endpoint=False)[1:]:
            P.append(self(t))
        P.append((self.x3,self.y3))
        return P

    def flatten_forward_iterative(self, n=50):
        """
        """
        h = 1.0 / n;
        fph = 3 * (self.p1 - self.p0) * h
        fpphh = (6 * self.p0 - 12 * self.p1 + 6 * self.p2) * h * h
        fppphhh = (-6 * self.p0 + 18 * self.p1 - 18 * self.p2 + 6 * self.p3) * h * h * h
        P = [(self.x0,self.y0)]
        p = np.array([self.x0,self.y0])
        for i in range(1,n-1):
            p += fph + fpphh/2. + fppphhh/6.
            P.append((p[0],p[1]))
            fph = fph + fpphh + fppphhh/2.
            fpphh = fpphh + fppphhh
        P.append((self.x3,self.y3))
        return P

    def flatten_forward_iterative_variable(self):
        """
        """

        d1 = np.sqrt(((self.p1-self.p0)**2).sum())
        d2 = np.sqrt(((self.p2-self.p1)**2).sum())
        d3 = np.sqrt(((self.p3-self.p2)**2).sum())
        n = int((d1+d2+d3)*0.05)


        h = 1.0 / n;
        fph = 3 * (self.p1 - self.p0) * h
        fpphh = (6 * self.p0 - 12 * self.p1 + 6 * self.p2) * h * h
        fppphhh = (-6 * self.p0 + 18 * self.p1 - 18 * self.p2 + 6 * self.p3) * h * h * h
        P = [(self.x0,self.y0)]
        p = np.array([self.x0,self.y0])
        for i in range(1,n-1):
            p += fph + fpphh/2. + fppphhh/6.
            P.append((p[0],p[1]))
            fph = fph + fpphh + fppphhh/2.
            fpphh = fpphh + fppphhh
        P.append((self.x3,self.y3))
        return P


    def flatten_forward_iterative_variable_2(self):
        """
        """

        d1 = np.sqrt(((self.p1-self.p0)**2).sum())
        d2 = np.sqrt(((self.p2-self.p1)**2).sum())
        d3 = np.sqrt(((self.p3-self.p2)**2).sum())
        d4 = np.sqrt(((self.p0-self.p3)**2).sum())
        d = d1+d2+d3+d4
        f = max((d-d4)/float(d),.5)
        n = 20 + int(80 * (f-0.5))

        h = 1.0 / n;
        fph = 3 * (self.p1 - self.p0) * h
        fpphh = (6 * self.p0 - 12 * self.p1 + 6 * self.p2) * h * h
        fppphhh = (-6 * self.p0 + 18 * self.p1 - 18 * self.p2 + 6 * self.p3) * h * h * h
        P = [(self.x0,self.y0)]
        p = np.array([self.x0,self.y0])
        for i in range(1,n-1):
            p += fph + fpphh/2. + fppphhh/6.
            P.append((p[0],p[1]))
            fph = fph + fpphh + fppphhh/2.
            fpphh = fpphh + fppphhh
        P.append((self.x3,self.y3))
        return P

    def flatten_recursive(self,flatness=0.125, angle=15):
        return curve4_bezier(self.p0, self.p1, self.p2, self.p3, flatness, angle)


    def flatten_iterative(self, flatness=0.125, angle=15):
        """
        Adapted from: Precise Flattening of Cubic Bézier Segments
                      Thomas F. Hain, Athar L. Ahmad, David D. Langan

        The idea is to split the curve at inflection points such that
        each part has monotonic curvature and can be approximated
        using parabolic flattening. Inflection point vicinity is
        approximated using a line since it is flat by definition
        (curvature=0).

        There are at most 2 inflection points for a cubic curve so
        a curve can be splitted into 5 parts at most (3 monotonic
        curves and 2 segments).

        """
    
        t1_minus, t1_plus = -1,-1
        t2_minus, t2_plus = +2,+2
        T = self.inflection_points()
        if T[0]:
            t1_minus, t1_plus = self.inflection_domain(T[0],flatness)
        if T[1]:
            t2_minus, t2_plus = self.inflection_domain(T[1], flatness)

        # Split the two domaisn if they overlap
        if t1_minus < t2_minus < t1_plus:
            t1_plus, t2_minus = t2_minus, t1_plus

        t1_out = t1_plus < 0 or t1_minus > 1
        t2_out = t2_plus < 0 or t2_minus > 1
        t1_t2_cross = t2_minus <  t1_plus

        # Make sure the possible out domain is always t1
        #  (this will save some specific tests below)
        if t1_out:
            t1_minus, t1_plus = t2_minus, t2_plus
            t1_out = t2_out
            t2_minus, t2_plus = +2, +2
            t2_out = True

        t1_in          = 0 < t1_minus < t1_plus < 1
        t1_cross_start = t1_minus < 0 < t1_plus < 1
        t1_cross_end   = 0 < t1_minus < 1 < t1_plus
        t1_cross       = t1_cross_start or t1_cross_end

        t2_out         = t2_plus < 0 or t2_minus > 1
        t2_in          = 0 < t2_minus < t2_plus < 1
        t2_cross_start = t2_minus < 0 < t2_plus < 1
        t2_cross_end   = 0 < t2_minus < 1 < t2_plus
        t2_cross       = t2_cross_start or t2_cross_end

        tmp = CubicBezier(self.x0, self.y0, self.x1, self.y1,
                          self.x2, self.y2, self.x3, self.y3)

        points = [(self.x0, self.y0),]

        # No inflection points
        if t1_out and t2_out:
            points += tmp.flatten(flatness,angle)

        # One inflection point
        elif (t1_in or t1_cross) and t2_out:
            if t1_cross_start:
                points.append( self(t1_plus) )
                _,right = self.split(t1_plus)
                points += right.flatten(flatness,angle)
            elif t1_cross_end:
                left,_ = self.split(t1_minus)
                points += left.flatten(flatness,angle)
                points.append( self(t1_minus) )
            else:
                left,_ = self.split(t1_minus)
                _,right = self.split(t1_plus)
                points += left.flatten(flatness, angle)
                points.append( self(t1_minus) )
                points.append( self(t1_plus) )
                points += right.flatten(flatness, angle)

        # Two inflection points
        elif (t1_in or t1_cross_start) and (t2_in or t2_cross_end):
            if not t1_cross_start:
                left,_ = self.split(t1_minus)
                points += left.flatten(flatness, angle)
                points.append( self(t1_minus) )
            if t1_t2_cross:
                points.append( self(t2_minus) )
                points.append( self(t1_plus) )
            else:
                points.append( self(t1_plus) )
                middle = self.extract(t1_plus,t2_minus)
                points += middle.flatten(flatness, angle)
                points.append( self(t2_minus) )

            if not t2_cross_end:
                points.append( self(t2_plus) )
                _,right = self.split(t2_plus)
                points += right.flatten(flatness, angle)

        points.append( (self.x3, self.y3) )
        return points

    def flatten_behdad_segment(self, flatness=0.125):
        """
        """
        P = [(self.x0,self.y0)]
        import behdad

	b = behdad.Bezier (behdad.Point(self.x0,self.y0),
			   behdad.Point(self.x1,self.y1),
			   behdad.Point(self.x2,self.y2),
			   behdad.Point(self.x3,self.y3))

	devfunc = behdad.MaxDeviationApproximatorExact ()
	errfunc = behdad.ArcBezierErrorApproximatorBehdad (devfunc)
	apprfunc = behdad.ArcBezierApproximatorMidpointTwoPart (errfunc)
	splinefunc = behdad.ArcsBezierApproximatorSpringSystem (apprfunc)

        arcs, error = splinefunc (b, flatness)
	for arc in arcs:
            P.append((arc.p1.x,arc.p1.y))
        return P

    def flatten_behdad_arc(self, flatness=0.125):
        """
        """
        A = []
        import behdad
	b = behdad.Bezier (behdad.Point(self.x0,self.y0),
			   behdad.Point(self.x1,self.y1),
			   behdad.Point(self.x2,self.y2),
			   behdad.Point(self.x3,self.y3))

	devfunc = behdad.MaxDeviationApproximatorExact ()
	errfunc = behdad.ArcBezierErrorApproximatorBehdad (devfunc)
	apprfunc = behdad.ArcBezierApproximatorMidpointTwoPart (errfunc)
	splinefunc = behdad.ArcsBezierApproximatorSpringSystem (apprfunc)

	arcs, error = splinefunc (b, flatness)
	for arc in arcs:
            A.append(arc.to_conventional())
        return A


    def plot(self, *args, **kwargs):

        kwargs['linewidth'] = kwargs.get('linewidth', 1)
        kwargs['facecolor'] = kwargs.get('facecolor', 'None')
        kwargs['edgecolor'] = kwargs.get('edgecolor', 'k')
        kwargs['alpha']     = kwargs.get('alpha', 1)

        # Actual curve
        verts = np.array([(self.x0,self.y0), (self.x1,self.y1),
                          (self.x2,self.y2), (self.x3,self.y3)])
        codes = [Path.MOVETO, Path.CURVE4,
                 Path.CURVE4, Path.CURVE4 ]
        path = Path(verts, codes)
        patch = patches.PathPatch(path, *args, **kwargs )
        patch.set_path_effects([PathEffects.Stroke(joinstyle='round',
                                                   capstyle='butt')])
        axes.add_patch(patch)

        plt.plot(verts[:2,0],verts[:2,1], color='k', lw=.5)
        plt.plot(verts[2:,0],verts[2:,1], color='k', lw=.5)
        plt.scatter(verts[:,0], verts[:,1], zorder=10,
                    s=20, edgecolor='k', facecolor='w')

        I = self.inflection_points()
        for t in self.inflection_points():
            if t:
                t_neg,t_pos = self.inflection_domain(t)
                if 0 <= t <= 1:
                    P = self(t)
                    plt.scatter([P[0],], [P[1],], s=50,
                                edgecolor='r', facecolor='None')
                    plt.scatter([P[0],], [P[1],], s=50, zorder=10,
                                edgecolor='None', facecolor='white',alpha=.5)
                if 0 <= t_neg <= 1:
                    P = self(t_neg)
                    plt.scatter([P[0],], [P[1],], s=10, color='k')
                if 0 <= t_pos <= 1:
                    P = self(t_pos)
                    plt.scatter([P[0],], [P[1],], s=10, color='k')
            
            

                    
def line_plot(verts, *args, **kwargs):

    kwargs['linewidth'] = kwargs.get('linewidth', 1)
    kwargs['facecolor'] = kwargs.get('facecolor', 'None')
    kwargs['edgecolor'] = kwargs.get('edgecolor', 'k')
    kwargs['alpha']     = kwargs.get('alpha', 1)

    # Actual curve
    verts = np.array(verts).reshape(len(verts),2)
    codes = [Path.MOVETO] + [Path.LINETO,]*(len(verts)-1)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, *args, **kwargs )
    patch.set_path_effects([PathEffects.Stroke(capstyle='butt')])
    axes.add_patch(patch)


def path_interpolate(P, n, start=None, end=None):

    X,Y = P[:,0], P[:,1]

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
    start = start or 0
    end  = (end or 1)*length
    T = np.linspace(start,end,n,endpoint=True)

    P_ = np.zeros((n,2))
    P_[:,0] = np.interp(T,C,X)
    P_[:,1] = np.interp(T,C,Y)
    return P_



def point_line_distance(x1,y1, x2,y2, x3,y3): # x3,y3 is the point
    px = x2-x1
    py = y2-y1

    something = px*px + py*py

    u =  ((x3 - x1) * px + (y3 - y1) * py) / float(something)

    if u > 1:
        u = 1
    elif u < 0:
        u = 0

    x = x1 + u * px
    y = y1 + u * py

    dx = x - x3
    dy = y - y3

    # Note: If the actual distance does not matter,
    # if you only want to compare what this function
    # returns to other results of this function, you
    # can just return the squared distance instead
    # (i.e. remove the sqrt) to gain a little performance

    dist = math.sqrt(dx*dx + dy*dy)

    return dist

def measure(bezier, segments):
    n = 128
    D = np.zeros(n)
    for i,t in enumerate(np.linspace(0,1,n,endpoint=True)):
        x,y = bezier(t)
        dmin = 1e9
        for k in range(len(segments)-1):
            x1,y1 = segments[k]
            x2,y2 = segments[k+1]
            dmin = min(dmin,point_line_distance(x1,y1,x2,y2,x,y))
        D[i] = dmin

    return D.mean(), D.std()


def measure_arc(bezier, arcs):
    n = 128
    D = np.zeros(n)
    for i,t in enumerate(np.linspace(0,1,n,endpoint=True)):
        x,y = bezier(t)
        dmin = 1e9
        for arc in arcs:
            # Do we really need to check angles ?
            center,radius,angle0,angle1,negative = arc
            dx = (center.x - x)
            dy = (center.y - y)
            d = abs(math.sqrt(dx*dx+dy*dy) - radius)
            dmin = min(dmin, d)
        D[i] = dmin

    return D.mean(), D.std()

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    import matplotlib
    import numpy as np
    matplotlib.rcParams['toolbar'] = 'None'
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import matplotlib.patheffects as PathEffects
    from matplotlib.patches import Arc

    from curves import curve4_bezier

    size = 800,800

    dpi = 72.0
    figsize= size[0]/float(dpi),size[1]/float(dpi)
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
    axes = fig.add_axes([0.0, 0.0, 1.0, 1.0], frameon=False)
    axes.set_xlim(0,size[0])
    axes.set_ylim(0,size[1])

    def on_key(event=None):
        if event and event.key == 'escape':
            sys.exit()
        plt.cla()
        plt.ion()
        axes.set_xlim(0,size[0])
        axes.set_ylim(0,size[1])

        # P = [275, 591, 143, 383, 348, 578, 560, 285]
        # P = [291,612, 533,579, 482,476, 332,175]
        # P = [332,175,482,476,533,579,291,612]
        # P = [145,207,293,238,565,453,564,103]
        # P = [551,198,637,684,417,500,566,573]
        # P = [351,529,291,642,460,333,324,353]
        # P = [573,157,699,546,101,594,157,589]
        # P = [200,400, 600,600, 200,600, 600, 400]
        P = np.random.randint(100,700,8)

        C = CubicBezier(*P)

        # P = C.flatten_brute_iterative(25)
        P = C.flatten_forward_iterative(25)
        mean,std = measure(C,P)
        # print "Error forward   = %.2f +/- %.3f (%d points)" % (mean,std,len(P))

        P = C.flatten_iterative(flatness=.125)
        mean,std = measure(C,P)
        print "Error iterative = %.2f +/- %.3f (%d points)" % (mean,std,len(P))

        C.plot()
        P = np.array(P)
        line_plot(P, lw=150, edgecolor='k', alpha=.1)
        #plt.plot(P[:,0], P[:,1], lw=150, color='k', alpha=.1)
        plt.scatter(P[:,0], P[:,1], s=5, color='r', zorder=40)

        # P = C.flatten_behdad_segment(flatness=.125)
        # mean,std = measure(C,P)
        # print "Error behdad    = %.2f +/- %.3f (%d points)" % (mean,std,len(P))

        # A = C.flatten_behdad_arc(flatness=.125)
        # for arc in A:
        #     center,radius,angle0,angle1,negative = arc
        #     angle0 = 180*angle0/math.pi
        #     angle1 = 180*angle1/math.pi
        #     if negative:
        #         angle0, angle1 =  angle1, angle0

        #     arc = Arc((center.x,center.y), 2*radius, 2*radius, 0,
        #               angle0, angle1, linewidth = 150, alpha=.25)
        #     fig.gca().add_artist(arc)
        # mean,std = measure_arc(C,A)
        # print "Error behdad arc= %.2f +/- %.3f (%d points)" % (mean,std,len(A))

        P = curve4_bezier(C.p0, C.p1, C.p2, C.p3)
        mean,std = measure(C,P)
        print "Error recursive = %.2f +/- %.3f (%d points)" % (mean,std,len(P))
        #print "Error recursive = %.2f (%d points)" % (measure(C,P), len(P))


        print "C = [%d,%d,%d,%d,%d,%d,%d,%d]" % (C.x0,C.y0,C.x1,C.y1,C.x2,C.y2,C.x3,C.y3)
        print

        plt.xticks([]),plt.yticks([])
        plt.ioff()

    on_key()
    fig.canvas.mpl_connect('key_press_event', on_key)
    plt.xticks([]),plt.yticks([])
    plt.show()
