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
from vec2 import *
from curves import curve4_bezier


class CubicBezier(object):
    """
    P(t) = (1-t)^3*P0 + 3*(1-t)^2*t*P1 + 3(1-t)*t*t*P2 + t^3*P3
    """

    def __init__(self, (x0,y0)=(0,0), (x1,y1)=(0,0), (x2,y2)=(0,0), (x3,y3)=(0,0)):
        """ """
        self.x0, self.y0 = float(x0), float(y0)
        self.x1, self.y1 = float(x1), float(y1)
        self.x2, self.y2 = float(x2), float(y2)
        self.x3, self.y3 = float(x3), float(y3)
 
    @property
    def p0(self):
        return vec2(self.x0,self.y0)

    @property
    def p1(self):
        return vec2(self.x1,self.y1)

    @property
    def p2(self):
        return vec2(self.x2,self.y2)

    @property
    def p3(self):
        return vec2(self.x3,self.y3)


    def __call__(self, t):
        """ Evaluate curve at t """
        m_t = 1.0 - t
        b = m_t * m_t
        c = t * t
        d = c * t
        a = b * m_t
        b *= 3. * t
        c *= 3. * m_t
        return vec2(a*self.x0 + b*self.x1 + c*self.x2 + d*self.x3,
                    a*self.y0 + b*self.y1 + c*self.y2 + d*self.y3)
        
    def split(self, t):
        """ Split curve at t into left and right cubic bezier curves """

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


    def segment(self, t1, t2):
        """ Segment a cubic bezier curve between t1 and t2 """

        _,right = self.split(t1)
        left,_= right.split((t2-t1)/(1-t1))
        return left


    def inflection_points(self):
        """ Find inflection points """

        A = -  self.p0 + 3*self.p1 - 3*self.p2 + self.p3
        B =  3*self.p0 - 6*self.p1 + 3*self.p2
        C = -3*self.p0 + 3*self.p1
        a = 3*cross(A,B)
        b = 3*cross(A,C)
        c =   cross(B,C)
        r = b*b - 4*a*c
        if r >= 0 and a:
            r = math.sqrt(r)
            ip1 = (-b + r) / (2*a)
            ip2 = (-b - r) / (2*a)
            return min(ip1,ip2), max(ip1,ip2)
        return None, None


    def inflection_domain(self, t, flatness=0.25):
        """ Determine the domain around an inflection point
            where the curve is flat. """

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


    def angle(self):
        """ Compute angle betwenn (p0,p1) and (p2,p3) """

        dx0 = self.x1-self.x0
        dy0 = self.y1-self.y0
        dx1 = self.x3-self.x2
        dy1 = self.y3-self.y2
        angle = math.atan2(abs(dx0*dy1-dy0*dx1), dx0*dx1+dy0*dy1)
        return angle


    def flatten(self, flatness=0.25, angle=15):
        angle *= math.pi/180.0
        P = []
        while 1:
            dx = self.x1 - self.x0
            dy = self.y1 - self.y0
            norm = math.hypot(dx,dy)
            if not norm:
                break
            s3 = abs((self.x2-self.x0)*dy-(self.y2-self.y0)*dx)/norm
            t = 2*math.sqrt(flatness /(3*s3))
            if t > 1:
                break

            # Check angle is below tolerance
            for i in xrange(20):
                left, right = self.split(t)
                if left.angle() > angle:
                   t /= 1.5
                else:
                    break

            self.x0, self.y0 = right.x0, right.y0
            self.x1, self.y1 = right.x1, right.y1
            self.x2, self.y2 = right.x2, right.y2
            self.x3, self.y3 = right.x3, right.y3
            P.append((self.x0, self.y0))
        return P


    def flatten_brute_iterative(self, n=50):
        """ Brute force segmentation """
        P = [(self.x0,self.y0)]
        for i in xrange(1,n-1):
            t = i / float(n)
            P.append(self(t))
        P.append((self.x3,self.y3))
        return P

    def flatten_forward_iterative(self, n=50):
        """ Dumb segmentation """

        h = 1.0 / n;
        fph = 3 * (self.p1 - self.p0) * h
        fpphh = (6 * self.p0 - 12 * self.p1 + 6 * self.p2) * h * h
        fppphhh = (-6 * self.p0 + 18 * self.p1 - 18 * self.p2 + 6 * self.p3) * h * h * h
        P = [(self.x0,self.y0)]
        p = vec2(self.x0,self.y0)
        for i in range(1,n-1):
            p += fph + fpphh/2. + fppphhh/6.
            P.append((p.x,p.y))
            fph = fph + fpphh + fppphhh/2.
            fpphh = fpphh + fppphhh
        P.append((self.x3,self.y3))
        return P

    def flatten_forward_iterative_variable(self):
        """  """

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

        tmp = CubicBezier( (self.x0,self.y0), (self.x1,self.y1),
                           (self.x2,self.y2), (self.x3,self.y3))
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
                middle = self.segment(t1_plus,t2_minus)
                points += middle.flatten(flatness, angle)
                points.append( self(t2_minus) )

            if not t2_cross_end:
                points.append( self(t2_plus) )
                _,right = self.split(t2_plus)
                points += right.flatten(flatness, angle)

        points.append( (self.x3, self.y3) )
        return points

    def flatten_behdad_segment(self, flatness=0.125):
        """ """
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


# I = self.inflection_points()
# for t in self.inflection_points():
#     if t:
#         t_neg,t_pos = self.inflection_domain(t)
#         if 0 <= t <= 1:
#             P = self(t)
#             plt.scatter([P[0],], [P[1],], s=50,
#                         edgecolor='r', facecolor='None')
#             plt.scatter([P[0],], [P[1],], s=50, zorder=10,
#                         edgecolor='None', facecolor='white',alpha=.5)
#         if 0 <= t_neg <= 1:
#             P = self(t_neg)
#             plt.scatter([P[0],], [P[1],], s=10, color='k')
#         if 0 <= t_pos <= 1:
#             P = self(t_pos)
#             plt.scatter([P[0],], [P[1],], s=10, color='k')
