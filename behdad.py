#!/usr/bin/python
# -*- coding: utf-8 -*-

#
# Copyright 2012,2013 Google, Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Google Author(s): Behdad Esfahbod
#
# Ported from GLyphy
#

import math


# returns tan (2 * atan (d))
def tan2atan(d): return 2 * d / (1. - d*d)

# returns sin (2 * atan (d))
def sin2atan(d): return 2 * d / (1. + d*d)


class Point:

    def __init__ (self, x, y):
        self.x = float (x)
        self.y = float (y)

    def __repr__ (self):
        return "Point(%g,%g)" % (self.x, self.y)

    def __eq__ (self, other):
        return self.x == other.x and self.y == other.y

    def __add__ (self, other):
        return Point (self.x + other.dx, self.y + other.dy)

    def __sub__ (self, other):
        if isinstance (other, Vector):
            return Point (self.x - other.dx, self.y - other.dy)
        else:
            return Vector (self.x - other.x, self.y - other.y)

    def midpoint(self, p):
        return self + (p - self)/2.

    def lerp (self, a, other):
        # The following two cases are special-cased to get better floating
        # point stability.  We require that points that are the same be
        # bit-equal.
        if a == 0  : return Point (self.x, self.y)
        if a == 1.0: return Point (other.x, other.y)
        return Point ((1-a) * self.x + a * other.x, (1-a) * self.y + a * other.y)

class Vector:

    def __init__ (self, dx, dy):
        self.dx = float (dx)
        self.dy = float (dy)

    def __repr__ (self):
        return "Vector(%g,%g)" % (self.dx, self.dy)

    def __eq__ (self, other):
        return self.dx == other.dx and self.dy == other.dy

    def __add__ (self, other):
        return Vector (self.dx + other.dx, self.dy + other.dy)

    def __sub__ (self, other):
        return Vector (self.dx - other.dx, self.dy - other.dy)

    def __mul__ (self, other):
        if isinstance (other, Vector):
            return self.dx * other.dx + self.dy * other.dy
        else:
            return Vector (self.dx * other, self.dy * other)

    def __div__ (self, s):
        return Vector(self.dx/s, self.dy/s)

    def len (self):
        return math.hypot (self.dx, self.dy)

    def angle (self):
        return math.atan2 (self.dy, self.dx)

    def perpendicular (self):
        return Vector (-self.dy, self.dx)

    def normalized (self):
        d = self.len ()
        if not d:
            d = 1.
        return Vector (self.dx / d, self.dy / d)

    def rebase (self, bx, by = None):
        if by is None:
            by = bx.perpendicular ()

        return Vector (self * bx, self * by)

class Arc:

    def __init__ (self, p0, p1, d):
        self.p0 = p0
        self.p1 = p1
        self.d = float (d)

    def to_conventional(self):
        radius = self.radius()
        center = self.center()
        angle0 = (self.p0 - center).angle ()
        angle1 = (self.p1 - center).angle ()
        negative = self.d < 0
        return (center.x,center.y),radius,angle0,angle1,negative,self

    def radius(self):
        return abs((self.p1-self.p0).len() / (2. * sin2atan (self.d)))

    def center(self):
        return (self.p0.midpoint(self.p1)) + \
               (self.p1-self.p0).perpendicular()/(2. * tan2atan(self.d))

    def tangents(self):
        dp = (self.p1 - self.p0) * .5
        pp = dp.perpendicular () * -tan2atan (self.d)
        return dp + pp, dp - pp

    def wedge_contains_point(self,p):
        # TODO this doesn't handle fabs(d) > 1.
        t0,t1 = self.tangents ()
        return (p - self.p0) * t0  >= 0 and (p - self.p1) * t1 <= 0

    def __repr__ (self):
        return "Arc(%s,%s,%g)" % (self.p0, self.p1, self.d)

    def distance_to_point (self, p):
        if self.wedge_contains_point (p):
            if abs (self.d) < 1e-5:
                # Return distance to line
                pp = (self.p1 - self.p0).perpendicular ().normalized ()
                return abs ((p - self.p) * pp)
            else:
                c = self.center ()
                r = self.radius ()
                return abs ((p - c).len () - r)

        return min ((p - self.p0).len (), (p - self.p1).len ())

    @staticmethod
    def from_three_points (p0, p1, pm, complement) :
        if p0 == pm or p1 == pm:
            d = 0.
        else:
            angle = ((p1-pm).angle () - (p0-pm).angle ()) / 2
            if not complement:
                angle -= math.pi / 2
            d = math.tan (angle)
        return Arc (p0, p1, d)

    def approximate_bezier (self):

        d = self.d

        dp = self.p1 - self.p0
        pp = dp.perpendicular ()

        error = dp.len () * (abs (d) ** 5) / (54 * (1 + d*d))

        p0s = self.p0 + dp * ((1 - d*d) / 3) - pp * (2 * d / 3)
        p1s = self.p1 - dp * ((1 - d*d) / 3) - pp * (2 * d / 3)

        return Bezier (self.p0, p0s, p1s, self.p1), error

class Bezier:

    def __init__ (self, p0, p1, p2, p3):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def __repr__ (self):
        return "Bezier(%s,%s,%s,%s)" % (self.p0, self.p1, self.p2, self.p3)

    def split (self, t):
        p0, p1, p2, p3 = self.p0, self.p1, self.p2, self.p3

        p01 = p0.lerp (t, p1)
        p12 = p1.lerp (t, p2)
        p23 = p2.lerp (t, p3)
        p012 = p01.lerp (t, p12)
        p123 = p12.lerp (t, p23)
        p0123 = p012.lerp (t, p123)
        return Bezier (p0, p01, p012, p0123), Bezier (p0123, p123, p23, p3)

    def segment (self, t0, t1):
        p0, p1, p2, p3 = self.p0, self.p1, self.p2, self.p3

        p01 = p0.lerp (t0, p1)
        p12 = p1.lerp (t0, p2)
        p23 = p2.lerp (t0, p3)
        p012 = p01.lerp (t0, p12)
        p123 = p12.lerp (t0, p23)
        p0123 = p012.lerp (t0, p123)

        q01 = p0.lerp (t1, p1)
        q12 = p1.lerp (t1, p2)
        q23 = p2.lerp (t1, p3)
        q012 = q01.lerp (t1, q12)
        q123 = q12.lerp (t1, q23)
        q0123 = q012.lerp (t1, q123)

        return Bezier (p0123,
                   p0123 + (p123 - p0123) * ((t1 - t0) / (1 - t0)),
                   q0123 + (q012 - q0123) * ((t1 - t0) / t1),
                   q0123)


class MaxDeviationApproximatorExact:
    
    def __call__ (self, d0, d1):
        """Returns 3 max(abs(d₀ t (1-t)² + d₁ t² (1-t)) for 0≤t≤1."""

        candidates = [0, 1]
        if d0 == d1:
            candidates.append (.5)
        else:
            delta = d0*d0 - d0*d1 + d1*d1
            t2 = 1. / (3 * (d0 - d1))
            t0 = (2 * d0 - d1) * t2
            if delta == 0:
                candidates.append (t0)
            elif delta > 0:
                # This code can be optimized to avoid the sqrt if
                # the solution is not feasible (ie. lies outside
                # (0,1)).  I have implemented that in
                # cairo-spline.c:_cairo_spline_bound().  Can be
                # reused here.
                t1 = math.sqrt (delta) * t2
                candidates.append (t0 - t1)
                candidates.append (t0 + t1)

        e = 0
        for t in candidates:
            if t < 0. or t > 1.:
                continue
            ee = abs (3 * t * (1-t) * (d0 * (1-t) + d1 * t))
            e = max (e, ee)
        
        return e


class ArcBezierErrorApproximatorBehdad:

    def __init__ (self, MaxDeviationApproximator):
        self.MaxDeviationApproximator = MaxDeviationApproximator

    def __call__ (self, b, a):
        """Returns upper bound for error between Bezier b and arc a."""
        b0 = b
        del b

        assert b0.p0 == a.p0
        assert b0.p3 == a.p1

        b1, ea = a.approximate_bezier ()

        assert b0.p0 == b1.p0
        assert b0.p3 == b1.p3

        v0 = b1.p1 - b0.p1
        v1 = b1.p2 - b0.p2

        b = (b0.p3 - b0.p0).normalized ()
        v0 = v0.rebase (b)
        v1 = v1.rebase (b)

        v = Vector (self.MaxDeviationApproximator (v0.dx, v1.dx),
                self.MaxDeviationApproximator (v0.dy, v1.dy))

        # Edge cases: If d*d is too close to being 1 default to a weak bound.
        if abs (a.d * a.d - 1) < 1e-4:
            return ea + v.len ()

        # We made sure that abs (a.d) != 1
        tan_half_alpha = 2 * abs (a.d) / (1 - a.d*a.d)

        # If v.dy == 0, perturb just a bit.
        if abs (v.dy) < 1e-6:
             # TODO Figure this one out.
             v.dy = 1e-6

        tan_v = v.dx / v.dy

        # TODO Double check and simplify these checks
        if abs (a.d) < 1e-6 or tan_half_alpha < 0 or \
           (-tan_half_alpha <= tan_v and tan_v <= tan_half_alpha):
            return ea + v.len ()

        c2 = (b1.p3 - b1.p0).len () / 2
        r = c2 * (a.d * a.d + 1) / (2 * abs (a.d))

        eb = Vector (c2 / tan_half_alpha + v.dy, c2 + v.dx).len () - r
        assert eb >= 0

        return ea + eb

class ArcBezierApproximatorMidpointTwoPart:

    def __init__ (self, ArcBezierErrorApproximator):
        self.ArcBezierErrorApproximator = ArcBezierErrorApproximator

    def __call__ (self, b, mid_t = .5):
        """Approximates Bezier b with arc passing through mid_t.
        Returns arc,error."""

        pair = b.split (mid_t)
        m = pair[1].p0

        a0 = Arc.from_three_points (b.p0, m, b.p3, True)
        a1 = Arc.from_three_points (m, b.p3, b.p0, True)

        e0 = self.ArcBezierErrorApproximator (pair[0], a0)
        e1 = self.ArcBezierErrorApproximator (pair[1], a1)
        error = max (e0, e1)

        return Arc.from_three_points (b.p0, b.p3, m, False), error

class ArcsBezierApproximatorSpringSystem:

    def __init__ (self, ArcBezierApproximator):
        self.ArcBezierApproximator = ArcBezierApproximator

    def __call__ (self, b, tolerance, max_segments = 100):

        tolerance = float (tolerance)

        # Technically speaking we can bsearch for n.
        for n in range (1, max_segments + 1):

            ts = [float (i) / n for i in range (n)]
            ts.append (1.0) # Do this out of the loop to get real 1.0, not .9999999999999998!

            arcs, errors = self.__calc_arcs (b, ts)

            if min (errors) <= tolerance:
                ts, arcs, errors = self.__jiggle (b, tolerance, ts, arcs, errors)

            if max (errors) <= tolerance:
                break

        return arcs, max (errors)


    def __calc_arcs (self, b, ts):

        arcs = []
        errors = []
        for i in range (len (ts) - 1):
            segment = b.segment (ts[i], ts[i + 1])
            arc, error = self.ArcBezierApproximator (segment)
            arcs.append (arc)
            errors.append (error)

        return arcs, errors

    def __jiggle (self, b, tolerance, ts, arcs, errors):

        n = len (ts) - 1
        max_jiggle = int (math.log (n) / math.log (2) + 1)

        for s in range (max_jiggle):

            total = 0.
            for i in range (n):
                l = ts[i + 1] - ts[i]
                k_inv = l * (errors[i] ** -.3)
                total += k_inv
                errors[i] = k_inv

            for i in range (n):
                k_inv = errors[i]
                l = k_inv / total
                ts[i + 1] = ts[i] + l
            ts[n] = 1.0 # Do this to get real 1.0, not .9999999999999998!

            arcs, errors = self.__calc_arcs (b, ts)

            if max (errors) <= tolerance or 2 * min (errors) - max (errors) > tolerance:
                break

        if s == max_jiggle:
            print "JIGGLE OVERFLOW n %d s %d" % (n, s)

        return ts, arcs, errors


if __name__ == "__main__":

    b = Bezier (Point(0,0), Point(0,1), Point(1,1), Point(1,0))

    errfunc = ArcBezierErrorApproximatorBehdad (MaxDeviationApproximatorExact ())
    apprfunc = ArcBezierApproximatorMidpointTwoPart (errfunc)
    splinefunc = ArcsBezierApproximatorSpringSystem (apprfunc)

    #print apprfunc (b)

    tolerance = .001
    arcs, error = splinefunc (b, tolerance)
    print len (arcs), error

    for arc in arcs:
        print arc.to_conventional()
