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

# returns cos (2 * atan (d))
def cos2atan(d): return (1. - d*d) / (1. + d*d)


class Point (object):

    __slots__ = ("x", "y")

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

class Vector (object):

    __slots__ = ("dx", "dy")

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

    def len2 (self):
        return self * self

    def angle (self):
        return math.atan2 (self.dy, self.dx)

    def ortho (self):
        return Vector (-self.dy, self.dx)

    def normalized (self):
        d = self.len ()
        if not d:
            d = 1.
        return Vector (self.dx / d, self.dy / d)

    def rebase (self, bx, by = None):
        if by is None:
            by = bx.ortho ()

        return Vector (self * bx, self * by)

class Arc (object):

    __slots__ = ("p0", "p1", "d")

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
        return (center.x,center.y),radius,angle0,angle1,negative

    def radius(self):
        return abs((self.p1-self.p0).len() / (2. * sin2atan (self.d)))

    def center(self):
        return (self.p0.midpoint(self.p1)) + \
               (self.p1-self.p0).ortho()/(2. * tan2atan(self.d))

    def tangents(self):
        x = self.p1.x - self.p0.x
        y = self.p1.y - self.p0.y

        d = (1 - self.d*self.d) * .5
        dpx = x * d
        dpy = y * d

        d = -self.d
        ppx = -y * d
        ppy = x * d

        return Vector(dpx+ppx, dpy+ppy), Vector(dpx-ppx, dpy-ppy)

        # Here's the original, slower code:
        dp = (self.p1 - self.p0) * .5
        pp = dp.ortho () * -sin2atan (self.d)
        dp = dp * cos2atan (self.d)
        return dp + pp, dp - pp


    def wedge_contains_point(self,p):
        t0,t1 = self.tangents ()
        if abs (self.d) <= 1.:
            return (p - self.p0) * t0  >= 0 and (p - self.p1) * t1 <= 0
        else:
            return (p - self.p0) * t0  >= 0 or (p - self.p1) * t1 <= 0

    def __repr__ (self):
        return "Arc(%s,%s,%g)" % (self.p0, self.p1, self.d)

    def distance_to_point (self, p):
        if self.wedge_contains_point (p):
            if abs (self.d) < 1e-5:
                # Return distance to line
                pp = (self.p1 - self.p0).ortho ().normalized ()
                return abs ((p - self.p0) * pp)
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
        pp = dp.ortho ()

        error = dp.len () * (abs (d) ** 5) / (54 * (1 + d*d))

        dp *= (1 - d*d) / 3
        pp *= 2 * d / 3

        p0s = self.p0 + dp - pp
        p1s = self.p1 - dp - pp

        return Bezier (self.p0, p0s, p1s, self.p1), error

class Bezier (object):

    __slots__ = ("p0", "p1", "p2", "p3")

    def __init__ (self, p0, p1, p2, p3):
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def __repr__ (self):
        return "Bezier(%s,%s,%s,%s)" % (self.p0, self.p1, self.p2, self.p3)

    def __call__ (self, t):
        p0, p1, p2, p3 = self.p0, self.p1, self.p2, self.p3

        p01 = p0.lerp (t, p1)
        p12 = p1.lerp (t, p2)
        p23 = p2.lerp (t, p3)
        p012 = p01.lerp (t, p12)
        p123 = p12.lerp (t, p23)
        p0123 = p012.lerp (t, p123)

        return p0123

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

        # Edge cases: If d*d is too large default to a weak bound.
        if a.d * a.d > 1. - 1e-4:
            return ea + v.len ()

        # If the wedge doesn't contain control points, default to weak bound.
        if not a.wedge_contains_point (b0.p1) or not a.wedge_contains_point (b0.p2):
            return ea + v.len ()

        # If straight line, return the max ortho deviation.
        if abs (a.d) < 1e-6:
            return ea + v.dy

        # We made sure that abs(a.d) < 1
        tan_half_alpha = abs (tan2atan (a.d))

        tan_v = v.dx / v.dy

        if abs (tan_v) <= tan_half_alpha:
            return ea + v.len ()

        c2 = (a.p1 - a.p0).len () * .5
        r = a.radius ()

        eb = Vector (c2 + v.dx, c2 / tan_half_alpha + v.dy).len () - r
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

        min_segments = 1
	ts = None

        # Technically speaking we can bsearch for n.
        for n in range (min_segments, max_segments + 1):

	    if not ts:
                ts = [float (i) / n for i in range (n)]
	    else:
		assert n == len(ts)

		new_ts = []
		for i in range (len (ts) - 1):
		    new_ts.append (ts[i] + (ts[i + 1] - ts[i]) * (n - 1 - i) / n)
		new_ts.insert(0, 0.0)
		ts = new_ts
            ts.append (1.0) # Do this out of the loop to get real 1.0, not .9999999999999998!

            arcs, errors = self.__calc_arcs (b, ts)

            if min (errors) <= tolerance:
                ts, arcs, errors = self.__jiggle (b, tolerance, ts, arcs, errors)

            if max (errors) <= tolerance:
                break

        return arcs, max (errors), ts


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

            if 0.0 in errors:
                bias = sum (errors) / len (errors)
                if bias == 0.0:
                    # All errors are zero!
                    return ts, arcs, errors
                errors = [e + bias for e in errors]

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

    test_beziers = [
        ((8, 21), (5, 17), (13, 27), (33, 48)),
        ((38, 23), (44, 23), (13, 25), (10, 2)),
        ((33, 8), (41, 30), (38, 0), (32, 32)),
        ((42, 23), (36, 24), (19, 46), (44, 36)),
        ((4, 43), (4, 18), (4, 31), (4, 36)),
        ((39, 24), (21, 24), (23, 24), (3, 24)),
        ((41, 28), (2, 11), (25, 11), (18, 15)),
        ((48, 13), (2, 47), (37, 43), (23, 35)),
        ((11, 17), (44, 49), (4, 30), (1, 44)),
        ((18, 4), (38, 23), (30, 21), (38, 35)),
        ((24, 46), (1, 40), (20, 20), (19, 4)),
        ((11, 8), (18, 17), (15, 4), (10, 45)),
        ((44, 6), (23, 6), (11, 33), (1, 48)),
        ((13, 47), (27, 32), (27, 11), (44, 4)),
        ((33, 14), (20, 20), (8, 32), (4, 27)),
        ((170, 660), (170, 660), (660, 510), (540, 600)),
        ((381, 278), (376, 354), (457, 568), (352, 590)),
        ((498, 662), (680, 315), (603, 578), (186, 241)),
        ((413, 615), (597, 143), (688, 126), (436, 397)),
        ((566, 115), (164, 296), (125, 467), (571, 382)),
        ((196, 369), (110, 379), (252, 576), (516, 240)),
        ((293, 672), (589, 252), (438, 230), (592, 366)),
        ((506, 550), (576, 349), (590, 451), (665, 541)),
        ((569, 372), (545, 460), (642, 236), (268, 258)),
        ((484, 541), (286, 301), (432, 384), (409, 369)),
        ((190, 533), (352, 538), (685, 298), (647, 326)),
        ((472, 371), (140, 106), (232, 152), (403, 346)),
        ((411, 163), (213, 676), (422, 586), (213, 393)),
        ((381, 278), (376, 354), (457, 568), (352, 590)),
        ((566, 115), (164, 296), (125, 467), (571, 382)),
        ((122, 226), (379, 481), (456, 255), (413, 695)),
        ((266, 388), (518, 379), (187, 343), (448, 686)),
        ((290, 402), (251, 421), (469, 353), (561, 664)),
        ((196, 369), (110, 379), (252, 576), (516, 240)),
        ((249, 481), (359, 621), (485, 480), (224, 272)),
        ((136, 467), (207, 304), (470, 339), (528, 535)),
        ((473, 323), (409, 458), (278, 231), (482, 173)),
        ((513, 166), (420, 610), (201, 485), (305, 404)),
        ((634, 433), (364, 213), (571, 438), (185, 122)),
        ((381, 278), (376, 354), (457, 568), (352, 590)),
        ((498, 662), (680, 315), (603, 578), (186, 241)),
        ((493, 107), (419, 634), (413, 613), (416, 309)),
        ((163, 647), (577, 645), (568, 695), (148, 388)),
        ((617, 107), (562, 434), (645, 639), (177, 130)),
        ((187, 591), (289, 460), (506, 321), (403, 129)),
        ((117, 407), (148, 304), (463, 343), (570, 263)),
        ((124, 551), (242, 447), (514, 530), (612, 505)),
        ((536, 508), (405, 538), (269, 436), (170, 473)),
    ]

    errfunc = ArcBezierErrorApproximatorBehdad (MaxDeviationApproximatorExact ())
    apprfunc = ArcBezierApproximatorMidpointTwoPart (errfunc)
    splinefunc = ArcsBezierApproximatorSpringSystem (apprfunc)

    tolerance = .125

    for points in test_beziers:
        b = Bezier (*(Point(x,y) for x,y in points))
        arcs, error, ts = splinefunc (b, tolerance)
        print len (arcs), error, ts
        N = 1000
        for i in range (N + 1):
            t = float (i) / N
            p = b(t)
            d = min (arc.distance_to_point (p) for arc in arcs)
            if d > tolerance:
                print t, p, d
                print arcs
                assert 0
