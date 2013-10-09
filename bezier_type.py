# Adapted from
# "http://jgt.akpeters.com/papers/Vincent02/BezType.html"
#
# Implementation of the algorithm described in:
#
#      Stephen Vincent.
#      Fast Detection of the Geometric Form of Two-Dimensional Cubic B&eacute;zier Curves.
#      Journal of Graphics Tools, 7(3):43-51, 2002
#
# See the paper for discussion of the algorithm.
#
import math
from vec2 import *

def enum(**enums):
    return type('Enum', (), enums)

CubicBezierType = enum(arch = 'Arch',
                       single = 'Single inflection',
                       double = 'Double inflection',
                       cusp   = 'Cusp',
                       loop   = 'Loop',
                       line   = 'Line')
epsilon = 1e-10


def sign(a):
    if a > epsilon:
        return 1
    elif a < -epsilon:
        return -1
    return 0


def evaluate_point(t, p0, p1, p2, p3):
    # Evaluates the point on a cubic Bezier whose control points are in the
    # array cpt at parameter value t and returns the resuls in pt

    return p0 * (1 - t)**3     + \
           p1 * 3*t * (1-t)**2 + \
           p2 * 3*t**2 *(1-t)  + \
           p3 * t**3


def cusp_value(det_012, det_013, det_023):
    """
    Distinguish between loop, cusp, and 2 inflection point curves.
    """

    a = 3*det_012 + det_023 - 2*det_013
    b = -3*det_012 + det_013
    c = det_012
    d = b*b - 4*a*c

    # Rather than test against epsilon, you could test against sqrt(d)/2*a for
    # a better approximation to a cusp.
    if d > epsilon:
        return CubicBezierType.double
    else:
        if d < -epsilon:
            return CubicBezierType.loop
        else:
            return CubicBezierType.cusp


def point_relative_to_line(p0, p1, p2):
    """
    Determine the position of a point p2 with respect to a line
    defined by p0 and p1. 'Find' the closest point on the line to p2.
    If it lies before p0 return -1 , if it lies after p1 return 1 :
    otherwise return 0
    """
    a = -(p1.x - p0.x) * (p2.x - p0.x)
    b = (p1.y - p0.y) * (p2.y - p0.y)
    if (a - b) >= 0:
        return -1
    else:
        a = - ( p0.x - p1.x ) * ( p2.x - p1.x )
        b = ( p0.y - p1.y ) * ( p2.y - p1.y )
        if ( a - b ) >= 0.0:
            return 1 
        else:
            return 0


def cubic_bezier_type(p0,p1,p2,p3):

    det_012 = det(p0, p1, p2)
    det_123 = det(p1, p2, p3)
    det_013 = det(p0, p1, p3)
    det_023 = det(p0, p2, p3)

    sign_012 = sign(det_012)
    sign_123 = sign(det_123)
    sign_013 = sign(det_013)
    sign_023 = sign(det_023)

    # First address the cases where 3 or more consecutive
    # control points are colinear
    # ------------------------------------------------------
    # Case E : all 4 points are colinear.
    if sign_012 == 0 and sign_123 == 0:
        # Could test for a single point here if necessary
        if sign_013 == 0:
            return CubicBezierType.straight_line
        # Points 1 and 2 coincident
        else:
            return CubicBezierType.arch

    # Case F : first 3 control points are colinear
    elif sign_012 == 0:
        if sign_013 == sign_123:
            return CubicBezierType.arch
        else:
            return CubicBezierType.single

    # Case F : second 3 control points are colinear
    elif sign_123 == 0:
        if sign_023 == sign_012:
            return CubicBezierType.arch
        else:
            return CubicBezierType.single

    # Case G : points 0,1,3 are colinear
    elif sign_013 == 0:
        k = point_relative_to_line (p0, p3, p1)
        if k == -1:
            return CubicBezierType.arch
        else:
            if k == 0:
                return CubicBezierType.single
            else:
                return cusp_value(det_012, det_013, det_023 )

    # Case G : points 0,2,3 are colinear
    elif sign_023 == 0:
        k = point_relative_to_line(p0, p3, p2)
        if k == 1:
            return CubicBezierType.arch
        else:
            if k == 0:
                return CubicBezierType.single
            else:
                return cusp_value(det_012, det_013, det_023)

    # On to the more interesting stuff. At this point it's
    # known that no 3 of the control points are colinear
    # ------------------------------------------------------

    # Case A : the control points zig-zag
    elif sign_012 != sign_123:
        return CubicBezierType.single

    # Case B : Convex control polygon
    elif sign_012 == sign_013 and sign_012 == sign_023:
        return CubicBezierType.arch

    # Case C : Self-intersecting control polygon
    elif sign_012 != sign_013 and sign_012 != sign_023:
        return cusp_value(det_012, det_013, det_023 )

    # Case D : Concave control polygon
    else:
        t = det_013 / (det_013 - det_023)
        pt = evaluate_point (t, p0,p1,p2,p3)
        if point_relative_to_line( p0, p3, pt) == 0:
            return cusp_value(det_012, det_013, det_023)
        else:
            return CubicBezierType.arch
