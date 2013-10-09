#!/usr/bin/python
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
import math

class vec2(list):

    def __init__ (self, x=0, y=0):
        self.x = float (x)
        self.y = float (y)

    def __repr__ (self):
        return "vec2(%g,%g)" % (self.x, self.y)

    def __eq__ (self, other):
        return self.x == other.x and self.y == other.y

    def __add__ (self, other):
        return vec2(self.x+other.x, self.y+other.y)

    def __iadd__ (self, other):
        self.x += other.x
        self.y += other.y
        return self

    def __getitem__ (self, index):
        if   index == 0: return self.x
        elif index == 1: return self.y
        raise IndexError

    def __setitem__ (self, index, value):
        if   index == 0: self.x = value
        elif index == 1: self.y = value
        else: raise IndexError

    def __len__ (self):
        return 2

    def __iter__(self):
        for v in (self.x, self.y):
            yield v

    def __sub__ (self, other):
        return vec2(self.x-other.x, self.y-other.y)

    def __isub__ (self, other):
        self.x -= other.x
        self.y -= other.y
        return self

    def __neg__ (self):
        return vec2(-self.x, -self.y)

    def __mul__ (self, value):
        return vec2(self.x*value, self.y*value)

    def __rmul__ (self, value):
        return vec2(self.x*value, self.y*value)

    def __div__ (self, value):
        return vec2(self.x/value, self.y/value)

    def __rdiv__ (self, value):
        return vec2(self.x/value, self.y/value)

    def length(self):
        return math.hypot(self.x, self.y)

    def angle(self):
        return math.atan2 (self.y, self.x)

    def ortho(self):
        return vec2(-self.y, self.x)

    def normalized(self):
        length = math.hypot(self.x, self.y)
        if not length:
            return vec2(self.x, self.y)
        return vec2(self.x/length, self.y /length)

    def rebase (self, bx, by = None):
        if by is None:
            by = bx.ortho()
        return vec2(self*bx,self*by)


def normalized(v):
    return vec2(v[0],v[1]).normalized()

def length(v):
    return vec2(v[0],v[1]).length()

def dot(v0, v1):
    return v0[0]*v1[0] + v0[1]*v1[1]

def cross(v0, v1):
    return v0[0]*v1[1] - v0[1]*v1[0] 

def mix(v0, v1, a):
    if a == 0.0: return vec2(v0[0], v0[1])
    if a == 1.0: return vec2(v1[0], v1[1])
    return vec2 ((1-a)*v0[0] + a*v1[0],
                 (1-a)*v0[1] + a*v1[1])

def middle(v0, v1):
    return mix(v0,v1,0.5)

def angle(v0,v1):
    return math.atan2(abs(cross(v0,v1)),dot(v0,v1))

def det(p0, p1, p2):
    return (p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0])
