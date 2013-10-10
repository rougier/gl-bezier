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
import time
import numpy as np
from cubic_bezier import CubicBezier


# Measure errors on 10,000 curves
np.random.seed(1)
n = 1000
curves = np.random.randint(100,700,(n,4,2))
flatness = 0.125
angle = 15
n1, n2 = 25, 50

# Forward method, n=25
# -------------------------------------
t_start = time.clock()
for i in range(n):
    C = CubicBezier(*curves[i])
    C.flatten_forward_iterative(n=n1)
t_end = time.clock()
t_mean = 1000 * (t_end - t_start)/float(n)
print "Forward iterative (n=%d): %.3f ms" % (n1, t_mean)

# Forward method, n=50
# -------------------------------------
t_start = time.clock()
for i in range(n):
    C = CubicBezier(*curves[i])
    C.flatten_forward_iterative(n=n2)
t_end = time.clock()
t_mean = 1000 * (t_end - t_start)/float(n)
print "Forward iterative (n=%d): %.3f ms" % (n2, t_mean)


# Smart iterative
# -------------------------------------
t_start = time.clock()
for i in range(n):
    C = CubicBezier(*curves[i])
    C.flatten_iterative(flatness=flatness, angle=angle)
t_end = time.clock()
t_mean = 1000 * (t_end - t_start)/float(n)
print "Smart iterative: %.3f ms" % (t_mean)

# Recursive
# -------------------------------------
t_start = time.clock()
for i in range(n):
    C = CubicBezier(*curves[i])
    C.flatten_recursive(flatness=flatness, angle=angle)
t_end = time.clock()
t_mean = 1000 * (t_end - t_start)/float(n)
print "Recursive: %.3f ms" % (t_mean)

# Arcs
# -------------------------------------
t_start = time.clock()
for i in range(n):
    C = CubicBezier(*curves[i])
    C.flatten_behdad_arc(flatness=flatness)
t_end = time.clock()
t_mean = 1000 * (t_end - t_start)/float(n)
print "Arc iterative: %.3f ms" % (t_mean)

