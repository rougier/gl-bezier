#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from behdad import Point, Bezier

np.random.seed(1)

n = 50000

# sampling precision
num_segs = 200
seg_ts = [i / float(num_segs) for i in range(0, num_segs + 1)]
chord_indices = ((0,3), (0,1), (2,3), (0, 2), (1,3), (1, 2))

def generate_test_data(n):
	numbers = np.random.randint(-1000,1000,(n,4,2))
	points = [[Point(x,y) for (x,y) in points] for points in numbers]
	beziers = [Bezier(*plist) for plist in points]
	return beziers, points, numbers

def calculate_bezier_lengths(beziers):
	n = len(beziers)
	lengths = np.zeros((n,1))
	for i, bez in enumerate(beziers):

		# Set up sampled length
		samples = [bez(t) for t in seg_ts]
		lengths[i] = sum((b-a).len() for a,b in zip(samples[1:], samples[:-1]))
	return lengths

def calculate_chord_lengths(points):
	n = len(points)
	chords = np.matrix(np.zeros((n,len(chord_indices))))
	for i, plist in enumerate(points):
		# Set up chord lengths
		chords[i] = [(plist[b]-plist[a]).len() for a,b in chord_indices]
	return chords

#train

beziers, points, _ = generate_test_data(n)
lengths = calculate_bezier_lengths(beziers)
chords = calculate_chord_lengths(points)
coeffs, residuals, _, _ = np.linalg.lstsq(chords, lengths)

print "Learned coefficients:"
print coeffs
print "Mean percentage error-squared-root (training):", (residuals[0] / n) ** .5 / (sum(lengths)[0] / n) * 100

lengths_new = np.array(chords * coeffs)
errors = (lengths_new - lengths) / lengths
print "Mean percentage error (training):", sum(abs(errors)) / len(errors) * 100

# test

beziers, points, _ = generate_test_data(n)
lengths = calculate_bezier_lengths(beziers)
chords = calculate_chord_lengths(points)

lengths_new = np.array(chords * coeffs)
errors = (lengths_new - lengths) / lengths
print "Mean percentage error (testing):", sum(abs(errors)) / len(errors) * 100
