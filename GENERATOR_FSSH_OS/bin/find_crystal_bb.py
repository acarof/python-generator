#!/usr/bin/env/python
""" Finds the crystal bounding box for a line in space."""
import math
import numpy as np
import sys


def abc_to_hmatrix(a, b, c, alpha, beta, gamma):
	""" Box vectors from box vector lengths and box vector angles (in degrees)."""
	alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
	result = np.zeros((3, 3))

	a = np.array((a, 0, 0))
	b = b * np.array((math.cos(gamma), math.sin(gamma), 0))
	bracket = (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
	c = c * np.array((math.cos(beta), bracket, math.sqrt(math.sin(beta) ** 2 - bracket ** 2)))

	result[:, 0] = a
	result[:, 1] = b
	result[:, 2] = c

	return result


def build_sphere_points(numpoints=1000):
	""" Generates a list of random points on a unit sphere."""
	random = np.random.normal(size=(numpoints, 3))
	norm = np.sqrt(np.sum(random * random, axis=1))
	return (random.T / norm).T


def build_circle_points(normal, numpoints=100):
	""" Generates equally spaced points on a unit circle around the origin given a normal vector."""
	base = np.identity(3)[np.argmax(np.abs(normal))]
	if np.linalg.norm(base - normal) < 10e-10:
		a, b = np.identity(3)[np.where(normal < 0.5)]
	else:
		a = np.cross(normal, base)
		b = np.cross(normal, a)
	a = a / np.linalg.norm(a)
	b = b / np.linalg.norm(b)
	ts = np.linspace(0, np.pi, numpoints, endpoint=False)
	return np.outer(np.cos(ts), a) + np.outer(np.sin(ts), b)



def find_crystal_bb( abc, start, length, vector, radius):

	hmat = abc_to_hmatrix(*abc)
	stepwidth = 0.5
	vector = vector / np.linalg.norm(vector)
	hinv = np.linalg.inv(hmat)

	# Find box
	intermediate = np.linspace(0., length, length/stepwidth + 1)
	sphere = build_sphere_points() * radius
	circle = build_circle_points(vector) * radius
	minfractionals = np.zeros(3)
	maxfractionals = np.zeros(3)
	for step in intermediate:
		circle_points = start + step * vector + circle
		fractional = hinv.dot(circle_points.T).T
		minfractionals = np.minimum(minfractionals, np.amin(fractional, axis=0))
		maxfractionals = np.maximum(maxfractionals, np.amax(fractional, axis=0))
	for step in (intermediate[0], intermediate[-1]):
		sphere_points = start + step * vector + sphere
		fractional = hinv.dot(sphere_points.T).T
		minfractionals = np.minimum(minfractionals, np.amin(fractional, axis=0))
		maxfractionals = np.maximum(maxfractionals, np.amax(fractional, axis=0))

	print 'Number of unit cells needed in direction of the three lattice vectors (assumes start point is in unit cell):'
	result = []
	coord_first = []
	for options in zip('abc', minfractionals, maxfractionals):
		print 'Vector %s needs to go from %f to %f' % options
		result.append( int(np.ceil(abs(options[2]))*np.copysign(1, options[2]) - np.ceil(abs(options[1]))*np.copysign(1, options[1])) )
		coord_first.append( int(np.ceil(abs(options[1]))))
	print 'Any fractional unit cells need to be fully included, i.e. -2.1 unit cells means 3 in negative direction and +2.1 means 3 in positive direction'
	print 'i.e.:   %s' % result
	print 'First molecule:', coord_first
	return result, coord_first


def _built_list_activated(coord_charge, result, nmol_unit):
		result =  (int((float(coord_charge[0]) - 1) * float(result[1] * result[2]) + \
				 (float(coord_charge[1]) - 1) * float(result[2]) + \
				 float(coord_charge[2])) - 1) * nmol_unit + 1
		print 'Charged molecule:' , result
		return result


def find_molecules(coord_first, size_crystal, length, vector, radius_aom = 0.0, psf_file = None, xyz_file = None, nmol_unit = -1 ):

		import MDAnalysis as mda
		u = mda.Universe(psf_file, xyz_file)

		def detect_molecules(universe):
			""" Returns a list of atom groups that are molecules as identified by their topology."""
			molecules = []
			bonds = [_.indices for _ in u.bonds]
			atommask = [True] * len(u.atoms)
			while True in atommask:
				first = atommask.index(True)
				found = True
				mol = [first]
				while found:
					found = False
					for bond in bonds:
						if bond[0] in mol and not bond[1] in mol:
							mol.append(bond[1])
							found = True
						if bond[1] in mol and not bond[0] in mol:
							mol.append(bond[0])
							found = True
				for atom in mol:
					atommask[atom] = False
				molecules.append(u.atoms[mol])
			return molecules
		molecules = detect_molecules(u)
		coms = [_.center_of_mass() for _ in molecules]


		stepwidth = 0.5
		vector = vector / np.linalg.norm(vector)
		intermediate = np.linspace(0., length, length / stepwidth + 1)

		first_mol =  _built_list_activated(coord_first, size_crystal, nmol_unit)
		start = coms[first_mol]
		active = [False] * len(molecules)
		for step in intermediate:
			this = start + step * vector
			for molecule in range(len(molecules)):
				delta = coms[molecule] - this
				if np.linalg.norm(delta) < radius_aom and np.dot(delta, vector) >= 0:
					active[molecule] = True

		active = [ (_[0]) for _ in enumerate(active) if _[1] == True]
		print 'Active molecule indices (0-based): ' + ' '.join(map(str, active))
		return (first_mol+1, [ x +1 for x in active])
