#!/usr/bin/env python

# Copyright 2012 by Silicos-it, a division of Imacosi BVBA
# Copyright 2020- by UAMC, the Medicinal Chemistry department of the University
# of Antwerp

__all__ = ['SpectrophoreCalculator']
__version__ = "1.0.1"

# Numba support
from numba import njit

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

# Numpy and SciPy
import numpy as np
import scipy
import scipy.linalg

# Math
import math



@njit()
def rotate(PARMS, COORD, ANGLES, ROTMAT, ENERGY, MINENERGY, BOXES, RESOLUTION):

    count = 0
    for step in ANGLES:

        # Rotate
        ca = math.cos(step[0])
        cb = math.cos(step[1])
        cc = math.cos(step[2])
        sa = math.sin(step[0])
        sb = math.sin(step[1])
        sc = math.sin(step[2])
        casb = ca*sb
        sasb = sa*sb
        cacc = ca*cc
        ROTMAT[0,0] = ca*cb
        ROTMAT[0,1] = casb*sc - sa*cc
        ROTMAT[0,2] = casb*cc + sa*sc
        ROTMAT[1,0] = sa*cb
        ROTMAT[1,1] = sasb*sc + cacc
        ROTMAT[1,2] = sasb*cc - ca*sc
        ROTMAT[2,0] = -sb
        ROTMAT[2,1] = cb*sc
        ROTMAT[2,2] = cacc
        for atom in range(len(COORD)):
            for i in range(3):
                for j in range(3):
                    COORD[atom][i+3] = 0.0
                    COORD[atom][i+3] += ROTMAT[i][j] * COORD[atom][j]


        # Update outer limits of molecule, taking into account atom radius and resolution
        for i in range(3):
            COORD[:,i+6] = COORD[:,i+3] - COORD[:,12] - RESOLUTION
            COORD[:,i+9] = COORD[:,i+3] + COORD[:,12] + RESOLUTION

        # Adopt box to fit around molecule
        mx = COORD[0,6]
        my = COORD[0,7]
        mz = COORD[0,8]
        px = COORD[0,6]
        py = COORD[0,7]
        pz = COORD[0,8]
        for i in range(1,len(COORD)):
            if COORD[i,6] < mx: mx = COORD[i,6]
            if COORD[i,7] < my: my = COORD[i,7]
            if COORD[i,8] < mz: mz = COORD[i,8]
        for i in range(1,len(COORD)):
            if COORD[i,6] > px: px = COORD[i,6]
            if COORD[i,7] > py: py = COORD[i,7]
            if COORD[i,8] > pz: pz = COORD[i,8]
        hx = (mx + px) / 2.0
        hy = (my + py) / 2.0
        hz = (mz + pz) / 2.0
        BOXES[0,-3:] =  hx,my,pz
        BOXES[1,-3:] =  px,hy,pz
        BOXES[2,-3:] =  hx,py,pz
        BOXES[3,-3:] =  mx,hy,pz
        BOXES[4,-3:] =  mx,my,hz
        BOXES[5,-3:] =  px,my,hz
        BOXES[6,-3:] =  mx,py,hz
        BOXES[7,-3:] =  px,py,hz
        BOXES[8,-3:] =  px,hy,mz
        BOXES[9,-3:] =  hx,my,mz
        BOXES[10,-3:] = mx,hy,mz
        BOXES[11,-3:] = hx,py,mz

        # Calculate energies
        # Empty the energy arrays
        for i in range(len(ENERGY)): ENERGY[i] = 0.0

        # Loop over each boxpoint (12 points)
        for boxPoint in range(PARMS[8]):

            # Loop over each atom
            for atom in range(len(COORD)):

                # Distance between boxpoint and atom
                d = math.sqrt(np.sum((BOXES[boxPoint][-3:] - COORD[atom][3:6])**2))

                # Loop over each probe
                for probe in range(PARMS[6]):

                    # Loop over each property (4 properties)
                    for prop in range(PARMS[7]):
                        index = PARMS[6]*prop+probe
                        ENERGY[index] += (COORD[atom][prop+13] * BOXES[boxPoint][probe]) / d

        if count == 0:
            MINENERGY = ENERGY
            count += 1
        else:
            MINENERGY = np.minimum(ENERGY, MINENERGY)

    # Finish off
    return(-100 * MINENERGY)






class SpectrophoreCalculator:
    """
    Class to calculate spectrophores

    Usage:

        SpectrophoreCalculator(parameters)

    Parameters:

        resolution = float_value [default = 3.0]
        accuracy = 1|2|3|4|5|6|9|10|12|15|18|20|30|36|45|60|90|180 [default = 20]
        stereo = none|unique|mirror|all [default = "none"]
        normalization = none|mean|std|all [default = "none"]

    Returns:

        numpy.ndarray
    """

    # #####################################
    # ####### VARIABLES DESCRIPTION #######
    # #####################################
    #
    # #####################################
    # self.PARMS[9]
    # #####################################
    #
    # self.PARMS[0]     Normalisation (0 = none, 1 = mean, 2 = std, 3 = all
    # self.PARMS[1]     Accuracy (integer value specifying the angular stepsize in degrees)
    # self.PARMS[2]     Stereo (0 = none, 1 = unique, 2 = mirror, 3 = all)
    # self.PARMS[3]     BeginProbe (integer, starting from 0)
    # self.PARMS[4]     EndProbe (integer, value itself is not included anymore)
    # self.PARMS[5]     Size of spectrophore (integer, equal to number of properties * number of probes)
    # self.PARMS[6]     Number of probes (integer, equal to EndProbe - BeginProbe)
    # self.PARMS[7]     Number of properties (integer, equal to 4)
    # self.PARMS[8]     Number of box points (integer, equal to 12)
    #
    #
    # #####################################
    # self.COORD[natoms][17]
    # #####################################
    #
    # self.COORD[i][0]      Original x-coordinate of atom i
    # self.COORD[i][1]      Original y-coordinate of atom i
    # self.COORD[i][2]      Original z-coordinate of atom i
    #
    # self.COORD[i][3]      Rotated x-coordinate of atom i
    # self.COORD[i][4]      Rotated y-coordinate of atom i
    # self.COORD[i][5]      Rotated z-coordinate of atom i
    #
    # self.COORD[i][6]      Rotated x-coordinate of atom i - radius of atom i
    # self.COORD[i][7]      Rotated y-coordinate of atom i - radius of atom i
    # self.COORD[i][8]      Rotated z-coordinate of atom i - radius of atom i
    #
    # self.COORD[i][9]      Rotated x-coordinate of atom i + radius of atom i
    # self.COORD[i][10]     Rotated y-coordinate of atom i + radius of atom i
    # self.COORD[i][11]     Rotated z-coordinate of atom i + radius of atom i
    #
    # self.COORD[i][12]     Radius of atom i
    #
    # self.COORD[i][13]     Atomic property 0 of atom i (atomic partial charges)
    # self.COORD[i][14]     Atomic property 1 of atom i (atomic lipophilicities)
    # self.COORD[i][15]     Atomic property 2 of atom i (atomic shape deviations)
    # self.COORD[i][16]     Atomic property 3 of atom i (atomic electrophilicities)
    #
    #
    # #####################################
    # self.SPHORE[size of spectrophore]
    # #####################################
    #
    # self.SPHORE[n]        Spectrophore of size self.PARMS[5]
    #
    #
    # #####################################
    # self.ENERGY[number of probes][number of properties]
    # self.MINENERGY[number of probes][number of properties]
    # #####################################
    #
    # self.ENERGY[i][j]     Energy for probe i and property j
    # self.MINENERGY[i][j]  Minimum energy for probe i and property j
    #
    #
    # #####################################
    # self.ANGLES[number of rotations][3]
    # #####################################
    #
    # self.ANGLES[i][0]     Angle alpha of rotation i
    # self.ANGLES[i][1]     Angle beta of rotation i
    # self.ANGLES[i][2]     Angle gamma of rotation i
    #
    #
    # #####################################
    # self.ROTMAT[3][3]
    # #####################################
    #
    # self.ROTMAT[0][0]     Rotation matrix element 0 x 0
    #
    #
    # #####################################
    # self.PROBES[48][12]
    # #####################################
    #
    # self.PROBES[i][j]     Probe value of the i'th probe (1-48) and the j'th box point (1-12)
    #
    #
    # #####################################
    # self.BOXES[number of box points][number of probes + 3]
    # #####################################
    #
    # self.BOXES[i][j]                   Probe value of the i'th box point (1-12) of the j'th probe (from BeginProbe to EndProbe)
    # self.BOXES[i][numberOProbes + 0]   x-coordinate of the i'th box point (1-12)
    # self.BOXES[i][numberOProbes + 1]   y-coordinate of the i'th box point (1-12)
    # self.BOXES[i][numberOProbes + 2]   z-coordinate of the i'th box point (1-12)




    ####################################################
    def __init__(self, resolution=3.0, accuracy=20, stereo='none', normalization='none'):

        # Initiate PARMS
        self.PARMS = np.array([
            0,      #  0 Normalisation
            20,     #  1 Accuracy
            0,      #  2 Stereo
            0,      #  3 BeginProbe
            12,     #  4 EndProbe
            4*12,   #  5 Size of spectrophore
            12,     #  6 Number of probes
            4,      #  7 Number of properties
            12      #  8 Number of box points
            ])


        # Initiate PROBES
        self.PROBES = np.array([
            #  1 / Dodecapole - non-stereo - probe 1
            [+1, +1, -1, -1, -1, +1, +1, -1, -1, -1, +1, +1],
            #  2 / Dodecapole - non-stereo - probe 2
            [+1, +1, -1, -1, +1, -1, -1, +1, -1, -1, +1, +1],
            #  3 / Dodecapole - non-stereo - probe 3
            [+1, +1, +1, -1, -1, -1, -1, -1, -1, +1, +1, +1],
            #  4 / Dodecapole - non-stereo - probe 4
            [+1, +1, +1, -1, -1, -1, -1, -1, +1, +1, -1, +1],
            #  5 / Dodecapole - non-stereo - probe 5
            [+1, +1, +1, -1, -1, +1, -1, +1, -1, -1, +1, -1],
            #  6 / Dodecapole - non-stereo - probe 6
            [+1, +1, +1, -1, +1, -1, +1, -1, -1, -1, +1, -1],
            #  7 / Dodecapole - non-stereo - probe 7
            [+1, +1, +1, -1, +1, -1, +1, -1, +1, -1, -1, -1],
            #  8 / Dodecapole - non-stereo - probe 8
            [+1, +1, +1, +1, -1, -1, -1, -1, +1, -1, +1, -1],
            #  9 / Dodecapole - non-stereo - probe 9
            [+1, +1, +1, +1, -1, -1, -1, -1, +1, +1, -1, -1],
            # 10 / Dodecapole - non-stereo - probe 10
            [+1, +1, +1, +1, +1, -1, -1, +1, -1, -1, -1, -1],
            # 11 / Dodecapole - non-stereo - probe 11
            [+1, +1, +1, +1, +1, +1, -1, -1, -1, -1, -1, -1],
            # 12 / Dodecapole - non-stereo - probe 12
            [+1, +1, +1, -1, -1, +1, -1, -1, -1, +1, -1, +1],
            # 13 / Dodecapole - mirror-stereo - probe 1
            [+1, +1, -1, -1, -1, -1, +1, +1, -1, +1, +1, -1],
            # 14 / Dodecapole - mirror-stereo - probe 2
            [+1, +1, +1, -1, -1, -1, -1, -1, +1, -1, +1, +1],
            # 15 / Dodecapole - mirror-stereo - probe 3
            [+1, +1, +1, -1, -1, -1, -1, +1, -1, +1, +1, -1],
            # 16 / Dodecapole - mirror-stereo - probe 4
            [+1, +1, +1, -1, -1, -1, +1, -1, -1, +1, -1, +1],
            # 17 / Dodecapole - mirror-stereo - probe 5
            [+1, +1, +1, -1, -1, -1, +1, -1, -1, +1, +1, -1],
            # 18 / Dodecapole - mirror-stereo - probe 6
            [+1, +1, +1, -1, -1, -1, +1, -1, +1, -1, +1, -1],
            # 19 / Dodecapole - mirror-stereo - probe 7
            [+1, +1, +1, -1, -1, -1, +1, -1, +1, +1, -1, -1],
            # 20 / Dodecapole - mirror-stereo - probe 8
            [+1, +1, +1, -1, -1, -1, +1, +1, -1, +1, -1, -1],
            # 21 / Dodecapole - mirror-stereo - probe 9
            [+1, +1, +1, -1, -1, -1, +1, +1, +1, -1, -1, -1],
            # 22 / Dodecapole - mirror-stereo - probe 10
            [+1, +1, +1, -1, -1, +1, +1, -1, -1, -1, -1, +1],
            # 23 / Dodecapole - mirror-stereo - probe 11
            [+1, +1, +1, -1, -1, +1, +1, -1, -1, -1, +1, -1],
            # 24 / Dodecapole - mirror-stereo - probe 12
            [+1, +1, +1, -1, -1, +1, +1, -1, -1, +1, -1, -1],
            # 25 / Dodecapole - mirror-stereo - probe 13
            [+1, +1, +1, -1, -1, +1, +1, +1, -1, -1, -1, -1],
            # 26 / Dodecapole - mirror-stereo - probe 14
            [+1, +1, +1, -1, +1, -1, -1, +1, +1, -1, -1, -1],
            # 27 / Dodecapole - mirror-stereo - probe 15
            [+1, +1, +1, -1, +1, -1, +1, -1, -1, -1, -1, +1],
            # 28 / Dodecapole - mirror-stereo - probe 16
            [+1, +1, +1, -1, +1, -1, +1, +1, -1, -1, -1, -1],
            # 29 / Dodecapole - mirror-stereo - probe 17
            [+1, +1, +1, +1, +1, -1, -1, -1, -1, -1, -1, +1],
            # 30 / Dodecapole - mirror-stereo - probe 18
            [+1, +1, +1, +1, +1, -1, -1, -1, -1, -1, +1, -1],
            # 31 / Dodecapole - unique-stereo - probe 1
            [+1, +1, -1, -1, +1, -1, +1, -1, +1, -1, -1, +1],
            # 32 / Dodecapole - unique-stereo - probe 2
            [+1, +1, +1, -1, -1, -1, -1, -1, +1, +1, +1, -1],
            # 33 / Dodecapole - unique-stereo - probe 3
            [+1, +1, +1, -1, -1, +1, -1, -1, -1, -1, +1, +1],
            # 34 / Dodecapole - unique-stereo - probe 4
            [+1, +1, +1, -1, +1, -1, -1, -1, -1, +1, -1, +1],
            # 35 / Dodecapole - unique-stereo - probe 5
            [+1, +1, +1, -1, +1, -1, -1, -1, -1, -1, +1, +1],
            # 36 / Dodecapole - unique-stereo - probe 6
            [+1, +1, +1, -1, +1, -1, -1, -1, +1, -1, +1, -1],
            # 37 / Dodecapole - unique-stereo - probe 7
            [+1, +1, +1, -1, +1, -1, -1, -1, +1, -1, -1, +1],
            # 38 / Dodecapole - unique-stereo - probe 8
            [+1, +1, +1, -1, +1, +1, -1, -1, -1, -1, -1, +1],
            # 39 / Dodecapole - unique-stereo - probe 9
            [+1, +1, +1, -1, +1, +1, -1, -1, +1, -1, -1, -1],
            # 40 / Dodecapole - unique-stereo - probe 10
            [+1, +1, +1, -1, +1, -1, -1, +1, -1, +1, -1, -1],
            # 41 / Dodecapole - unique-stereo - probe 11
            [+1, +1, +1, -1, +1, -1, -1, +1, -1, -1, +1, -1],
            # 42 / Dodecapole - unique-stereo - probe 12
            [+1, +1, +1, -1, +1, -1, -1, +1, -1, -1, -1, +1],
            # 43 / Dodecapole - unique-stereo - probe 13
            [+1, +1, +1, -1, +1, +1, -1, +1, -1, -1, -1, -1],
            # 44 / Dodecapole - unique-stereo - probe 14
            [+1, +1, +1, -1, +1, +1, -1, -1, -1, -1, +1, -1],
            # 45 / Dodecapole - unique-stereo - probe 15
            [+1, +1, +1, -1, +1, -1, +1, -1, -1, +1, -1, -1],
            # 46 / Dodecapole - unique-stereo - probe 16
            [+1, +1, +1, -1, +1, +1, +1, -1, -1, -1, -1, -1],
            # 47 / Dodecapole - unique-stereo - probe 17
            [+1, +1, +1, +1, +1, -1, -1, -1, +1, -1, -1, -1],
            # 48 / Dodecapole - unique-stereo - probe 18
            [+1, +1, +1, +1, +1, -1, -1, -1, -1, +1, -1, -1]
        ])
        print("Probes initialised: %d number of probes in total" % (len(self.PROBES)))


        # Initiate resolution
        if resolution > 0: self.RESOLUTION = resolution
        else: raise ValueError('Resolution should be larger than 0')


        # Initiate the type of normalization
        if   normalization.lower() == 'none': self.PARMS[0] = 0
        elif normalization.lower() == 'mean': self.PARMS[0] = 1
        elif normalization.lower() == 'std':  self.PARMS[0] = 2
        elif normalization.lower() == 'all':  self.PARMS[0] = 3
        else: raise ValueError(
            'The normalization flag should be "none", "mean", "std" or "all"')


        # Initiate accuracy
        if (180 % int(accuracy)) == 0: self.PARMS[1] = int(accuracy)
        else: raise ValueError('(180 modus accuracy) should be equal to 0')
        self.ANGLES = []
        for a in range(0, 360, self.PARMS[1]):
	        for b in range(0, 360, self.PARMS[1]):
		        for c in range(0, 180, self.PARMS[1]):
			        self.ANGLES.append([math.radians(a), math.radians(b), math.radians(c)])
        self.ANGLES = np.array(self.ANGLES)


        # Initiate stereo
        if   stereo.lower() == 'none':   self.PARMS[2:7] = [0, 0,12,4*12,12]
        elif stereo.lower() == 'unique': self.PARMS[2:7] = [1,12,30,4*18,18]
        elif stereo.lower() == 'mirror': self.PARMS[2:7] = [2,30,48,4*18,18]
        elif stereo.lower() == 'all':    self.PARMS[2:7] = [3,12,48,4*36,36]
        else: raise ValueError('The stereo flag should be "none", "unique", "mirror" or "all"')


        # Setup the boxes
        self.BOXES = np.array(np.zeros(self.PARMS[8] * (self.PARMS[6] + 3))).reshape(self.PARMS[8], self.PARMS[6] + 3)
        print("Only using %d probes" % (self.PARMS[6]))
        for probe in range(self.PARMS[3], self.PARMS[4]):   # from beginProbe to endProbe
            for boxpoint in range(self.PARMS[8]):           # loop over all box points
                self.BOXES[boxpoint][probe] = self.PROBES[probe][boxpoint]




    ####################################################
    def resolution(self, resolution=None):
        if resolution is None: return self.RESOLUTION
        elif resolution > 0: self.RESOLUTION = resolution
        else: raise ValueError('Resolution should be larger than 0')





    ####################################################
    def normalization(self, normalization=None):
        if normalization is None:
            if   self.PARMS[0] == 0: return 'none'
            elif self.PARMS[0] == 1: return 'mean'
            elif self.PARMS[0] == 2: return 'std'
            elif self.PARMS[0] == 3: return 'all'
        else:
            if   normalization.lower() == 'none': self.PARMS[0] = 0
            elif normalization.lower() == 'mean': self.PARMS[0] = 1
            elif normalization.lower() == 'std':  self.PARMS[0] = 2
            elif normalization.lower() == 'all':  self.PARMS[0] = 3
            else: raise ValueError(
            'The normalization flag should be "none", "mean", "std" or "all"')




    ####################################################
    def accuracy(self, accuracy=None):
        if accuracy is None: return self.PARMS[1]
        elif (180 % accuracy) == 0: self.PARMS[1] = int(accuracy)
        else: raise ValueError('(180 modus accuracy) should be equal to 0')
        self.ANGLES = []
        for a in range(0, 360, self.PARMS[1]):
	        for b in range(0, 360, self.PARMS[1]):
		        for c in range(0, 180, self.PARMS[1]):
			        self.ANGLES.append([math.radians(a), math.radians(b), math.radians(c)])
        self.ANGLES = np.array(self.ANGLES)




    ####################################################
    def stereo(self, stereo=None):
        if stereo is None:
            if   self.PARMS[2] == 0: return 'none'
            elif self.PARMS[2] == 1: return 'unique'
            elif self.PARMS[2] == 2: return 'mirror'
            elif self.PARMS[2] == 3: return 'all'
        else:
            if   stereo.lower() == 'none':   self.PARMS[2:7] = [0, 0,12,4*12,12]
            elif stereo.lower() == 'unique': self.PARMS[2:7] = [1,12,30,4*18,18]
            elif stereo.lower() == 'mirror': self.PARMS[2:7] = [2,30,48,4*18,18]
            elif stereo.lower() == 'all':    self.PARMS[2:7] = [3,12,48,4*36,36]
            else: raise ValueError('The stereo flag should be "none", "unique", "mirror" or "all"')

        # Setup the boxes
        self.BOXES = np.array(np.zeros(self.PARMS[8] * (self.PARMS[6] + 3))).reshape(self.PARMS[8], self.PARMS[6] + 3)
        print("Only using %d probes" % (self.PARMS[6]))
        for probe in range(self.PARMS[3], self.PARMS[4]):   # from beginProbe to endProbe
            for boxpoint in range(self.PARMS[8]):           # loop over all box points
                self.BOXES[boxpoint][probe] = self.PROBES[probe][boxpoint]





    ####################################################
    def calculate(self, mol, confID=0):

        # Wrong conformation flag
        wrong3d = False

        # Check number of atoms after adding the hydrogens
        nAtoms = mol.GetNumAtoms()
        if nAtoms < 3: raise ValueError( '>=3 atoms are needed in molecule, only %d given' % (nAtoms))

        # Create the COORD array
        self.COORD = np.array(np.zeros(nAtoms * 17)).reshape(nAtoms, 17)

        # Atomic properties
        # [0]: atomic partial charges -> conformation dependent
        # [1]: atomic lipophilicities
        # [2]: atomic shape deviations -> conformation dependent
        # [3]: atomic electrophilicities -> conformation dependent

        chi = np.zeros(nAtoms)
        eta = np.zeros(nAtoms)
        A = np.zeros((nAtoms + 1, nAtoms + 1))
        B = np.zeros(nAtoms + 1)

        a = 0
        for atom in mol.GetAtoms():
            n = atom.GetAtomicNum()
            if   n ==  1:   # H
                self.COORD[a][12] = +1.20
                eta[a] = +0.65971
                chi[a] = +0.20606
                if atom.GetTotalValence():
                    neighbors = atom.GetNeighbors()
                    self.COORD[a][14] = -0.018
                    for neighbor in neighbors:
                        an = neighbor.GetAtomicNum()
                        if an != 1 and an != 6:
                            self.COORD[a][14] = -0.374
                            break
                else:
                    prop[a][1] = -0.175
            elif n ==  3:   # Li
                self.COORD[a][12] = 1.82
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n ==  5:   # B
                self.COORD[a][12] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.32966
                self.COORD[a][14] = -0.175
            elif n ==  6:   # C
                self.COORD[a][12] = 1.70
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = +0.271
            elif n ==  7:   # N
                self.COORD[a][12] = 1.55
                eta[a] = +0.34519
                chi[a] = +0.49279
                self.COORD[a][14] = -0.137
            elif n ==  8:   # O
                self.COORD[a][12] = 1.52
                eta[a] = +0.54428
                chi[a] = +0.73013
                self.COORD[a][14] = -0.321
            elif n ==  9:   # F
                self.COORD[a][12] = 1.47
                eta[a] = +0.72664
                chi[a] = +0.72052
                self.COORD[a][14] = +0.217
            elif n == 11:   # Na
                self.COORD[a][12] = 2.27
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 12:   # Mg
                self.COORD[a][12] = 1.73
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 14:   # Si
                self.COORD[a][12] = 2.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 15:   # P
                self.COORD[a][12] = 1.80
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 16:   # S
                self.COORD[a][12] = 1.80
                eta[a] = +0.20640
                chi[a] = +0.62020
                self.COORD[a][14] = +0.385
            elif n == 17:   # Cl
                self.COORD[a][12] = 1.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = +0.632
            elif n == 19:   # K
                self.COORD[a][12] = 2.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 20:   # Ca
                self.COORD[a][12] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 26:   # Fe
                self.COORD[a][12] = 1.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 29:   # Cu
                self.COORD[a][12] = 1.40
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 30:   # Zn
                self.COORD[a][12] = 1.39
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.COORD[a][14] = -0.175
            elif n == 35:   # Br
                self.COORD[a][12] = 1.85
                eta[a] = +0.54554
                chi[a] = +0.70052
                self.COORD[a][14] = +0.815
            elif n == 53:   # I
                self.COORD[a][12] = 1.98
                eta[a] = +0.30664
                chi[a] = +0.68052
                self.COORD[a][14] = +0.198
            else:
                self.COORD[a][12] = 1.50
                eta[a] = +0.65971
                chi[a] = +0.20606
                self.COORD[a][14] = -0.175
            a += 1

        # Conformers
        if mol.GetNumConformers() < confID + 1:
            raise ValueError(
                'At least %d conformation(s) should be present, %d found' %
                (confID + 1, mol.GetNumConformers()))
        conf = mol.GetConformer(confID)

        # Coordinates
        for r in range(nAtoms):
            c = conf.GetAtomPosition(r)
            for i in range(3): self.COORD[r][i] = c[i]
            A[r][r] = 2 * eta[r]

        # Complete A matrix
        for r in range(nAtoms):
            for i in range(r + 1, nAtoms):
                d =  (self.COORD[r][0] - self.COORD[i][0])**2
                d += (self.COORD[r][1] - self.COORD[i][1])**2
                d += (self.COORD[r][2] - self.COORD[i][2])**2

                if d == 0: return(np.zeros(self.PARMS[6] * self.PARMS[7]))
                d = np.sqrt(d)
                A[r][i] = 0.529176 / d    # Angstrom to au
                A[i][r] = A[r][i]

        # Property [0]: partial atomic charges
        for i in range(nAtoms):
            A[i][nAtoms] = -1
            A[nAtoms][i] = +1
            B[i] = -chi[i]
        A[nAtoms][nAtoms] = 0
        B[nAtoms] = Chem.GetFormalCharge(mol)
        X = scipy.linalg.solve(A, B)
        chi2 = X[nAtoms] * X[nAtoms]
        for a in range(nAtoms): self.COORD[a][13] = X[a]

        # Property [2]: atomic shape deviations
        cog = np.mean(self.COORD[:,0:3],0)
        d = np.zeros(nAtoms)
        for a in range(nAtoms):
            d[a] = math.sqrt((self.COORD[a][0] - cog[0])**2 + (self.COORD[a][1] - cog[1])**2 + (self.COORD[a][2] - cog[2])**2)
        avg_d = np.average(d)
        std_d = np.std(d)
        for a in range(nAtoms): self.COORD[a][15] = avg_d + ((d[a] - avg_d) / std_d)

        # Property [3]: atomic electrophilicities
        B = np.ones(nAtoms + 1)
        B[nAtoms] = 0
        for i in range(nAtoms):
            A[i][nAtoms] = 0
            A[nAtoms][i] = +1
        A[nAtoms][nAtoms] = -1
        X = scipy.linalg.solve(A, B)
        for a in range(nAtoms): self.COORD[a][16] = X[a] * chi2

        # Orient molecule to its center of gravity and orient in standard way
        # 1) Center molecule around its center of gravity
        self.COORD[:,0:3] -= cog

        # 2) Determine atom that is furthest away from origin
        d = self.COORD[:,0:3]**2
        d = np.sqrt(d.sum(axis=1))
        maxAtom = np.argmax(d)

        # 3) Rotate all atoms along z-axis
        angle = -np.arctan2(self.COORD[maxAtom][1], self.COORD[maxAtom][0])
        c = np.cos(angle)
        s = np.sin(angle)
        for i in range(nAtoms):
            x = c * self.COORD[i][0] - s * self.COORD[i][1]
            y = s * self.COORD[i][0] + c * self.COORD[i][1]
            self.COORD[i][0] = x
            self.COORD[i][1] = y

        # 4) Rotate all atoms along y-axis to place the maxAtom on z
        angle = -np.arctan2(self.COORD[maxAtom][0], self.COORD[maxAtom][2])
        c = np.cos(angle)
        s = np.sin(angle)
        for i in range(nAtoms):
            x = c * self.COORD[i][0] + s * self.COORD[i][2]
            z = c * self.COORD[i][2] - s * self.COORD[i][0]
            self.COORD[i][0] = x
            self.COORD[i][2] = z

        # 5) Center molecule again around its COG
        cog = np.mean(self.COORD[:,0:3],0)
        self.COORD[:,0:3] -= cog

        # Rotate
        self.ENERGY = np.zeros(self.PARMS[6] * self.PARMS[7])
        self.MINENERGY = np.zeros(self.PARMS[6] * self.PARMS[7])
        self.ROTMAT = np.ndarray(shape=(3,3))
        sphore = rotate(self.PARMS, self.COORD, self.ANGLES, self.ROTMAT, self.ENERGY, self.MINENERGY, self.BOXES, self.RESOLUTION)

        # Normalise
        if self.PARMS[0] == 0: return(sphore)
        else:
            t = sphore.reshape(self.PARMS[7], self.PARMS[6])
            m = np.mean(t,1)
            s = np.std(t,1)

            if self.PARMS[0] == 1:
                for r in range(self.PARMS[7]): t[r,:] = t[r,:] - m[r]
            elif self.PARMS[0] == 2:
                for r in range(self.PARMS[7]): t[r,:] = t[r,:] / s[r]
            elif self.PARMS[0] == 3:
                for r in range(self.PARMS[7]): t[r,:] = (t[r,:] - m[r]) / s[r]
            return(t.flatten())



# #############################################################
# Main
# #############################################################

import progressbar
import argparse


# Function to calculate spectrophores
def __calculateSpectrophore(molecule, calculator, label):
	spec = calculator.calculate(molecule)
	if not np.all(spec): return None
	specString = np.array2string(spec, max_line_width=1000, suppress_small=True, formatter={'float':lambda x: "%.5f" % x})
	specString = label + " " + specString[1:-1]
	return specString
 

# Function to check the value of the resolution command-line argument
def __check_resolution(value):
	v = float(value)
	if v <= 0: raise argparse.ArgumentTypeError("should be larger than 0")
	return v


# Function to parse the command-line arguments
def __processCommandline():
    parser = argparse.ArgumentParser(description = "Calculate spectrophores",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--norm', 
                        dest = 'norm', 
                        help = 'normalization setting', 
                        default = 'none', 
                        choices = ['none','mean','all','std'])
    parser.add_argument('-s', '--stereo', 
                        dest = 'stereo', 
                        help = 'stereo setting', 
                        default = 'none', 
                        choices = ['none','unique','mirror','all'])
    parser.add_argument('-a', '--accuracy', 
                        dest = 'accuracy', 
                        help = 'accuracy setting', 
                        type = int, 
                        default = 20, 
                        choices = [1, 2, 3, 4, 5, 6, 9, 10, 12, 15, 18, 20, 30, 36, 45, 60, 90, 180])
    parser.add_argument('-r', '--resolution', 
                        dest = 'resolution', 
                        help = 'resolution setting (>0)', 
                        default = 3, 
                        type = __check_resolution)
    parser.add_argument('-p', '--np',
                        dest = 'max_workers',
                        help = 'number of processors to use; -1 is all processors',
                        default = -1,
                        type = int)
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i', '--in', 
                               dest = 'infile', 
                               help = 'input sdf file', 
                               required = True)
    requiredNamed.add_argument('-o', '--out', 
                               dest = 'outfile', 
                               help = 'output spectrophore file', 
                               required = True)
    return parser.parse_args()



# Main
if __name__ == "__main__":
    args = __processCommandline()

    calculator = spectrophore.SpectrophoreCalculator(normalization = args.norm, 
                                        stereo = args.stereo, 
                                        accuracy = args.accuracy, 
                                        resolution = args.resolution)
    print("Normalization: ", args.norm)
    print("Stereo:        ", args.stereo)
    print("Accuracy:      ", args.accuracy)
    print("Resolution:    ", args.resolution)
    print("Input file:    ", args.infile)
    print("Output file:   ", args.outfile)
    print("Processors:    ", args.max_workers)
    supplier = Chem.SDMolSupplier(args.infile, removeHs=False)
    of = open(args.outfile, 'w')
    mw = None
    if args.max_workers > 0: mw = args.max_workers
    with futures.ProcessPoolExecutor(max_workers = mw) as executor:
	    jobs = []
	    for mol in supplier:
	    	if mol:
		    	label = mol.GetProp("_Name")
		    	job = executor.submit(__calculateSpectrophore, mol, calculator, label)
		    	jobs.append(job)
	
	    widgets = ["Generating spectrophores; ", progressbar.Percentage(), " ", progressbar.ETA(), " ", progressbar.Bar()]
	    pbar = progressbar.ProgressBar(widgets = widgets, maxval = len(jobs))
	    for job in pbar(futures.as_completed(jobs)):
		    spec = job.result()
		    if spec is not None: of.write("%s\n" % (spec))

    of.close()

