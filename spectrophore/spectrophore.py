#!/usr/bin/env python

import math
import numpy as np

import scipy.linalg

from numba import njit
from rdkit import Chem


@njit(fastmath=True)
def rotate(xx, xy, xz, yx, yy, yz, zx, zy, zz, c):
    """
    Rotate a tuple of input coordinates using rotation matrix coefficients
    """
    xr = c[0] * xx + c[1] * xy + c[2] * xz
    yr = c[0] * yx + c[1] * yy + c[2] * yz
    zr = c[0] * zx + c[1] * zy + c[2] * zz
    return (xr, yr, zr)


@njit(parallel=False, fastmath=True)
def calculateSpectrophore(
    COORD_ORI,
    COORD_ROT,
    STEPSIZE,
    BOX,
    PROBES,
    RADIUS,
    PROP,
    RESOLUTION,
):
    """
    Calculate a spectrophore fingerprint
    """
    nFramesA = int(360 / STEPSIZE)
    nFramesB = int(360 / STEPSIZE)
    nFramesG = int(180 / STEPSIZE) + 1

    # Initialise a numpy array with the correct dimensions to save the sampled values
    OUTPUT_ARRAY = np.zeros(
        (nFramesA * nFramesB * nFramesG, len(PROBES) * 4), dtype=np.float32
    )

    # Initialise outer limits of molecule, taking into account atom radius and resolution
    px, py, pz = COORD_ROT[0] + RADIUS[0] + RESOLUTION
    mx, my, mz = COORD_ROT[0] - RADIUS[0] - RESOLUTION

    for atom in range(1, len(COORD_ROT)):
        x, y, z = COORD_ROT[atom] + RADIUS[atom] + RESOLUTION

        px = max(x, px)
        py = max(y, py)
        pz = max(z, pz)

        x, y, z = COORD_ROT[atom] - RADIUS[atom] - RESOLUTION

        mx = min(x, mx)
        my = min(y, my)
        mz = min(z, mz)

    hx = (mx + px) / 2.0
    hy = (my + py) / 2.0
    hz = (mz + pz) / 2.0

    BOX[0] = hx, my, pz
    BOX[1] = px, hy, pz
    BOX[2] = hx, py, pz
    BOX[3] = mx, hy, pz
    BOX[4] = mx, my, hz
    BOX[5] = px, my, hz
    BOX[6] = mx, py, hz
    BOX[7] = px, py, hz
    BOX[8] = px, hy, mz
    BOX[9] = hx, my, mz
    BOX[10] = mx, hy, mz
    BOX[11] = hx, py, mz

    frameIdx = 0

    # Loop over all angles
    for ia in range(0, nFramesA):
        for ib in range(0, nFramesB):
            for ig in range(0, nFramesG):
                # Rotate
                alpha = math.radians(float(ia * STEPSIZE))
                beta = math.radians(float(ib * STEPSIZE))
                gamma = math.radians(float(ig * STEPSIZE))

                ca = math.cos(alpha)
                cb = math.cos(beta)
                cc = math.cos(gamma)
                sa = math.sin(alpha)
                sb = math.sin(beta)
                sc = math.sin(gamma)

                casb = ca * sb
                sasb = sa * sb
                cacc = ca * cc

                xx = ca * cb
                xy = casb * sc - sa * cc
                xz = casb * cc + sa * sc
                yx = sa * cb
                yy = sasb * sc + cacc
                yz = sasb * cc - ca * sc
                zx = -sb
                zy = cb * sc
                zz = cacc

                for a in range(COORD_ORI.shape[0]):
                    COORD_ROT[a] = rotate(
                        xx, xy, xz, yx, yy, yz, zx, zy, zz, COORD_ORI[a]
                    )

                # Update outer limits of molecule, taking into account atom radius and resolution
                px, py, pz = COORD_ROT[0] + RADIUS[0] + RESOLUTION
                mx, my, mz = COORD_ROT[0] - RADIUS[0] - RESOLUTION

                for atom in range(1, len(COORD_ROT)):
                    x, y, z = COORD_ROT[atom] + RADIUS[atom] + RESOLUTION

                    px = max(x, px)
                    py = max(y, py)
                    pz = max(z, pz)

                    x, y, z = COORD_ROT[atom] - RADIUS[atom] - RESOLUTION

                    mx = min(x, mx)
                    my = min(y, my)
                    mz = min(z, mz)

                hx = (mx + px) / 2.0
                hy = (my + py) / 2.0
                hz = (mz + pz) / 2.0

                BOX[0] = hx, my, pz
                BOX[1] = px, hy, pz
                BOX[2] = hx, py, pz
                BOX[3] = mx, hy, pz
                BOX[4] = mx, my, hz
                BOX[5] = px, my, hz
                BOX[6] = mx, py, hz
                BOX[7] = px, py, hz
                BOX[8] = px, hy, mz
                BOX[9] = hx, my, mz
                BOX[10] = mx, hy, mz
                BOX[11] = hx, py, mz

                # Initialise array to store spectrophore data for this frame
                FRAME_ARRAY = np.zeros(OUTPUT_ARRAY.shape[-1], dtype=np.float32)

                # Calculate energies for this frame
                # Loop over each boxpoint (12 points)
                for boxPoint in range(12):
                    # Loop over each atom
                    for atom in range(len(COORD_ROT)):
                        # Distance between boxpoint and atom
                        d = math.sqrt(np.sum((BOX[boxPoint] - COORD_ROT[atom]) ** 2))

                        # Loop over each probe
                        for probe in range(len(PROBES)):
                            # Loop over each property (4 properties)
                            for prop in range(4):
                                index = (len(PROBES) * prop) + probe
                                FRAME_ARRAY[index] += (
                                    PROP[atom][prop] * PROBES[probe][boxPoint]
                                ) / d

                OUTPUT_ARRAY[frameIdx] = FRAME_ARRAY
                frameIdx += 1

    return -100 * OUTPUT_ARRAY


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
        mode = classic|full

    Returns:

        numpy.ndarray
    """

    # #####################################
    # ####### VARIABLES DESCRIPTION #######
    # #####################################
    #
    # #####################################
    # self.PARAMS[9]
    # #####################################
    #
    # self.PARAMS[0]     normalization (0 = none, 1 = mean, 2 = std, 3 = all)
    # self.PARAMS[1]     Mode (0 = classic, 1 = full)
    # self.PARAMS[2]     Accuracy (integer value specifying the angular stepsize in degrees)
    # self.PARAMS[3]     Stereo (0 = none, 1 = unique, 2 = mirror, 3 = all)
    # self.PARAMS[4]     BeginProbe (integer, starting from 0)
    # self.PARAMS[5]     EndProbe (integer, value itself is not included anymore)
    # self.PARAMS[6]     Size of spectrophore (integer, equal to number of properties * number of probes)
    # self.PARAMS[7]     Number of probes (integer, equal to EndProbe - BeginProbe)
    #
    #
    # #####################################
    # self.COORD_ORI[natoms][3]
    # self.COORD_ROT[natoms][3]
    # #####################################
    #
    # self.COORD_ORI[i][0]      Original x-coordinate of atom i
    # self.COORD_ORI[i][1]      Original y-coordinate of atom i
    # self.COORD_ORI[i][2]      Original z-coordinate of atom i
    #
    # self.COORD_ROT[i][1]      Rotated x-coordinate of atom i
    # self.COORD_ROT[i][1]      Rotated y-coordinate of atom i
    # self.COORD_ROT[i][2]      Rotated z-coordinate of atom i
    #
    #
    # #####################################
    # self.RADIUS[natoms]
    # #####################################
    #
    # self.RADIUS[i]        Radius of atom i
    #
    #
    # #####################################
    # self.PROP[natoms][4]
    # #####################################
    #
    # self.PROP[i][0]      Atomic property 0 of atom i (atomic partial charges)
    # self.PROP[i][1]      Atomic property 1 of atom i (atomic lipophilicities)
    # self.PROP[i][2]      Atomic property 2 of atom i (atomic shape deviations)
    # self.PROP[i][3]      Atomic property 3 of atom i (atomic electrophilicities)
    #
    #
    # #####################################
    # self.PROBES[48][12]
    # #####################################
    #
    # self.PROBES[i][j]     Probe value of the i'th probe and the j'th box point (1-12)
    #
    #
    # #####################################
    # self.BOX[number of box points][3]
    # #####################################
    #
    # self.BOX[i][0]      x-coordinate of the i'th box point (1-12)
    # self.BOX[i][1]      y-coordinate of the i'th box point (1-12)
    # self.BOX[i][2]      z-coordinate of the i'th box point (1-12)
    ####################################################

    def __init__(
        self,
        resolution=3.0,
        accuracy=20,
        stereo="none",
        normalization="all",
        mode="classic",
    ):
        # Initiate PARAMS
        self.PARAMS = np.array(
            [
                3,  #  0 normalization
                0,  #  1 Mode
                20,  #  2 Accuracy
                0,  #  3 Stereo
                0,  #  4 BeginProbe
                12,  #  5 EndProbe
                4 * 12,  #  6 Size of spectrophore
                12,  #  7 Number of probes
            ]
        )

        # Initiate PROBES
        self.PROBES_TEMPLATE = np.array(
            [
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
                [+1, +1, +1, +1, +1, -1, -1, -1, -1, +1, -1, -1],
            ]
        )
        print(
            f"Probes initialised: {len(self.PROBES_TEMPLATE)} number of probes in total"
        )

        # Initiate resolution
        if resolution > 0:
            self.RESOLUTION = resolution
        else:
            raise ValueError("Resolution should be larger than 0")

        # Initiate the type of normalization
        if normalization.lower() == "none":
            self.PARAMS[0] = 0
        elif normalization.lower() == "mean":
            self.PARAMS[0] = 1
        elif normalization.lower() == "std":
            self.PARAMS[0] = 2
        elif normalization.lower() == "all":
            self.PARAMS[0] = 3
        else:
            raise ValueError(
                'The normalization flag should be "none", "mean", "std" or "all"'
            )

        # Initiate the calculation mode
        if mode.lower() == "classic":
            self.PARAMS[1] = 0
        elif mode.lower() == "full":
            self.PARAMS[1] = 1
        else:
            raise ValueError('The mode flag should be "classic" or "full"')

        # Initiate accuracy
        if (180 % int(accuracy)) == 0:
            self.PARAMS[2] = int(accuracy)
        else:
            raise ValueError("(180 modus accuracy) should be equal to 0")

        # Initiate stereo
        if stereo.lower() == "none":
            self.PARAMS[3:8] = [0, 0, 12, 4 * 12, 12]
        elif stereo.lower() == "unique":
            self.PARAMS[3:8] = [1, 12, 30, 4 * 18, 18]
        elif stereo.lower() == "mirror":
            self.PARAMS[3:8] = [2, 30, 48, 4 * 18, 18]
        elif stereo.lower() == "all":
            self.PARAMS[3:8] = [3, 12, 48, 4 * 36, 36]
        else:
            raise ValueError(
                'The stereo flag should be "none", "unique", "mirror" or "all"'
            )
        self.PROBES = self.PROBES_TEMPLATE[self.PARAMS[4] : self.PARAMS[5]]
        print(f"{self.PARAMS[7]} probes are used due to the imposed stereo flag")

        # Setup the box
        self.BOX = np.zeros(12 * 3).reshape(12, 3)

    ####################################################
    def resolution(self, resolution=None):
        if resolution is None:
            return self.RESOLUTION
        elif resolution > 0:
            self.RESOLUTION = resolution
        else:
            raise ValueError("Resolution should be larger than 0")

    ####################################################
    def normalization(self, normalization=None):
        if normalization is None:
            if self.PARAMS[0] == 0:
                return "none"
            elif self.PARAMS[0] == 1:
                return "mean"
            elif self.PARAMS[0] == 2:
                return "std"
            elif self.PARAMS[0] == 3:
                return "all"
        else:
            if normalization.lower() == "none":
                self.PARAMS[0] = 0
            elif normalization.lower() == "mean":
                self.PARAMS[0] = 1
            elif normalization.lower() == "std":
                self.PARAMS[0] = 2
            elif normalization.lower() == "all":
                self.PARAMS[0] = 3
            else:
                raise ValueError(
                    'The normalization flag should be "none", "mean", "std" or "all"'
                )

    ####################################################
    def mode(self, mode=None):
        if mode is None:
            if self.PARAMS[1] == 0:
                return "classic"
            elif self.PARAMS[1] == 1:
                return "full"
        else:
            if mode.lower() == "classic":
                self.PARAMS[1] = 0
            elif mode.lower() == "full":
                self.PARAMS[1] = 1
            else:
                raise ValueError('The mode flag should be "classic" or "full"')

    ####################################################
    def accuracy(self, accuracy=None):
        if accuracy is None:
            return self.PARAMS[2]
        elif (180 % accuracy) == 0:
            self.PARAMS[2] = int(accuracy)
        else:
            raise ValueError("(180 modus accuracy) should be equal to 0")

    ####################################################
    def stereo(self, stereo=None):
        if stereo is None:
            if self.PARAMS[3] == 0:
                return "none"
            elif self.PARAMS[3] == 1:
                return "unique"
            elif self.PARAMS[3] == 2:
                return "mirror"
            elif self.PARAMS[3] == 3:
                return "all"
        else:
            if stereo.lower() == "none":
                self.PARAMS[3:8] = [0, 0, 12, 4 * 12, 12]
            elif stereo.lower() == "unique":
                self.PARAMS[3:8] = [1, 12, 30, 4 * 18, 18]
            elif stereo.lower() == "mirror":
                self.PARAMS[3:8] = [2, 30, 48, 4 * 18, 18]
            elif stereo.lower() == "all":
                self.PARAMS[3:8] = [3, 12, 48, 4 * 36, 36]
            else:
                raise ValueError(
                    'The stereo flag should be "none", "unique", "mirror" or "all"'
                )
            self.PROBES = self.PROBES_TEMPLATE[self.PARAMS[4] : self.PARAMS[5]]
            print(f"{self.PARAMS[7]} probes are used due to the imposed stereo flag")

    ####################################################
    def calculate(self, mol, confID=0):
        # Check number of atoms after adding the hydrogens
        nAtoms = mol.GetNumAtoms()
        if nAtoms < 3:
            raise ValueError(f">=3 atoms are needed in molecule, only {nAtoms} given")

        # Check if conformer exists and is 3D
        if confID not in [x.GetId() for x in mol.GetConformers()]:
            raise ValueError(f"Conformer ID {confID} is not a valid")
        elif not mol.GetConformer(confID).Is3D():
            raise ValueError("The input molecule doesn't have a valid 3D conformation")

        conf = mol.GetConformer(confID)

        # Create the PROP and RADIUS array
        self.PROP = np.zeros(nAtoms * 4).reshape(nAtoms, 4)
        self.RADIUS = np.zeros(nAtoms)

        # Create the COORD arrays
        self.COORD_ORI = np.zeros(nAtoms * 3, dtype=np.float32).reshape(nAtoms, 3)
        self.COORD_ROT = np.zeros(nAtoms * 3, dtype=np.float32).reshape(nAtoms, 3)

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
            if n == 1:  # H
                self.RADIUS[a] = +1.20
                eta[a] = +0.65971
                chi[a] = +0.20606
                if atom.GetTotalValence():
                    neighbors = atom.GetNeighbors()
                    self.PROP[a][1] = -0.018
                    for neighbor in neighbors:
                        an = neighbor.GetAtomicNum()
                        if an != 1 and an != 6:
                            self.PROP[a][1] = -0.374
                            break
                else:
                    self.PROP[a][1] = -0.175
            elif n == 3:  # Li
                self.RADIUS[a] = 1.82
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 5:  # B
                self.RADIUS[a] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.32966
                self.PROP[a][1] = -0.175
            elif n == 6:  # C
                self.RADIUS[a] = 1.70
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = +0.271
            elif n == 7:  # N
                self.RADIUS[a] = 1.55
                eta[a] = +0.34519
                chi[a] = +0.49279
                self.PROP[a][1] = -0.137
            elif n == 8:  # O
                self.RADIUS[a] = 1.52
                eta[a] = +0.54428
                chi[a] = +0.73013
                self.PROP[a][1] = -0.321
            elif n == 9:  # F
                self.RADIUS[a] = 1.47
                eta[a] = +0.72664
                chi[a] = +0.72052
                self.PROP[a][1] = +0.217
            elif n == 11:  # Na
                self.RADIUS[a] = 2.27
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 12:  # Mg
                self.RADIUS[a] = 1.73
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 14:  # Si
                self.RADIUS[a] = 2.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 15:  # P
                self.RADIUS[a] = 1.80
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 16:  # S
                self.RADIUS[a] = 1.80
                eta[a] = +0.20640
                chi[a] = +0.62020
                self.PROP[a][1] = +0.385
            elif n == 17:  # Cl
                self.RADIUS[a] = 1.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = +0.632
            elif n == 19:  # K
                self.RADIUS[a] = 2.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 20:  # Ca
                self.RADIUS[a] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 26:  # Fe
                self.RADIUS[a] = 1.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 29:  # Cu
                self.RADIUS[a] = 1.40
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 30:  # Zn
                self.RADIUS[a] = 1.39
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 35:  # Br
                self.RADIUS[a] = 1.85
                eta[a] = +0.54554
                chi[a] = +0.70052
                self.PROP[a][1] = +0.815
            elif n == 53:  # I
                self.RADIUS[a] = 1.98
                eta[a] = +0.30664
                chi[a] = +0.68052
                self.PROP[a][1] = +0.198
            else:
                self.RADIUS[a] = 1.50
                eta[a] = +0.65971
                chi[a] = +0.20606
                self.PROP[a][1] = -0.175
            a += 1

        # Coordinates
        for r in range(nAtoms):
            c = conf.GetAtomPosition(r)
            for i in range(3):
                self.COORD_ORI[r][i] = c[i]
            A[r][r] = 2 * eta[r]

        # Complete A matrix
        for r in range(nAtoms):
            for i in range(r + 1, nAtoms):
                d = np.sqrt(np.sum((self.COORD_ORI[r] - self.COORD_ORI[i]) ** 2))
                if d == 0:
                    return np.zeros(self.PARAMS[7] * 4)
                A[r][i] = 0.529176 / d  # Angstrom to au
                A[i][r] = A[r][i]

        # Property [0]: partial atomic charges
        A[:-1, nAtoms] = -1
        A[nAtoms, :-1] = +1
        A[nAtoms][nAtoms] = 0
        B[:-1] = -chi
        B[nAtoms] = Chem.GetFormalCharge(mol)
        X = scipy.linalg.solve(A, B)
        chi2 = X[nAtoms] * X[nAtoms]
        self.PROP[:, 0] = X[:-1]

        # Property [2]: atomic shape deviations
        cog = np.mean(self.COORD_ORI, 0)
        d = np.sqrt(np.sum((self.COORD_ORI - cog) ** 2, 1))
        avg_d = np.average(d)
        std_d = np.std(d)
        self.PROP[:, 2] = avg_d + ((d - avg_d) / std_d)

        # Property [3]: atomic electrophilicities
        B = np.ones(nAtoms + 1)
        B[nAtoms] = 0
        A[:-1, nAtoms] = 0
        A[nAtoms, :-1] = 1
        A[nAtoms][nAtoms] = -1
        X = scipy.linalg.solve(A, B)
        for a in range(nAtoms):
            self.PROP[a][3] = X[a] * chi2

        # Orient molecule to its center of gravity and orient in standard way
        # 1) Center molecule around its center of gravity
        self.COORD_ORI -= cog

        # 2) Determine atom that is furthest away from origin
        d = self.COORD_ORI**2
        d = np.sqrt(d.sum(axis=1))
        maxAtom = np.argmax(d)

        # 3) Rotate all atoms along z-axis
        angle = -np.arctan2(self.COORD_ORI[maxAtom][1], self.COORD_ORI[maxAtom][0])
        c = np.cos(angle)
        s = np.sin(angle)
        for i in range(nAtoms):
            x = c * self.COORD_ORI[i][0] - s * self.COORD_ORI[i][1]
            y = s * self.COORD_ORI[i][0] + c * self.COORD_ORI[i][1]
            self.COORD_ORI[i][0] = x
            self.COORD_ORI[i][1] = y

        # 4) Rotate all atoms along y-axis to place the maxAtom on z
        angle = -np.arctan2(self.COORD_ORI[maxAtom][0], self.COORD_ORI[maxAtom][2])
        c = np.cos(angle)
        s = np.sin(angle)
        for i in range(nAtoms):
            x = c * self.COORD_ORI[i][0] + s * self.COORD_ORI[i][2]
            z = c * self.COORD_ORI[i][2] - s * self.COORD_ORI[i][0]
            self.COORD_ORI[i][0] = x
            self.COORD_ORI[i][2] = z

        # 5) Center molecule again around its COG
        cog = np.mean(self.COORD_ORI, 0)
        self.COORD_ORI -= cog

        sphore = calculateSpectrophore(
            self.COORD_ORI,
            self.COORD_ROT,
            self.PARAMS[2],
            self.BOX,
            self.PROBES,
            self.RADIUS,
            self.PROP,
            self.RESOLUTION,
        )

        # Post-process to the correct mode
        if self.PARAMS[1] == 0:
            # NOTE: we use np.max here instead of np.min, since we already multiplied the spectrophore with -100
            sphore = np.max(sphore, axis=0)

        # Normalize
        if not self.PARAMS[0] == 0:
            t = sphore.reshape(-1, 4, self.PARAMS[7])
            m = np.mean(t, axis=-1)
            s = np.std(t, axis=-1)

            if self.PARAMS[0] == 1:
                for r in range(4):
                    t[:, r, :] = t[:, r, :] - m[:, r, np.newaxis]
            elif self.PARAMS[0] == 2:
                for r in range(4):
                    t[:, r, :] = t[:, r, :] / s[:, r, np.newaxis]
            elif self.PARAMS[0] == 3:
                for r in range(4):
                    t[:, r, :] = (t[:, r, :] - m[:, r, np.newaxis]) / s[
                        :, r, np.newaxis
                    ]

            # reshape back
            sphore = t.reshape(-1, self.PARAMS[6]).squeeze()

        return sphore

    def calculate_string(self, mol, confID=0, sep=" "):
        """
        Calculate the spectrophore for a given molecule in String format
        """
        if mol.HasProp("_Name"):
            label = mol.GetProp("_Name")
            if label == "":
                label = Chem.MolToSmiles(Chem.RemoveHs(mol))
        else:
            label = Chem.MolToSmiles(Chem.RemoveHs(mol))

        try:
            spec = self.calculate(mol, confID=confID)

            if not np.all(spec):
                print("Something went wrong with", label)
                return None

            specString = np.array2string(
                spec.flatten(),
                max_line_width=1000,
                suppress_small=True,
                formatter={"float": lambda x: f"{x:.5f}"},
                separator=sep,
            )
            specString = label + " " + specString[1:-1]
            return specString

        except Exception as error:
            print("Error molecule ", label, "  ", type(error).__name__, ": ", error)
            return None
