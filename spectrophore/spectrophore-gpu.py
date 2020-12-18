#!/usr/bin/env python
# coding: utf-8

# In[1]:


from numba import cuda, float32
import numpy as np
import math
import scipy
import scipy.linalg
from rdkit import Chem
from rdkit.Chem import AllChem


# In[2]:


print(cuda.gpus)


# In[18]:


@cuda.jit(device=True)
def rotationMatrix(rm,c):
    xr = c[0]*rm[0][0] + c[1]*rm[0][1] + c[2]*rm[0][2]
    yr = c[0]*rm[1][0] + c[1]*rm[1][1] + c[2]*rm[1][2]
    zr = c[0]*rm[2][0] + c[1]*rm[2][1] + c[2]*rm[2][2]
    return(xr,yr,zr)


# In[19]:


@cuda.jit()
def calculateEnergy(STEPSIZE, CO, RA, RES, PROBES, PROP, ENERGY):

    nFramesA = int(360 / STEPSIZE)
    nFramesB = int(360 / STEPSIZE)
    nFramesG = int(180 / STEPSIZE)

    probe, prop = cuda.grid(2)
    if probe > PROBES.shape[0]: return
    if prop > PROP.shape[0]: return
    
    # Local arrays
    BOX = cuda.local.array((12,3), dtype=float32)
    RM = cuda.local.array((3,3), dtype=float32)

    # First run to initiate ENERGY array
    xr,yr,zr = CO[0]
    px = xr + float(RA[0]) + float(RES)
    py = yr + float(RA[0]) + float(RES)
    pz = zr + float(RA[0]) + float(RES)
    mx = xr - RA[0] - RES
    my = yr - RA[0] - RES
    mz = zr - RA[0] - RES
    for a in range(1, CO.shape[0]):
        xr,yr,zr = rotationMatrix(RM,CO[a])
        x = xr + float(RA[a]) + float(RES)
        y = yr + float(RA[a]) + float(RES)
        z = zr + float(RA[a]) + float(RES)
        if x > px: px = x
        if y > py: py = y
        if z > pz: pz = z
        x = xr - float(RA[a]) - float(RES)
        y = yr - float(RA[a]) - float(RES)
        z = zr - float(RA[a]) - float(RES)
        if x < mx: mx = x
        if y < my: my = y
        if z < mz: mz = z
    hx = (mx + px) / 2.0
    hy = (my + py) / 2.0
    hz = (mz + pz) / 2.0
    BOX[0] =  hx,my,pz
    BOX[1] =  px,hy,pz
    BOX[2] =  hx,py,pz
    BOX[3] =  mx,hy,pz
    BOX[4] =  mx,my,hz
    BOX[5] =  px,my,hz
    BOX[6] =  mx,py,hz
    BOX[7] =  px,py,hz
    BOX[8] =  px,hy,mz
    BOX[9] =  hx,my,mz
    BOX[10] = mx,hy,mz
    BOX[11] = hx,py,mz

    # Calculate energies
    # Loop over each boxpoint (12 points)
    v = 0.0
    for bp in range(len(BOX)):

        # Loop over each atom
        for a in range(len(CO)):
            
            # Rotate atom
            xr,yr,zr = rotationMatrix(RM,CO[a])

            # Distance between boxpoint and atom
            d = ((BOX[bp][0]-xr)**2 + (BOX[bp][1]-yr)**2 + (BOX[bp][2]-zr)**2)**0.5
                        
            # Energy
            v += PROP[a][prop] * PROBES[probe][bp] / d
                        
    # Store in array
    index = PROBES.shape[0] * prop + probe
    ENERGY[index] = v
    
    # Loop over all angles
    for alpha in range(0, nFramesA, STEPSIZE):
        for beta in range(0, nFramesB, STEPSIZE):
            for gamma in range(0, nFramesG, STEPSIZE):
                
                # Rotation matrix
                ca = math.cos(math.radians(float(alpha)))
                cb = math.cos(math.radians(float(beta)))
                cc = math.cos(math.radians(float(gamma)))
                sa = math.sin(math.radians(float(alpha)))
                sb = math.sin(math.radians(float(beta)))
                sc = math.sin(math.radians(float(gamma)))
                casb = ca*sb
                sasb = sa*sb
                cacc = ca*cc
                RM[0][0] = ca*cb
                RM[0][1] = casb*sc - sa*cc
                RM[0][2] = casb*cc + sa*sc
                RM[1][0] = sa*cb
                RM[1][1] = sasb*sc + cacc
                RM[1][2] = sasb*cc - ca*sc
                RM[2][0] = -sb
                RM[2][1] = cb*sc
                RM[2][2] = cacc

                # Box limits
                xr,yr,zr = rotationMatrix(RM,CO[0])
                px = xr + float(RA[0]) + float(RES)
                py = yr + float(RA[0]) + float(RES)
                pz = zr + float(RA[0]) + float(RES)
                mx = xr - RA[0] - RES
                my = yr - RA[0] - RES
                mz = zr - RA[0] - RES
                for a in range(1, CO.shape[0]):
                    xr,yr,zr = rotationMatrix(RM,CO[a])
                    x = xr + float(RA[a]) + float(RES)
                    y = yr + float(RA[a]) + float(RES)
                    z = zr + float(RA[a]) + float(RES)
                    if x > px: px = x
                    if y > py: py = y
                    if z > pz: pz = z
                    x = xr - float(RA[a]) - float(RES)
                    y = yr - float(RA[a]) - float(RES)
                    z = zr - float(RA[a]) - float(RES)
                    if x < mx: mx = x
                    if y < my: my = y
                    if z < mz: mz = z
                hx = (mx + px) / 2.0
                hy = (my + py) / 2.0
                hz = (mz + pz) / 2.0
                BOX[0] =  hx,my,pz
                BOX[1] =  px,hy,pz
                BOX[2] =  hx,py,pz
                BOX[3] =  mx,hy,pz
                BOX[4] =  mx,my,hz
                BOX[5] =  px,my,hz
                BOX[6] =  mx,py,hz
                BOX[7] =  px,py,hz
                BOX[8] =  px,hy,mz
                BOX[9] =  hx,my,mz
                BOX[10] = mx,hy,mz
                BOX[11] = hx,py,mz

                # Calculate energies
                # Loop over each boxpoint (12 points)
                v = 0.0
                for bp in range(len(BOX)):

                    # Loop over each atom
                    for a in range(len(CO)):
            
                        # Rotate atom
                        xr,yr,zr = rotationMatrix(RM,CO[a])

                        # Distance between boxpoint and atom
                        d = ((BOX[bp][0]-xr)**2 + (BOX[bp][1]-yr)**2 + (BOX[bp][2]-zr)**2)**0.5
                        
                        # Energy
                        v += PROP[a][prop] * PROBES[probe][bp] / d
                        
                # Store in array
                index = PROBES.shape[0] * prop + probe
                cuda.atomic.min(ENERGY, index, v)


# In[45]:


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
    # self.ENERGY[number of probes][number of properties]
    # #####################################
    #
    # self.ENERGY[i][j]     Energy for probe i and property j
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
    
    def __printParams(self):
        print()
        print("SETTINGS:")
        print("Accuracy:      %-3d" % (self.accuracy()))
        print("Resolution:    %-4.2f" % (self.resolution()))
        print("Normalization: %s" % (self.normalization()))
        print("Stereo:        %s" % (self.stereo()))
        print()
        
    
    
    def __init__(self, resolution=3.0, accuracy=20, stereo='none', normalization='all'):

        # Initiate PARMS
        self.PARMS = np.array([
            3,      #  0 Normalisation
            20,     #  1 Accuracy
            0,      #  2 Stereo
            0,      #  3 BeginProbe
            12,     #  4 EndProbe
            4*12,   #  5 Size of spectrophore
            12,     #  6 Number of probes
            ])

        # Initiate PROBES
        self.PROBES_TEMPLATE = np.array([
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
        print("Probes initialised: %d number of probes in total" % (len(self.PROBES_TEMPLATE)))

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

        # Initiate stereo
        if   stereo.lower() == 'none':   self.PARMS[2:7] = [0, 0,12,4*12,12]
        elif stereo.lower() == 'unique': self.PARMS[2:7] = [1,12,30,4*18,18]
        elif stereo.lower() == 'mirror': self.PARMS[2:7] = [2,30,48,4*18,18]
        elif stereo.lower() == 'all':    self.PARMS[2:7] = [3,12,48,4*36,36]
        else: raise ValueError('The stereo flag should be "none", "unique", "mirror" or "all"')
        self.PROBES = self.PROBES_TEMPLATE[self.PARMS[3]:self.PARMS[4]]
        self.PROBES = self.PROBES.astype(np.float32)
        print("%d probes are used due to the imposed stereo flag" % (self.PARMS[6]))

        # Setup the box
        self.__printParams()



        
    ####################################################
    def settings(self):
        self.__printParams()

        
        

    ####################################################
    def resolution(self, resolution=None):
        if resolution is None: return self.RESOLUTION
        elif resolution > 0:
            self.RESOLUTION = resolution
            self.__printParams()
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
            self.__printParams()




    ####################################################
    def accuracy(self, accuracy=None):
        if accuracy is None: return self.PARMS[1]
        elif (180 % accuracy) == 0:
            self.PARMS[1] = int(accuracy)
            self.__printParams()
        else: raise ValueError('(180 modus accuracy) should be equal to 0')



    ####################################################
    def stereo(self, stereo=None):
        if stereo is None:
            if   self.PARMS[2] == 0: return 'none'
            elif self.PARMS[2] == 1: return 'unique'
            elif self.PARMS[2] == 2: return 'mirror'
            elif self.PARMS[2] == 3: return 'all'
        else:
            if stereo.lower() == 'none':
                self.PARMS[2:7] = [0, 0,12,4*12,12]
            elif stereo.lower() == 'unique':
                self.PARMS[2:7] = [1,12,30,4*18,18]
            elif stereo.lower() == 'mirror':
                self.PARMS[2:7] = [2,30,48,4*18,18]
            elif stereo.lower() == 'all':
                self.PARMS[2:7] = [3,12,48,4*36,36]
            else:
                raise ValueError('The stereo flag should be "none", "unique", "mirror" or "all"')
            self.PROBES = self.PROBES_TEMPLATE[self.PARMS[3]:self.PARMS[4]]
            print("%d probes are used due to the imposed stereo flag" % (self.PARMS[6]))
            self.__printParams()





    ####################################################
    def calculate(self, mol, confID=0):
        
        ###################################
        #self.__printParams()
        ###################################

        # Wrong conformation flag
        wrong3d = False

        # Check number of atoms after adding the hydrogens
        nAtoms = mol.GetNumAtoms()
        if nAtoms < 3: raise ValueError( '>=3 atoms are needed in molecule, only %d given' % (nAtoms))

        # Create the PROP and RADIUS array
        self.PROP = np.zeros(nAtoms * 4).reshape(nAtoms, 4)
        self.PROP = self.PROP.astype(np.float32)
        self.RADIUS = np.zeros(nAtoms, dtype=np.float32)
       
        # Create the COORD array
        self.COORD_ORI = np.zeros(nAtoms * 3, dtype=np.float32).reshape(nAtoms, 3)
 
        # Atomic properties
        # [0]: atomic partial charges -> conformation dependent
        # [1]: atomic lipophilicities
        # [2]: atomic shape deviations -> conformation dependent
        # [3]: atomic electrophilicities -> conformation dependent

        chi = np.zeros(nAtoms, dtype=np.float32)
        eta = np.zeros(nAtoms, dtype=np.float32)
        A = np.zeros((nAtoms + 1, nAtoms + 1), dtype=np.float32)
        B = np.zeros(nAtoms + 1, dtype=np.float32)

        a = 0
        for atom in mol.GetAtoms():
            n = atom.GetAtomicNum()
            if   n ==  1:   # H
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
                    prop[a][1] = -0.175
            elif n ==  3:   # Li
                self.RADIUS[a] = 1.82
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n ==  5:   # B
                self.RADIUS[a] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.32966
                self.PROP[a][1] = -0.175
            elif n ==  6:   # C
                self.RADIUS[a] = 1.70
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = +0.271
            elif n ==  7:   # N
                self.RADIUS[a] = 1.55
                eta[a] = +0.34519
                chi[a] = +0.49279
                self.PROP[a][1] = -0.137
            elif n ==  8:   # O
                self.RADIUS[a] = 1.52
                eta[a] = +0.54428
                chi[a] = +0.73013
                self.PROP[a][1] = -0.321
            elif n ==  9:   # F
                self.RADIUS[a] = 1.47
                eta[a] = +0.72664
                chi[a] = +0.72052
                self.PROP[a][1] = +0.217
            elif n == 11:   # Na
                self.RADIUS[a] = 2.27
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 12:   # Mg
                self.RADIUS[a] = 1.73
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 14:   # Si
                self.RADIUS[a] = 2.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 15:   # P
                self.RADIUS[a] = 1.80
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 16:   # S
                self.RADIUS[a] = 1.80
                eta[a] = +0.20640
                chi[a] = +0.62020
                self.PROP[a][1] = +0.385
            elif n == 17:   # Cl
                self.RADIUS[a] = 1.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = +0.632
            elif n == 19:   # K
                self.RADIUS[a] = 2.75
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 20:   # Ca
                self.RADIUS[a] = 2.00
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 26:   # Fe
                self.RADIUS[a] = 1.10
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 29:   # Cu
                self.RADIUS[a] = 1.40
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 30:   # Zn
                self.RADIUS[a] = 1.39
                eta[a] = +0.32966
                chi[a] = +0.36237
                self.PROP[a][1] = -0.175
            elif n == 35:   # Br
                self.RADIUS[a] = 1.85
                eta[a] = +0.54554
                chi[a] = +0.70052
                self.PROP[a][1] = +0.815
            elif n == 53:   # I
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

        # Conformers
        if mol.GetNumConformers() < confID + 1:
            raise ValueError(
                'At least %d conformation(s) should be present, %d found' %
                (confID + 1, mol.GetNumConformers()))
        conf = mol.GetConformer(confID)

        # Coordinates
        for r in range(nAtoms):
            c = conf.GetAtomPosition(r)
            for i in range(3): self.COORD_ORI[r][i] = c[i]
            A[r][r] = 2 * eta[r]

        # Complete A matrix
        for r in range(nAtoms):
            for i in range(r + 1, nAtoms):
                d = np.sqrt(np.sum((self.COORD_ORI[r] - self.COORD_ORI[i])**2))
                if d == 0: return(np.zeros(self.PARMS[6] * 4))
                A[r][i] = 0.529176 / d    # Angstrom to au
                A[i][r] = A[r][i]

        # Property [0]: partial atomic charges
        A[:-1,nAtoms] = -1
        A[nAtoms,:-1] = +1
        A[nAtoms][nAtoms] = 0
        B[:-1] = -chi
        B[nAtoms] = Chem.GetFormalCharge(mol)
        X = scipy.linalg.solve(A, B)
        chi2 = X[nAtoms] * X[nAtoms]
        self.PROP[:,0] = X[:-1]

        # Property [2]: atomic shape deviations
        cog = np.mean(self.COORD_ORI,0)
        d = np.sqrt(np.sum((self.COORD_ORI - cog)**2, 1))
        avg_d = np.average(d)
        std_d = np.std(d)
        self.PROP[:,2] = avg_d + ((d - avg_d) / std_d)

        # Property [3]: atomic electrophilicities
        B = np.ones(nAtoms + 1)
        B[nAtoms] = 0
        A[:-1,nAtoms] = 0
        A[nAtoms,:-1] = 1
        A[nAtoms][nAtoms] = -1
        X = scipy.linalg.solve(A, B)
        for a in range(nAtoms): self.PROP[a][3] = X[a] * chi2
        # NOG TE VEREENVOUDIGEN

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
        cog = np.mean(self.COORD_ORI,0)
        self.COORD_ORI -= cog
             
        # Calculate energies
        self.TPB = (16,4)
        bpg0 = int(math.ceil(len(self.PROBES) / self.TPB[0]))
        bpg1 = int(math.ceil(len(self.PROP) / self.TPB[1]))
        self.BPG = (bpg0, bpg1)
        self.GPU_CO = cuda.to_device(self.COORD_ORI)
        self.GPU_RADIUS = cuda.to_device(self.RADIUS)
        self.GPU_PROBES = cuda.to_device(self.PROBES)
        self.GPU_PROPS = cuda.to_device(self.PROP)
        self.ENERGY = cuda.to_device(np.zeros((4 * self.PARMS[6]), dtype=np.float32))
        
        calculateEnergy[self.BPG, self.TPB](
            np.int32(self.PARMS[1]), 
            self.GPU_CO, 
            self.GPU_RADIUS, 
            self.RESOLUTION,
            self.GPU_PROBES,
            self.GPU_PROPS,
            self.ENERGY)
        self.SPHORE = self.ENERGY.copy_to_host() * -100.0
 
        # Normalise
        if self.PARMS[0] == 0: return(self.SPHORE)
        else:
            t = self.SPHORE.reshape(4, self.PARMS[6])
            m = np.mean(t,1)
            s = np.std(t,1)

            if self.PARMS[0] == 1:
                for r in range(4): t[r,:] = t[r,:] - m[r]
            elif self.PARMS[0] == 2:
                for r in range(4): t[r,:] = t[r,:] / s[r]
            elif self.PARMS[0] == 3:
                for r in range(4): t[r,:] = (t[r,:] - m[r]) / s[r]
            return(t.flatten())


# In[46]:


mol = Chem.MolFromSmiles("FC(Cl)(I)Br")
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, randomSeed=1)
calculator = SpectrophoreCalculator(accuracy = 20, normalization = 'none')
spec = calculator.calculate(mol)
np.set_printoptions(precision=3, suppress=True)
print(spec)
get_ipython().run_line_magic('timeit', 'spec = calculator.calculate(mol)')


# In[27]:


spec = calculator.calculate(mol)
print(spec)
#%timeit spec = calculator.calculate(mol)


# In[ ]:




