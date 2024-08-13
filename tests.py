#!/usr/bin/env python

import unittest

from spectrophore import spectrophore
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

TOLERANCE = 1E-5

class Tests(unittest.TestCase):
    
    def test_spectrophore_creation_default(self):
        """
        Test that a spectrophore calculator can be created
        """
        calculator = spectrophore.SpectrophoreCalculator()
        self.assertIsInstance(calculator, spectrophore.SpectrophoreCalculator)
       
        
    def test_spectrophore_creation_normalization_1(self):
        """
        Test the normalization default setting
        """
        calculator = spectrophore.SpectrophoreCalculator()
        self.assertEqual(calculator.normalization(), "all")
       
        
    def test_spectrophore_creation_normalization_2(self):
        """
        Test the normalization 'all' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(normalization='all')
        self.assertEqual(calculator.normalization(), "all")
       
        
    def test_spectrophore_creation_normalization_3(self):
        """
        Test the normalization 'none' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(normalization='none')
        self.assertEqual(calculator.normalization(), "none")
       
        
    def test_spectrophore_creation_normalization_4(self):
        """
        Test the normalization 'std' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(normalization='std')
        self.assertEqual(calculator.normalization(), "std")
       
        
    def test_spectrophore_creation_normalization_5(self):
        """
        Test the normalization 'mean' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(normalization='mean')
        self.assertEqual(calculator.normalization(), "mean")
       
        
    def test_spectrophore_creation_stereo_1(self):
        """
        Test the stereo default setting
        """
        calculator = spectrophore.SpectrophoreCalculator()
        self.assertEqual(calculator.stereo(), "none")
       
        
    def test_spectrophore_creation_stereo_2(self):
        """
        Test the stereo 'none' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(stereo="none")
        self.assertEqual(calculator.stereo(), "none")
       
        
    def test_spectrophore_creation_stereo_3(self):
        """
        Test the stereo 'unique' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(stereo="unique")
        self.assertEqual(calculator.stereo(), "unique")
       
        
    def test_spectrophore_creation_stereo_4(self):
        """
        Test the stereo 'mirror' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(stereo="mirror")
        self.assertEqual(calculator.stereo(), "mirror")
       
        
    def test_spectrophore_creation_stereo_5(self):
        """
        Test the stereo 'all' setting
        """
        calculator = spectrophore.SpectrophoreCalculator(stereo="all")
        self.assertEqual(calculator.stereo(), "all")
       
        
    def test_spectrophore_creation_resolution_1(self):
        """
        Test the resolution default setting
        """
        calculator = spectrophore.SpectrophoreCalculator()
        self.assertEqual(calculator.resolution(), 3.0)
       
        
    def test_spectrophore_creation_resolution_2(self):
        """
        Test the resolution setting
        """
        calculator = spectrophore.SpectrophoreCalculator(resolution=5.0)
        self.assertEqual(calculator.resolution(), 5.0)
       
        
    def test_spectrophore_creation_accuracy_1(self):
        """
        Test the accuracy default setting
        """
        calculator = spectrophore.SpectrophoreCalculator()
        self.assertEqual(calculator.accuracy(), 20)
       
        
    def test_spectrophore_creation_accuracy_2(self):
        """
        Test the accuracy setting
        """
        calculator = spectrophore.SpectrophoreCalculator(accuracy=30)
        self.assertEqual(calculator.accuracy(), 30)
       
        
    def test_spectrophore_run_1(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-1.1031914 , -1.0750767 , -1.3654453 , -0.42381215,  0.59894806,
        0.5699494 , -0.55262715, -0.57715636,  0.4276119 ,  1.0509295 ,
        2.2181115 ,  0.23175779, -1.6306345 , -1.4963043 , -0.92958903,
       -0.32305473,  0.4869356 ,  0.55584854, -0.5126002 , -0.16577469,
        0.8719677 ,  1.4856232 ,  1.5202019 ,  0.13738069, -1.5233064 ,
       -1.2432761 , -0.27087742, -0.07694667,  0.31983256,  0.36887285,
       -0.885755  , -0.791568  ,  0.14561523,  1.1058402 ,  2.1657825 ,
        0.6857875 , -1.6640162 , -1.4057074 , -0.42076126, -0.08906693,
        0.00726358,  0.6733182 , -0.780978  , -0.28257385,  0.0892069 ,
        1.3344021 ,  1.9459364 ,  0.59297746], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_2(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=10, normalization="all", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-1.0970671 , -1.0719464 , -1.3717929 , -0.44025284,  0.5940364 ,
        0.57050735, -0.5535874 , -0.5920576 ,  0.43895695,  1.0622    ,
        2.2018838 ,  0.25912088, -1.6827124 , -1.5546383 , -0.8647143 ,
       -0.25011283,  0.4348574 ,  0.53342175, -0.57047194, -0.12031432,
        0.847524  ,  1.4454637 ,  1.49618   ,  0.28551823, -1.5919987 ,
       -1.3039782 , -0.28732124,  0.12322673,  0.27585822,  0.23473941,
       -1.0049644 , -0.56039035,  0.15567124,  0.840308  ,  2.115431  ,
        1.0034201 , -1.7006627 , -1.4457755 , -0.527023  ,  0.13723296,
       -0.06602039,  0.63088846, -0.8454364 , -0.15673935,  0.0565788 ,
        1.3564981 ,  1.7587292 ,  0.8017287 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_3(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="none", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([ 2.5149038,  2.5953748,  1.76427  ,  4.459451 ,  7.386838 ,
        7.303837 ,  4.0907516,  4.020543 ,  6.8964324,  8.680518 ,
       12.021275 ,  6.3358507,  3.5912237,  4.195427 ,  6.744452 ,
        9.472579 , 13.115829 , 13.425793 ,  8.620024 , 10.180008 ,
       14.847663 , 17.60782  , 17.763351 , 11.543569 , 13.475672 ,
       15.828707 , 23.999565 , 25.629124 , 28.963175 , 29.37525  ,
       18.83288  , 19.624313 , 27.499264 , 35.56783  , 44.474297 ,
       32.038216 ,  1.5704744,  1.8956152,  3.135395 ,  3.5529082,
        3.6741621,  4.512544 ,  2.68198  ,  3.3093355,  3.7773066,
        5.3446693,  6.114425 ,  4.411417 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_4(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="mean", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-3.1576002e+00, -3.0771291e+00, -3.9082341e+00, -1.2130527e+00,
        1.7143340e+00,  1.6313329e+00, -1.5817523e+00, -1.6519608e+00,
        1.2239285e+00,  3.0080142e+00,  6.3487706e+00,  6.6334677e-01,
       -7.3344212e+00, -6.7302179e+00, -4.1811929e+00, -1.4530659e+00,
        2.1901846e+00,  2.5001478e+00, -2.3056211e+00, -7.4563694e-01,
        3.9220181e+00,  6.6821756e+00,  6.8377066e+00,  6.1792374e-01,
       -1.2800018e+01, -1.0446983e+01, -2.2761250e+00, -6.4656639e-01,
        2.6874847e+00,  3.0995598e+00, -7.4428101e+00, -6.6513767e+00,
        1.2235737e+00,  9.2921391e+00,  1.8198606e+01,  5.7625256e+00,
       -2.0945449e+00, -1.7694041e+00, -5.2962422e-01, -1.1211109e-01,
        9.1428757e-03,  8.4752488e-01, -9.8303938e-01, -3.5568380e-01,
        1.1228728e-01,  1.6796501e+00,  2.4494059e+00,  7.4639773e-01],
      dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_5(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="std", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([0.87864834, 0.906763  , 0.6163945 , 1.5580276 , 2.580788  ,
       2.551789  , 1.4292126 , 1.4046834 , 2.4094517 , 3.0327692 ,
       4.199951  , 2.2135975 , 0.7984234 , 0.9327537 , 1.4994689 ,
       2.1060033 , 2.9159937 , 2.9849064 , 1.9164578 , 2.2632833 ,
       3.3010256 , 3.9146812 , 3.9492598 , 2.5664387 , 1.6037147 ,
       1.8837451 , 2.8561437 , 3.0500743 , 3.4468536 , 3.495894  ,
       2.241266  , 2.335453  , 3.2726364 , 4.2328615 , 5.2928033 ,
       3.8128085 , 1.2476672 , 1.5059761 , 2.4909222 , 2.8226166 ,
       2.918947  , 3.5850017 , 2.1307054 , 2.6291096 , 3.0008903 ,
       4.2460856 , 4.85762   , 3.5046608 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_6(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="unique", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-2.0392199 , -1.0451725 , -0.28498945, -1.2991787 , -0.47926974,
       -1.2253278 , -0.8954297 ,  0.1003789 ,  1.0659237 ,  0.19154277,
        0.02370401,  0.499685  ,  1.3921181 ,  1.2063197 ,  0.0598392 ,
        0.52852434,  0.42327705,  1.777275  , -1.9831187 , -0.93960583,
       -1.0471339 , -0.58259386, -0.6017486 , -1.3884957 , -0.664746  ,
        0.03728981,  1.0724077 ,  0.20224081, -0.45646247,  0.00250024,
        1.5450943 ,  0.95032066,  0.47626516,  0.68226475,  1.3573786 ,
        1.3381474 , -1.8927505 , -0.37909365, -0.8576261 , -0.3310822 ,
       -0.54383785, -1.4202882 , -1.1178166 ,  0.52769566,  1.6502874 ,
        0.45543864, -0.9291761 ,  0.2235955 ,  1.9028864 ,  0.20087405,
        0.10442823,  0.65212655,  0.64576787,  1.1085755 , -2.140659  ,
       -0.16367035, -1.1340805 , -0.1692811 , -0.6985873 , -1.1200222 ,
       -1.2796654 ,  0.14967908,  0.9537448 ,  0.92984784, -0.7686497 ,
        0.01253645,  1.7815126 ,  0.22271025,  0.5704576 ,  0.9154165 ,
        0.8240188 ,  1.1146964 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_7(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="mirror", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-1.9092162 , -0.9825288 , -0.90366286, -0.9143438 , -0.7486751 ,
       -0.8012937 , -0.9516235 ,  0.28901854,  1.3555537 ,  0.6454737 ,
       -0.549702  , -0.2771643 ,  1.6498784 ,  0.88433164,  0.14506489,
        0.8713128 ,  0.7724476 ,  1.4251299 , -2.2960248 , -0.8243914 ,
       -0.8502815 , -0.52425385, -0.4224533 , -1.4565605 , -0.6545931 ,
        0.23557016,  1.1703192 ,  0.11799809, -0.36495772,  0.27398908,
        1.6857136 ,  0.3137505 ,  0.48427516,  0.76461554,  0.89003146,
        1.4572529 , -2.2340105 , -0.7796255 , -0.73147553, -0.17699519,
       -0.45110518, -1.2328252 , -0.95921105,  0.31719577,  1.718044  ,
        0.29338098, -0.65229195,  0.41400814,  1.7445465 ,  0.28271115,
        0.37067434,  0.65412444,  0.09685306,  1.3260028 , -2.5056496 ,
       -0.70145863, -0.9819242 , -0.22952928, -0.69301647, -1.0164218 ,
       -0.7544724 ,  0.30910036,  1.0584307 ,  0.16536063, -0.36628118,
        0.8994829 ,  1.7892605 ,  0.434241  ,  0.19778353,  0.55352926,
        0.41452837,  1.4270407 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_8(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="all", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array([-1.8868657 , -1.0286354 , -0.37231636, -1.2479366 , -0.54005206,
       -1.1841761 , -0.89935195, -0.03960101,  0.7940211 ,  0.03910712,
       -0.10579979,  0.3051478 ,  1.0756475 ,  0.91523474, -0.07460175,
        0.33004677,  0.23917948,  1.4081802 , -1.9848713 , -0.9601767 ,
       -0.8729698 , -0.88478035, -0.7015905 , -0.759774  , -0.9260028 ,
        0.44585046,  1.6251832 ,  0.8400046 , -0.48157376, -0.18021232,
        1.9506359 ,  1.1041243 ,  0.28667217,  1.0897286 ,  0.9804073 ,
        1.7021178 , -2.0215762 , -0.92931026, -1.041862  , -0.55561876,
       -0.5756684 , -1.3991723 , -0.641609  ,  0.09322588,  1.1767044 ,
        0.2658834 , -0.4235945 ,  0.05681095,  1.6714748 ,  1.0489135 ,
        0.5527102 ,  0.76833403,  1.4749892 ,  1.4548595 , -2.2305825 ,
       -0.83562994, -0.86017096, -0.5511313 , -0.45463514, -1.4348592 ,
       -0.67467904,  0.16910143,  1.0551445 ,  0.05765557, -0.4001354 ,
        0.20551851,  1.543684  ,  0.24320813,  0.40484744,  0.6705804 ,
        0.78946143,  1.3271272 , -1.9042444 , -0.3714665 , -0.85604393,
       -0.3228485 , -0.5382918 , -1.425814  , -1.1195213 ,  0.54677784,
        1.6835508 ,  0.47360805, -0.92849785,  0.23883615,  1.9393407 ,
        0.21582766,  0.1181635 ,  0.67278063,  0.6663416 ,  1.1349957 ,
       -2.217496  , -0.7819456 , -0.7344192 , -0.18711956, -0.45767972,
       -1.2292762 , -0.9592055 ,  0.30067146,  1.6833782 ,  0.27716509,
       -0.6562611 ,  0.3962301 ,  1.7095375 ,  0.26663342,  0.35345745,
        0.6332368 ,  0.08318226,  1.2964141 , -2.1362412 , -0.10708191,
       -1.1031003 , -0.11284073, -0.6561148 , -1.0886708 , -1.2525269 ,
        0.21453649,  1.0398207 ,  1.0152931 , -0.72802603,  0.0737748 ,
        1.8894325 ,  0.28949487,  0.64641887,  1.0004809 ,  0.90667135,
        1.2050197 , -2.48907   , -0.74067366, -1.0124658 , -0.2833388 ,
       -0.7324926 , -1.0458965 , -0.792048  ,  0.23863365,  0.9647909 ,
        0.0993391 , -0.4158616 ,  0.81075853,  1.6730196 ,  0.35990432,
        0.13075931,  0.4755036 ,  0.3408013 ,  1.3220016 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
       
        
    def test_spectrophore_run_9(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="none", resolution=5.0)
        spec = calculator.calculate(mol)
        refe = np.array([-1.1769704 , -1.1675407 , -1.3753582 , -0.3884973 ,  0.62028074,
        0.50084645, -0.46682152, -0.477354  ,  0.44754717,  1.0362803 ,
        2.187679  ,  0.25990626, -1.7077427 , -1.6246742 , -0.87070256,
       -0.20418099,  0.6459274 ,  0.47249973, -0.5150156 , -0.12306376,
        0.88154876,  1.3598802 ,  1.424343  ,  0.2611806 , -1.6462182 ,
       -1.4651036 , -0.1785752 ,  0.06534798,  0.33775648,  0.34625006,
       -0.85368073, -0.64889866,  0.2699157 ,  1.000516  ,  2.0210416 ,
        0.75164795, -1.7681832 , -1.5752606 , -0.35909587,  0.07542755,
        0.0941598 ,  0.56239283, -0.7233517 , -0.16622925,  0.11416941,
        1.2129686 ,  1.8534594 ,  0.6795451 ], dtype=np.float32)
        np.testing.assert_allclose(spec, refe, rtol=TOLERANCE)
        

    def test_spectrophore_valueerror_natoms(self):
        """
        Test the error raised when less than 3 atoms are present in the molecule
        """
        mol = Chem.MolFromSmiles("Br")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="none", resolution=5.0)
        with self.assertRaises(ValueError) as error:
            spec = calculator.calculate(mol)
        self.assertEqual(str(error.exception), ">=3 atoms are needed in molecule, only 2 given")


    def test_spectrophore_valueerror_conformers(self):
        """
        Test the error raised when no conformation is available
        """
        mol = Chem.MolFromSmiles("c1ccc(CCN)cc1")
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="none", resolution=5.0)
        with self.assertRaises(ValueError) as error:
            spec = calculator.calculate(mol)
        self.assertEqual(str(error.exception), "At least 1 conformation(s) should be present, 0 found")


if __name__ == '__main__':
    unittest.main()
