#!/usr/bin/env python

import unittest

from spectrophore import spectrophore
from rdkit import Chem
from rdkit.Chem import AllChem

import numpy as np


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
        refe = np.array(
        [-1.10319109,-1.07507620,-1.36544508,-0.42381262, 0.59894784, 0.56994993,
         -0.55262759,-0.57715614, 0.42761176, 1.05092946, 2.21811155, 0.23175817,
         -1.63061581,-1.49627106,-0.92949955,-0.32289585, 0.48718168, 0.55379588,
         -0.51246375,-0.16559990, 0.87225420, 1.48597513, 1.52055110, 0.13758795,
         -1.52330889,-1.24327825,-0.27087133,-0.07694924, 0.31983545, 0.36887293,
         -0.88575512,-0.79157141, 0.14561522, 1.10584918, 2.16577390, 0.68578756,
         -1.66401726,-1.40570784,-0.42076108,-0.08906513, 0.00725977, 0.67331657,
         -0.78097254,-0.28257457, 0.08920660, 1.33439735, 1.94594202, 0.59297611])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_2(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=10, normalization="all", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-1.09706685,-1.07194684,-1.37179309,-0.44025359, 0.59403602, 0.57050793,
         -0.55358716,-0.59205748, 0.43895745, 1.06219895, 2.20188373, 0.25912093,
         -1.68271362,-1.55463831,-0.86471250,-0.25011503, 0.43485924, 0.53341893,
         -0.57047323,-0.12031136, 0.84752376, 1.44546425, 1.49617920, 0.28551867,
         -1.59199956,-1.30398008,-0.28732486, 0.12322232, 0.27585675, 0.23474463,
         -1.00496457,-0.56038389, 0.15567317, 0.84030459, 2.11543153, 1.00341998,
         -1.70066565,-1.44577692,-0.52701561, 0.13723203,-0.06601790, 0.63088907,
         -0.84543821,-0.15673933, 0.05657921, 1.35649453, 1.75872797, 0.80173082])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_3(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="none", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [ 2.51490484, 2.59537645, 1.76427065, 4.45944989, 7.38683701, 7.30383802,
          4.09075036, 4.02054376, 6.89643184, 8.68051745,12.02127412, 6.33585163,
          3.59122485, 4.19542870, 6.74443436, 9.47258186,13.11583520,13.41542677,
          8.62001669,10.18000672,14.84766542,17.60782202,17.76332445,11.54356774,
         13.47567667,15.82870969,23.99962104,25.62910378,28.96319414,29.37524476,
         18.83289419,19.62429819,27.49926162,35.56788651,44.47418983,32.03820547,
          1.57047230, 1.89561423, 3.13539671, 3.55291250, 3.67415958, 4.51254555,
          2.68198757, 3.30933627, 3.77730855, 5.34466806, 6.11443796, 4.41141857])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_4(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="mean", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-3.15759900e+00,-3.07712739e+00,-3.90823319e+00,-1.21305395e+00,
          1.71433317e+00, 1.63133419e+00,-1.58175347e+00,-1.65196008e+00,
          1.22392800e+00, 3.00801362e+00, 6.34877029e+00, 6.63347792e-01,
         -7.33355305e+00,-6.72934920e+00,-4.18034354e+00,-1.45219604e+00,
          2.19105731e+00, 2.49064887e+00,-2.30476121e+00,-7.44771182e-01,
          3.92288752e+00, 6.68304412e+00, 6.83854655e+00, 6.18789846e-01,
         -1.28000138e+01,-1.04469808e+01,-2.27606945e+00,-6.46586715e-01,
          2.68750365e+00, 3.09955427e+00,-7.44279631e+00,-6.65139230e+00,
          1.22357113e+00, 9.29219602e+00, 1.81984993e+01, 5.76251498e+00,
         -2.09454919e+00,-1.76940726e+00,-5.29624776e-01,-1.12108989e-01,
          9.13809562e-03, 8.47524067e-01,-9.83033918e-01,-3.55685215e-01,
          1.12287064e-01, 1.67964657e+00, 2.44941647e+00, 7.46397080e-01])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_5(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="std", stereo="none", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [0.87864881,0.90676371,0.61639482,1.55802728,2.58078774,2.55178984,
         1.42921231,1.40468377,2.40945166,3.03276937,4.19995146,2.21359808,
         0.79850899,0.93285374,1.49962525,2.10622895,2.91630648,2.98292068,
         1.91666105,2.26352490,3.30137900,3.91509993,3.94967590,2.56671275,
         1.60371843,1.88374908,2.85615599,3.05007809,3.44686277,3.49590025,
         2.24127220,2.33545591,3.27264254,4.23287650,5.29280122,3.81281489,
         1.24766371,1.50597313,2.49091989,2.82261584,2.91894074,3.58499754,
         2.13070843,2.62910640,3.00088757,4.24607832,4.85762299,3.50465708])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_6(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="unique", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-2.03921899,-1.04517342,-0.28499011,-1.29917805,-0.47926903,-1.22532822,
         -0.89542927, 0.10037913, 1.06592355, 0.19154231, 0.02370325, 0.49968536,
          1.39211832, 1.20631878, 0.05983865, 0.52852503, 0.42327740, 1.77727531,
         -2.00962557,-0.94052749,-1.05069265,-0.57476055,-0.59438612,-1.40042448,
         -0.65892811, 0.06031672, 1.12082133, 0.22931550,-0.47580130, 0.02467334,
          1.60509384, 0.99574236, 0.49127697, 0.72111261, 1.28871180, 1.16808180,
         -1.89353923,-0.37359974,-0.85411562,-0.32538559,-0.53902770,-1.41911805,
         -1.11538883, 0.53695328, 1.66420002, 0.46439284,-0.92596599, 0.23159390,
          1.91784546, 0.20877391, 0.11192487, 0.66190005, 0.52831513, 1.12024128,
         -2.14066262,-0.16366780,-1.13408345,-0.16928469,-0.69858810,-1.12002024,
         -1.27966436, 0.14968130, 0.95374096, 0.92985097,-0.76864761, 0.01253743,
          1.78150756, 0.22270588, 0.57046401, 0.91541664, 0.82401849, 1.11469565])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_7(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="mirror", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-1.90921713,-0.98252904,-0.90366253,-0.91434369,-0.74867565,-0.80129392,
         -0.95162298, 0.28901830, 1.35555346, 0.64547353,-0.54970129,-0.27716388,
          1.64987831, 0.88433158, 0.14506460, 0.87131347, 0.77244895, 1.42512791,
         -2.31299993,-0.82483996,-0.85101572,-0.52132848,-0.41838436,-1.46410667,
         -0.65312759, 0.24702886, 1.19228061, 0.12813404,-0.36024007, 0.28587791,
          1.71346271, 0.32608845, 0.47657839, 0.78201641, 0.90884304, 1.34573237,
         -2.23400905,-0.77962677,-0.73147309,-0.17699684,-0.45110149,-1.23282715,
         -0.95921007, 0.31719927, 1.71804647, 0.29338296,-0.65229095, 0.41400049,
          1.74454209, 0.28271240, 0.37067729, 0.65412822, 0.09683844, 1.32600779,
         -2.50565022,-0.70145335,-0.98193163,-0.22953142,-0.69302100,-1.01641870,
         -0.75447023, 0.30909326, 1.05843658, 0.16536189,-0.36628141, 0.89948652,
          1.78925164, 0.43424018, 0.19778535, 0.55352457, 0.41453461, 1.42704334])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_8(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="all", resolution=3.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-1.88686588,-1.02863665,-0.37231710,-1.24793665,-0.54005171,-1.18417691,
         -0.89935202,-0.03960081, 0.79402138, 0.03910675,-0.10580048, 0.30514826, 
          1.07564820, 0.91523443,-0.07460224, 0.33004757, 0.23917990, 1.40818123,
         -1.98487173,-0.96017665,-0.87296917,-0.88477998,-0.70159080,-0.75977402,
         -0.92600195, 0.44585017, 1.62518262, 0.84000428,-0.48157282,-0.18021172,
          1.95063541, 1.10412407, 0.28667186, 1.08972916, 0.98040869, 1.70211530,
         -2.04385246,-0.93105785,-1.04572572,-0.55034117,-0.57076887,-1.40975190,
         -0.63794884, 0.11069319, 1.21454310, 0.28659935,-0.44733722, 0.07359299,
          1.71860896, 1.08435187, 0.55926780, 0.79849735, 1.38929565, 1.26373522,
         -2.25678837,-0.83561989,-0.86061731,-0.54577138,-0.44746143,-1.44610916,
         -0.67163738, 0.18799734, 1.09069733, 0.07445472,-0.39193458, 0.22509754,
          1.58841773, 0.26349796, 0.40721340, 0.69890172, 0.82001908, 1.23724125,
         -1.90478601,-0.36869716,-0.85431853,-0.31997074,-0.53588280,-1.42532411,
         -1.11836776, 0.55153048, 1.69075421, 0.47819908,-0.92693231, 0.24292665,
          1.94709464, 0.21986420, 0.12198614, 0.67780480, 0.54280054, 1.14101591,
         -2.21871209,-0.78006430,-0.73243159,-0.18395404,-0.45509324,-1.22836166,
         -0.95770475, 0.30489552, 1.69058741, 0.28133687,-0.65410609, 0.40064949,
          1.71679638, 0.27078175, 0.35779498, 0.63817935, 0.08691871, 1.30279007,
         -2.13624368,-0.10707831,-1.10310208,-0.11284343,-0.65611452,-1.08866776,
         -1.25252470, 0.21453969, 1.03981754, 1.01529712,-0.72802282, 0.07377676,
          1.88942791, 0.28949131, 0.64642637, 1.00048188, 0.90667184, 1.20501963,
         -2.48907143,-0.74066938,-1.01247387,-0.28334165,-0.73249780,-1.04589442,
         -0.79204671, 0.23862606, 0.96479593, 0.09933958,-0.41586263, 0.81076140,
          1.67301052, 0.35990281, 0.13076035, 0.47549833, 0.34080663, 1.32200352])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
       
        
    def test_spectrophore_run_9(self):
        """
        Test the output data
        """
        mol = Chem.MolFromSmiles("ClC(F)(Br)I")
        AllChem.EmbedMolecule(mol, randomSeed = 1)
        calculator = spectrophore.SpectrophoreCalculator(accuracy=20, normalization="all", stereo="none", resolution=5.0)
        spec = calculator.calculate(mol)
        refe = np.array(
        [-1.17696999,-1.16754104,-1.37535841,-0.38849771, 0.62028062, 0.50084625,
         -0.46682031,-0.47735307, 0.44754847, 1.03628022, 2.18767882, 0.25990616,
         -1.70774136,-1.62467284,-0.87070400,-0.20418283, 0.64592735, 0.47249982,
         -0.51501546,-0.12306470, 0.88154932, 1.35987348, 1.42435165, 0.26117956,
         -1.64621639,-1.46510089,-0.17857477, 0.06535166, 0.33775689, 0.34624943,
         -0.85367692,-0.64890728, 0.26991872, 1.00050675, 2.02104921, 0.75164359,
         -1.76818132,-1.57525892,-0.35909169, 0.07542518, 0.09415803, 0.56238827,
         -0.72335774,-0.16622585, 0.11416589, 1.21297226, 1.85346106, 0.67954482])
        np.testing.assert_allclose(spec, refe, rtol=1e-5)
        
        
        

if __name__ == '__main__':
    unittest.main()