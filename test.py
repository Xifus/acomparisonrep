import unittest
import numpy as np
import a_comparison as acn
import sampling_method as sm
import math as math


class TestAComparison(unittest.TestCase):
    
    def test_testfunction(self):
        #k needs to be the length of X and a. The length of X and a should be equal
        X = np.array([0, 1, 2])
        a = np.array([0, 1, 2])
        b = np.array([0, 1])
        #a cannot be -1
        c = np.array([-1, -1, -1])
        self.assertEqual(acn.testfunction(3, X, a), 8.)
        self.assertFalse(acn.testfunction(2, X, a))
        self.assertFalse(acn.testfunction(3, X, b))
        self.assertFalse(acn.testfunction(3, X, c))
        
    def test_vana(self):
        self.assertFalse(acn.vana(-1))
        self.assertEqual(acn.vana(1.), 1./12.)
        
    def test_analyticalvalues(self):
        a = np.array([0, 1])
        b = np.array([-1., -1.])
        v1 = np.array([1./3., 1./12.])
        v2 = np.array([13./36., 1./9.])
        vv = 4./9.
        self.assertFalse(acn.analyticalvalues(1, a))
        self.assertFalse(acn.analyticalvalues(-1, a))
        self.assertFalse(acn.analyticalvalues(2, b))
        v1_e, v2_e, vv_e = acn.analyticalvalues(2, a)
        for i in range(2):
            self.assertAlmostEqual(v1[i], v1_e[i])
            self.assertAlmostEqual(v2[i], v2_e[i])
        self.assertAlmostEqual(vv, vv_e)
        
    def test_analyticalSAvalues(self):
        a = np.array([0, 1])
        S1 = np.array([3./4., 3./16.])
        S2 = np.array([39./48., 1./4.])
        S1_e, S2_e = acn.analyticalSAvalues(2, a)
        for i in range(2):
            self.assertAlmostEqual(S1[i], S1_e[i])
            self.assertAlmostEqual(S2[i], S2_e[i])
        
    def test_a_e_s_st(self):
        S1 = np.array([[1, 2], [3, 4]])
        S_ana = np.array([1, 2])
        AES = np.array([1, 1])
        AES_e = acn.a_e_s_st(2, S1, S_ana, 2)
        self.assertAlmostEqual(AES[0], AES_e[0])
        self.assertAlmostEqual(AES[1], AES_e[1])
        
    def test_m_a_e_s_st(self):
        a = np.array([0, 1])
        self.assertFalse(acn.m_a_e_s_st(1, a))
        self.assertAlmostEqual(acn.m_a_e_s_st(2, a), math.log(1./2.))

    def test_e_var(self):
        X = np.array([[1, 1], [2, 2]])
        a = np.array([1, 2])
        result = 44./3.
        self.assertAlmostEqual(acn.e_var(2, X, a, 1), result)
        
    def test_e_si(self):
        X = np.array([[1, 1], [2, 2],[3, 3], [4, 4],[3, 1], [4, 2],[1, 3], [2, 4]])
        a = np.array([1, 2])
        result = np.array([136., 346./3.])
        result_e = acn.e_si(2, X, a, 2, 2)
        self.assertAlmostEqual(result[0], result_e[0])
        self.assertAlmostEqual(result[1], result_e[1])
        
    def test_e_sti_1(self):
        X = np.array([[1, 1], [2, 2],[3, 3], [4, 4],[3, 1], [4, 2],[1, 3], [2, 4]])
        a = np.array([1, 2])
        result = np.array([160./9., 116./9.])
        result_e = acn.e_sti_1(2, X, a, 2, 2)
        self.assertAlmostEqual(result[0], result_e[0])
        self.assertAlmostEqual(result[1], result_e[1])
        
    def test_e_sti_2(self):
        X = np.array([[1, 1], [2, 2],[3, 3], [4, 4],[3, 1], [4, 2],[1, 3], [2, 4]])
        a = np.array([1, 2])
        result = np.array([- 248./9., - 214./9.])
        result_e = acn.e_sti_2(2, X, a, 2, 2)
        self.assertAlmostEqual(result[0], result_e[0])
        self.assertAlmostEqual(result[1], result_e[1])
        
    def test_less_eq_integer_matrix(self):
        X = np.array([[1, 1], [2, 2]])
        z = np.array([[1.1, 1.1], [2.1, 2.1]])
        result_1 = acn.less_eq_integer_matrix(X, X)
        result_2 = acn.less_eq_integer_matrix(X, z)
        result = np.array([[2, 2], [4, 4]])
        for i in range(2):
            for j in range(2):
                self.assertEqual(result[i,j], result_1[i,j])
                self.assertEqual(result[i,j], result_2[i,j])
        
    def test_shift(self):
        X = np.array([[1, 1], [2, 2]])
        z = np.array([[1.1, 1.1], [2.1, 2.1]])
        result = np.array([[.1, .1], [.1, .1]])
        result_1 = acn.shift(X, z)
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(result[i,j], result_1[i,j])
        
class TestSamplingMethod(unittest.TestCase):
    
    def test_create_sample_matrix(self):
        X = np.array([[1, 1, 3, 4,], [2, 2, 4, 4]])
        result = np.array([[1, 1], [2, 2],[3, 4], [4, 4],[3, 1], [4, 2],[1, 4], [2, 4]])
        result_1 = sm.create_sample_matrix(X)
        for i in range(8):
            for j in range(2):
                self.assertEqual(result[i,j], result_1[i,j])
    
    def test_choose_sampling_method(self):
        N = 2
        k = 2
        lower = 0
        upper = 1
        result_1 = sm.choose_sampling_method(N, k, "S")
        result_2 = sm.choose_sampling_method(N, k, "L")
        row, col = result_1.shape
        for i in range(row):
            for j in range(col):
                self.assertGreaterEqual(result_1[i, j], lower)
                self.assertLessEqual(result_1[i, j], upper)
                self.assertGreaterEqual(result_2[i, j], lower)
                self.assertLessEqual(result_2[i, j], upper)
                
    def test_create_coefficient_a(self):
        k = 10
        a1 = np.array([0, 0, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52])
        a2 = np.array([6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 0, 0])
        a3 = np.array([0, 0, 0, 0, 0, 6.52, 6.52, 6.52, 6.52, 6.52])
        a4 = np.array([6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52, 6.52])
        a1_e = sm.create_coefficient_a(k, "A1-1")
        a2_e = sm.create_coefficient_a(k, "A1-2")
        a3_e = sm.create_coefficient_a(k, "A2")
        a4_e = sm.create_coefficient_a(k, "B")
        for i in range(k):
            self.assertEqual(a1[i], a1_e[i])
            self.assertEqual(a2[i], a2_e[i])
            self.assertEqual(a3[i], a3_e[i])
            self.assertEqual(a4[i], a4_e[i])
            
if __name__ == '__main__':
    unittest.main()