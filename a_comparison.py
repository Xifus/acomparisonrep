import numpy as np
import math as math

def testfunction(k, X, a, rule = "V-function"):
    #build test function
    if rule == "V-function":
        if k != len(X) or k != len(a) or len(X) != len(a):
            return False
        Y = 1.
        for i in range(k):
            if a[i] + 1. == 0:
                return False
            v = (abs(4. * X[i] - 2.) + a[i])/(1. + a[i])
            Y = Y * v
        return Y

def vana(a):
    #analytical value for V_i
    if a + 1.0 == 0:
        return False
    return (1./3.)/(math.pow((1. + a),2))

def analyticalvalues(k, a):
    if k <= 0 or k != len(a):
        print "something is wrong"
        return False
    #v-function analytical values
    var = 1.
    v_ana = np.zeros((k))
    vt_ana = np.ones((k))
    #define v_ana, vt_ana array
    for i in range(k):
        if a[i] + 1.0 == 0:
            return False
        v_ana[i] = vana(a[i])
        var = var * (1. + v_ana[i])
        for j in range(k):
            if j != i:
                vt_ana[i] = vt_ana[i] * (1. + vana(a[j]))
        vt_ana[i] = vt_ana[i]*v_ana[i]
    return v_ana, vt_ana, (var - 1.)

def analyticalSAvalues(k, a):
    #calculate analytical sensitivity indices
    V_ana, Vt_ana, var = analyticalvalues(k, a)
    S_ana = V_ana/var
    St_ana = Vt_ana/var
    return S_ana, St_ana

#def functions for convergence
def a_e_s_st(R, S, S_ana, k):
    #mean absolute error of each S_i/S_ti, return an array
    AES = np.zeros((k))
    for i in range(k):
        for r in range(R):
            AES[i] = abs(S[r, i] - S_ana[i]) + AES[i]
        AES[i] = AES[i]/R
    return AES

def m_a_e_s_st(k, AES):
    #average of the AES/AEST values over all inputs to give a single error measure
    #return a number
    if k <= 0 or k != len(AES):
        return False
    MAES = 0.
    MAES = float(np.sum(AES))
    return math.log(MAES/k)
    
def e_var(k, X, a, N, rule = "V-function"):
    #estimated variance
    if k != len(a):
        return False
    A = 0.
    B = 0.
    for j in range(N):
        f_A = testfunction(k, X[j, :], a, rule)
        f_B = testfunction(k, X[N + j, :], a, rule)
        A = math.pow(f_A, 2) + A
        B = f_A * f_B + B
    return abs(A/N - B/N)
        
def e_si(k, X, a, N, EVar, rule = "V-function"):
    #estimated S_i, return an array
    S_i = np.zeros((k))
    for i in range(k):
        AB = 0.
        for j in range(N):
            f_B = testfunction(k, X[N + j, :], a, rule)
            f_AB = testfunction(k, X[(2 + i) * N + j, :], a, rule)
            f_A = testfunction(k, X[j, :], a, rule)
            AB = f_B * (f_AB - f_A) + AB
        S_i[i] = AB/(N * EVar)
    return S_i
        
def e_sti_1(k, X, a, N, EVar, rule = "V-function"):
    #estimated S_ti by first way, return an array
    S_ti = np.zeros((k))
    for i in range(k):
        AB = 0.
        for j in range(N):
            f_A = testfunction(k, X[j, :], a, rule)
            f_AB = testfunction(k, X[(2 + i) * N + j, :], a, rule)
            AB = math.pow((f_A - f_AB), 2) + AB
        S_ti[i] = AB/(2. * N * EVar)
    return S_ti

def e_sti_2(k, X, a, N, Evar, rule = "V-function"):
    #estimated S_ti by second way, return an array
    S_ti = np.zeros((k))
    for i in range(k):
        AB = 0.
        for j in range(k):
            f_AB = testfunction(k, X[(2 + i) * N + j, :], a, rule)
            f_A = testfunction(k, X[j, :], a, rule)
            AB = f_A * (f_A - f_AB) + AB
        S_ti[i] = AB/(N * Evar)
    return S_ti
     
def less_eq_integer_matrix(X, z):
    N, k = X.shape
    Y = np.zeros((N, k))
    for i in range(N):
        for j in range(k):
            Y[i, j] = math.floor(X[i, j] + z[i, j])
    return Y
       
def shift(X, z):
    #shift the matrix
    N, k = X.shape
    Y = np.zeros((N, k))
    Y = less_eq_integer_matrix(X, z)
    X = X + z - Y
    #X = X + z + Y
    return X
        
        