import a_comparison as acn
import numpy as np
import matplotlib.pyplot as plt
import sampling_method as sm

def n_calculate(k, N, a, R, rl, rule = "V-function"):
    leng = len(N)
    #er_si = np.zeros((R, k))
    #er_sti = np.zeros((R, k))
    er_si_N = np.zeros(((leng * R), k))
    er_sti_N = np.zeros(((leng * R), k))
    for i in range(leng):
        c = i*R
        er_si, er_sti = calculate(k, N[i], a, R, rl, rule)
        er_si_N[c:c+R, :] = er_si
        er_sti_N[c:c+R, :] = er_sti
    return er_si_N, er_sti_N

def calculate(k, N, a, R, rl = "S", rule = "V-function"):
    #sampling by chaospy method. rule = sampling method, R = replicated times
    er_si = np.zeros((R, k))
    er_sti = np.zeros((R, k))
    
    
    for i in range(R):
        #create a Nx2k sample matrix
        x = sm.choose_sampling_method(N, 2*k, rl)
    
        #create the A, B, A_i_B matrix
        X = sm.create_sample_matrix(x)
        rows, cols = X.shape
    
        #generate a shift matrix
        #z = cp_sampling(lower, upper, rows, cols, rl)
        z_1 = sm.choose_sampling_method(1, cols, rl)
        z = np.tile(z_1, (rows, 1))
        X_shift = acn.shift(X, z)
        
        estimated_var = acn.e_var(k, X_shift, a, N, rule)
        
        estimated_si = np.zeros((k))
        estimated_sti_1 = np.zeros((k))
        estimated_si = acn.e_si(k, X_shift, a, N, estimated_var, rule)
        estimated_sti_1 = acn.e_sti_1(k, X_shift, a, N, estimated_var, rule)
        #estimated_Sti_2 = acn.e_sti_2(k, X_shift, a, N, estimated_var, rule)
        
        er_si[i, :] = estimated_si
        er_sti[i, :] = estimated_sti_1
               
    return er_si, er_sti

def aes_aest_N(N, k, R, a, er_si_N, er_sti_N):
    s_ana, st_ana = acn.analyticalSAvalues(k, a)
    leng = len(N)
    aes_N = np.zeros((leng, k))
    aest_N = np.zeros((leng, k))
    for i in range(leng):
        y = i*R
        z = (i+1)*R
        e_si = er_si_N[y:z,:]
        e_sti = er_sti_N[y:z,:]
        aes, aest = aes_aest(R, e_si, e_sti, s_ana, st_ana, k)
        
        aes_N[i, :] = aes
        aest_N[i, :]= aest
    return aes_N, aest_N

def maes_maest_N(N, k, aes_N, aest_N):
    leng = len(N)
    maes = np.zeros(leng)
    maest = np.zeros(leng)
    for i in range(leng):
        #aes = aes_N[i, :].flatten()
        #aest = aest_N[i, :].flatten()
        aes = aes_N[i, :]
        aest = aest_N[i, :]
        maes[i] = acn.m_a_e_s_st(k, aes)
        maest[i] = acn.m_a_e_s_st(k, aest)
    return maes, maest
        
        

def aes_aest(R, e_si, e_sti, s_ana, st_ana, k):
    aes = acn.a_e_s_st(R, e_si, s_ana, k)
    aest = acn.a_e_s_st(R, e_sti, st_ana, k)
    return aes, aest

def draw_aes_plot(k, a, R, e_si_1, e_si_2, e_si_3, e_sti_1, e_sti_2, e_sti_3):
    s_ana, st_ana = acn.analyticalSAvalues(k, a)
    #ESi_1, ESti_1 = calculate(k, N, a, R, rl1, lower, upper)
    #ESi_2, ESti_2 = calculate(k, N, a, R, rl2, lower, upper)
    aes_1, aest_1 = aes_aest(R, e_si_1, e_sti_1, s_ana, st_ana, k)
    aes_2, aest_2 = aes_aest(R, e_si_2, e_sti_2, s_ana, st_ana, k)
    aes_3, aest_3 = aes_aest(R, e_si_3, e_sti_3, s_ana, st_ana, k)
    
    kk = np.arange(1, k+1)
    plt.figure()
    
    plt.xlabel('input variable index')
    plt.ylabel('AES at max(N)')
    plt.plot(kk, aes_1, 'bx-', label = "Sobol")
    plt.plot(kk, aes_2, 'rx-', label = "LHS")
    plt.plot(kk, aes_3, 'gx-', label = "Random")
    plt.legend(loc='best')
    plt.show()
    
    plt.xlabel('input variable index')
    plt.ylabel('AEST at max(N)')
    plt.plot(kk, aest_1, 'bx-', label = "Sobol")
    plt.plot(kk, aest_2, 'rx-', label = "LHS")
    plt.plot(kk, aest_3, 'gx-', label = "Random")
    plt.legend(loc='best')
    plt.show()
    
def draw_maes_plot(N, maes_1, maes_2, maes_3, maest_1, maest_2, maest_3):
    
    plt.figure()
    
    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAES)')
    plt.xscale('log')
    plt.plot(N, maes_1, 'bx-', label = "Sobol")
    plt.plot(N, maes_2, 'rx-', label = "LHS")
    plt.plot(N, maes_3, 'gx-', label = "Random")
    plt.legend(loc='best')
    plt.show()
    

    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAEST)')
    plt.xscale('log')
    plt.plot(N, maest_1, 'bx-', label = "Sobol")
    plt.plot(N, maest_2, 'rx-', label = "LHS")
    plt.plot(N, maest_3, 'gx-', label = "Random")
    plt.legend(loc='best')
    plt.show()
    
    
    
    
    
    
    
    