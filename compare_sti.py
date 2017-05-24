import a_comparison as acn
import numpy as np
import matplotlib.pyplot as plt
import sampling_method as sm

def n_calculate_compare_sti(k, N, a, R, rl, rule = "V-function"):
    leng = len(N)
    #er_si = np.zeros((R, k))
    #er_sti = np.zeros((R, k))
    er_si_N = np.zeros(((leng * R), k))
    er_sti_N_1 = np.zeros(((leng * R), k))
    er_sti_N_2 = np.zeros(((leng * R), k))
    for i in range(leng):
        c = i*R
        er_si, er_sti_1, er_sti_2 = calculate_compare_sti(k, N[i], a, R, rl, rule)
        er_si_N[c:c+R, :] = er_si
        er_sti_N_1[c:c+R, :] = er_sti_1
        er_sti_N_2[c:c+R, :] = er_sti_2
    return er_si_N, er_sti_N_1, er_sti_N_2

def calculate_compare_sti(k, N, a, R, rl = "S", rule = "V-function"):
    #sampling by chaospy method. rule = sampling method, R = replicated times
    er_si = np.zeros((R, k))
    er_sti_1 = np.zeros((R, k))
    er_sti_2 = np.zeros((R, k))
    
    
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
        estimated_sti_2 = np.zeros((k))
        estimated_si = acn.e_si(k, X_shift, a, N, estimated_var, rule)
        estimated_sti_1 = acn.e_sti_1(k, X_shift, a, N, estimated_var, rule)
        estimated_sti_2 = acn.e_sti_2(k, X_shift, a, N, estimated_var, rule)
        
        er_si[i, :] = estimated_si
        er_sti_1[i, :] = estimated_sti_1
        er_sti_2[i, :] = estimated_sti_2
               
    return er_si, er_sti_1, er_sti_2

def aes_aest_N(N, k, R, a, er_si_N, er_sti_N_1, er_sti_N_2):
    s_ana, st_ana = acn.analyticalSAvalues(k, a)
    leng = len(N)
    aes_N = np.zeros((leng, k))
    aest_N_1 = np.zeros((leng, k))
    aest_N_2 = np.zeros((leng, k))
    for i in range(leng):
        y = i*R
        z = (i+1)*R
        e_si = er_si_N[y:z,:]
        e_sti_1 = er_sti_N_1[y:z,:]
        e_sti_2 = er_sti_N_2[y:z,:]
        aes, aest_1, aest_2 = aes_aest(R, e_si, e_sti_1, e_sti_2, s_ana, st_ana, k)
        
        aes_N[i, :] = aes
        aest_N_1[i, :]= aest_1
        aest_N_2[i, :]= aest_2        
    return aes_N, aest_N_1, aest_N_2

def maes_maest_N(N, k, aes_N, aest_N_1, aest_N_2):
    leng = len(N)
    maes = np.zeros(leng)
    maest_1 = np.zeros(leng)
    maest_2 = np.zeros(leng)
    for i in range(leng):
        #aes = aes_N[i, :].flatten()
        #aest = aest_N[i, :].flatten()
        aes = aes_N[i, :]
        aest_1 = aest_N_1[i, :]
        aest_2 = aest_N_2[i, :]
        maes[i] = acn.m_a_e_s_st(k, aes)
        maest_1[i] = acn.m_a_e_s_st(k, aest_1)
        maest_2[i] = acn.m_a_e_s_st(k, aest_2)
    return maes, maest_1, maest_2
        
        

def aes_aest(R, e_si, e_sti_1, e_sti_2, s_ana, st_ana, k):
    aes = acn.a_e_s_st(R, e_si, s_ana, k)
    aest_1 = acn.a_e_s_st(R, e_sti_1, st_ana, k)
    aest_2 = acn.a_e_s_st(R, e_sti_2, st_ana, k)
    return aes, aest_1, aest_2

def draw_aest_compare_plot(k, a, R, e_si, e_sti_1, e_sti_2, rule):
    s_ana, st_ana = acn.analyticalSAvalues(k, a)
    aes, aest_1, aest_2 = aes_aest(R, e_si, e_sti_1, e_sti_2, s_ana, st_ana, k)
    
    kk = np.arange(1, k+1)
    plt.figure()
    
    plt.xlabel('input variable index')
    plt.ylabel('AEST at max(N)')
    plt.title(rule)
    plt.plot(kk, aest_1, 'bx-', label = "Jansen 1999")
    plt.plot(kk, aest_2, 'rx-', label = "Sobol' 2007")
    plt.legend(loc='best')
    plt.show()
    
def draw_maest_compare_plot(N, maest_1, maest_2, rule):
    
    plt.figure()    
    plt.title(rule)
    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAEST)')
    plt.xscale('log')
    plt.plot(N, maest_1, 'bx-', label = "Jansen 1999")
    plt.plot(N, maest_2, 'rx-', label = "Sobol' 2007")
    plt.legend(loc='best')
    plt.show()
    
def draw_maest_plot(N, maes_1, maes_2, maes_3, maest_1, maest_2, maest_3):
    
    plt.figure()
    
    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAEST)')
    plt.xscale('log')
    plt.plot(N, maes_1, 'bx-', label = "Sobol")
    plt.plot(N, maes_2, 'rx-', label = "LHS")
    plt.plot(N, maes_3, 'gx-', label = "Random")
    plt.title("Jansen 1999")
    plt.legend(loc='best')
    plt.show()
    

    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAEST)')
    plt.xscale('log')
    plt.plot(N, maest_1, 'bx-', label = "Sobol")
    plt.plot(N, maest_2, 'rx-', label = "LHS")
    plt.plot(N, maest_3, 'gx-', label = "Random")
    plt.title("Sobol' 2007")
    plt.legend(loc='best')
    plt.show()    
    
    
    
    
    
    