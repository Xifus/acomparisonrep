import a_comparison as acn
import sampling_method as sm
import numpy as np

def create_a(kk):
#choose dimension
    a1 = sm.create_coefficient_a(kk, "A1-1")
    a2 = sm.create_coefficient_a(kk, "A1-2")
    a3 = sm.create_coefficient_a(kk, "A2")
    a4 = sm.create_coefficient_a(kk, "B")
    a5 = sm.create_coefficient_a(kk)
    return a1, a2, a3, a4, a5

def sum_si(k, a1, a2, a3, a4, a5):
    s_ana_1, st_ana_1 = acn.analyticalSAvalues(k, a1)
    s_ana_2, st_ana_2 = acn.analyticalSAvalues(k, a2)
    s_ana_3, st_ana_3 = acn.analyticalSAvalues(k, a3)
    s_ana_4, st_ana_4 = acn.analyticalSAvalues(k, a4)
    s_ana_5, st_ana_5 = acn.analyticalSAvalues(k, a5)
    s_sum_1 = np.sum(s_ana_1)
    s_1_st_1_1 = s_ana_1[0]/st_ana_1[0]
    s_k_st_k_1 = s_ana_1[k-1]/st_ana_1[k-1]
    s_sum_2 = np.sum(s_ana_2)
    s_1_st_1_2 = s_ana_2[0]/st_ana_2[0]
    s_k_st_k_2 = s_ana_2[k-1]/st_ana_2[k-1]
    s_sum_3 = np.sum(s_ana_3)
    s_1_st_1_3 = s_ana_3[0]/st_ana_3[0]
    s_k_st_k_3 = s_ana_3[k-1]/st_ana_3[k-1]
    s_sum_4 = np.sum(s_ana_4)
    s_1_st_1_4 = s_ana_4[0]/st_ana_4[0]
    s_k_st_k_4 = s_ana_4[k-1]/st_ana_4[k-1]
    s_sum_5 = np.sum(s_ana_5)
    s_1_st_1_5 = s_ana_5[0]/st_ana_5[0]
    s_k_st_k_5 = s_ana_5[k-1]/st_ana_5[k-1]
    print "k = ", k
    print "Type A1-1"
    print "sum(S_i) = ", s_sum_1
    print "S_1/S_Tot_1 = ", s_1_st_1_1
    print "S_k/S_Tot_k = ", s_k_st_k_1
    print "---------------------------"
    print "Type A1-2"
    print "sum(S_i) = ", s_sum_2
    print "S_1/S_Tot_1 = ", s_1_st_1_2
    print "S_k/S_Tot_k = ", s_k_st_k_2
    print "---------------------------"
    print "Type A2"
    print "sum(S_i) = ", s_sum_3
    print "S_1/S_Tot_1 = ", s_1_st_1_3
    print "S_k/S_Tot_k = ", s_k_st_k_3
    print "---------------------------"
    print "Type B"
    print "sum(S_i) = ", s_sum_4
    print "S_1/S_Tot_1 = ", s_1_st_1_4
    print "S_k/S_Tot_k = ", s_k_st_k_4
    print "---------------------------"
    print "Type C"
    print "sum(S_i) = ", s_sum_5
    print "S_1/S_Tot_1 = ", s_1_st_1_5
    print "S_k/S_Tot_k = ", s_k_st_k_5
    
    
    
    
    
    
    
    
    
    
    
    
    
    