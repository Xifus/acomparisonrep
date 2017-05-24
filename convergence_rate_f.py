from math import log
import numpy as np
import matplotlib.pyplot as plt

def convergence_rate_N(N, Y, rule_name):
    if len(Y) != len(rule_name):
        print "you may miss something."
        return
    leng = len(Y)
    for i in range(leng):
        convergence_rate(N, Y[i], rule_name[i])


def convergence_rate(N, Y, rule_name):
    N_star = [float(i) for i in N]
    x_star = [log(i) for i in N_star]
    z = np.polyfit(x_star, Y, 1)
    print rule_name, "convergence rate = ", abs(z[0])
    print "-------------------------------------"
    return z

def check_trend_line(z, N, maest, name):
    N_star = [float(i) for i in N]
    x_star = [log(i) for i in N_star]
    p = np.poly1d(z)
    plt.figure()
    
    plt.xlabel('N(log scale)')
    plt.ylabel('log(MAES)')
    plt.xscale('log')
    plt.plot(N, maest, 'bx-', label = name)
    plt.plot(N, p(x_star), 'r--', label = "trendline")
    plt.legend(loc='best')
    plt.show()