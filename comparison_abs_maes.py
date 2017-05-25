import matplotlib.pyplot as plt
import numpy as np

def draw_comparison_abs_maes(N, maes_d_1, maes_d_2, name1, name2):
    
    
    A = 0.85 * max(max(maes_d_1), max(maes_d_2))
    if A <= 0:
        A = 0.05
    B = 0.85 * min(min(maes_d_1), min(maes_d_2))
    label_y = 'log(MAES, '+ name1+ ') - log(MAES,'+ name2+ ')'
    index = np.arange(len(N))
    width = 0.3
    text_pos = len(N)/4 + width
 
    plt.bar(index + 0.5*width, maes_d_1, width,
                 label='k = 10')
 
    plt.bar(index + 1.5*width, maes_d_2, width,
                 label='k = 19')
 
    plt.xlabel('Number of model runs')
    plt.ylabel(label_y)
    plt.xticks(index + width, N)
    plt.text(text_pos, A,'LSS advantage')
    plt.text(text_pos, B,'Sobol advantage')
    plt.legend()
    plt.show()