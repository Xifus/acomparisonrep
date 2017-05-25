These codes are trying to reproduce the results from Stefano Tarantola, William Becker, and Dirk Zeitz's paper "A comparison of two sampling methods for global sensitivity analysis".
In S. Tarantola et al, they used Sobol and LSS(Latin Supercube Sampling). For here, I used Sobol, Latin Hypercube Sampling, and Random by using Chaospy library (https://github.com/jonathf/chaospy).
The main code files are:
  a_comparison.py: function to calculate test function; analytical values for corresponding test function; function for estimated variance; function for estimated first order effects and total effects (in two ways: Jansen 1999 e_sti_1; Sobol' 2007 e_sti_2); function for measuring convergence AES, MAES, AEST, MAEST. And also there are two tool functions for shifting sample matrix (less_eq_integer_matrix and shift).
  calculation_plot.py: This is the main calculation code. It stores value of R replicated times of estimated si and sti for multiple N into a big matrix. After that, the function passes the matrix to calculate AES, MAES, AEST, MAEST. Finally, I use the value of AES, MAEST, AEST, MAEST to draw convergence and i_wise error plot for each case.
  compare_sti.py: This is a modified calculation_plot.py code for comparing the convergence rate and efficience of Jansen 1999 method and Sobol' 2007 method, according to A. Saltelli et al. "Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index". I keep sampling method, sampling matrix all the same and just change the total effects calculation method.
  comparison_abs_maes.py: This is a function to graph Fig.7 comparison of absolute MAES values for different test case function. Because Chaospy's Sobol sampling is only up to 40, so I only compared k = 10 and k = 19.
  convergence_rate_f.py: This contains functions to use trend line to estimate convergence rates of AES, MAES, AEST, MAEST. Then I draw the trend line and the estimating plot to check the reliability of the trend line.
  sampling_method.py: I put the sampling method, sampling matrix, and coefficient (different text case A1-1, B, C etc) creating functions in this module. I can easily add more rules in the future.
  test.py: This is a simple unittest code for testing a_comparison.py and sampling_method.py. It needs update for new modules.
  test_ana.py: This module tests if the analytical values are the same as indicated on the paper. The result is in Cehck Analytical Values.ipynb.
  
I use ipython notebook to visualize the whole process and final results.
A comparison of sampling methods "test case" k = 10 notebooks: 
  for checking the result for one k (dimension). The expected result is Sobol > LHS > Random. From what I and my supervisor learned from Chaospy library, LHS sampling has a similar behavior as Random, so the results of LHS and Random are pretty close at some cases.
A comparsion of sampling methods "test case" k = 10, k = 19 notebooks:
  This notebook contains the code to calculate k = 10 case and repeat for k = 19 case. In addition, this also contains convergence rates and trendline checking plots. At the very end, I did the comparison of absolute MAES values inbetween the three sampling methods. There are currently some typo on the bar graphs. I will fix them as soon as possible in a later date.
Compare Jansen 1999 and Sobol' 2007 "test case" k = 10 notebooks:
  This uses convergence_rate_f.py and compare_sti.py to compare the convergence rate and efficience of two different sti method Jansen 1999 and Sobol' 2007. The expected result is that Jansen 1999 is better than Sobol' 2007 according to A. Saltelli et al.
