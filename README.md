# ELEC5882M-project
This is a readMe file that documents all the source code files and images

1)  Gauss_newton5.m : This code performs an iterative non-linear Gauss-newton algebraic method for 3-D localisation of an initial guesstimate tag location near to its true original location within reasonable accuracy.
Reference paper used: Q. Luo, K. Yang, X. Yan, J. Li, C. Wang, and Z. Zhou, “An Improved Trilateration Positioning Algorithm with Anchor Node Combination and K-Means Clustering,” Sensors, vol. 22, no. 16, p. 6085, Aug. 2022, doi: 10.3390/s22166085.

2) Trilateration_3D.m: This code performs a simple linear algebraic algorithm to determine a position of tag based on ideal conditions using Cramer's rule.

3) Chan2d_tdoa.m : This code performs a non iterative linearised 2-D space TDoA algorithm to estimate an unknown tag's position.
   Reference papers used: Y. T. Chan and K. C. Ho, "A simple and efficient estimator for hyperbolic location," in IEEE Transactions on Signal Processing, vol. 42,       no. 8, pp. 1905-1915, Aug. 1994, doi: 10.1109/78.301830. 

    S. H. Olsen, “Radio Tracking of Open Range Sheep Methods for Radio Location in a Sub-GHz Base Station Network,” Department of Electronics and Telecommunications,     Jun. 2014. Accessed: Aug. 01, 2023.

4) ADC_plot.m : Plots ADC erros vs ideal characteristics
5) STM32L476RG_DW3000_code.txt : A C program file that documents all the register address and programming done for a simple transmit case.
