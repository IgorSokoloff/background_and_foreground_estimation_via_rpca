# background_and_foreground_estimation_via_rpca

Here I store the Geometric data analysis methods course project which I did on the 9th term at MIPT.

The goal of the project is to implement a  background/foreground estimation of the given set of video frames using robust PCA.
The algorithms were taken from the following papers:

1. ["A Nonconvex Projection Method for Robust PCA"](https://arxiv.org/abs/1805.07962)\
(Aritra Dutta, Filip Hanzely, Peter Richtarik)\
In the source code I will call it "rpca-alt"(Alternating projection method for RPCA) - The Algorithm 2 on the page 5 of the paper

2. ["Robust PCA using Matrix Factorization for Background/Foreground Separation"](https://www.researchgate.net/publication/324051198_Robust_PCA_using_Matrix_Factorization_for_BackgroundForeground_Separation)\
(Shuqin Wang, Yongli Wang, Yongyong Chen, Peng Pan, Zhipeng Sun, and Guoping He)\
Later I will call it "rpca-ALF" (RPCA with Augmented Lagrangian Function)
