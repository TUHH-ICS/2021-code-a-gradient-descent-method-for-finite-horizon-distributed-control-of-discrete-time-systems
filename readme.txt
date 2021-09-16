This folder contains the files to reproduce the results presented in

S. Heinke and H. Werner, "A Gradient Descent Method for Finite Horizon 
Distributed Control of Discrete Time Systems", Conference on Decision 
and Control, 2021.

The figures can be reproduced using the file 'generate_plots.m'. 
The data for Figures 1 and 2 is generated using 'main_sf.m', the data 
for Figure 3, 4 and 5 is generated using 'main_sf_comp.m', 'main_of.m' 
and 'main_of_comp.m', respectively.

In order to run the code, the following Matlab packages are required:
- yalmip (https://yalmip.github.io/download/)
- CVX (http://cvxr.com/cvx/download/)

The code has been tested using:
- Matlab R2020a
- yalmip version 21-Nov-2017
- CVX version 2.1

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.