% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m 
% examples :
 
CO + #1 <-> CO#1                                       [  0.626   0.7300 ] 
CO + #2 <-> CO#2                                       [  0.593   0.9100 ] 
H2 + #1 <-> H2#1                                       [  0.549   0.9500 ] 
H2#1 + #2 <-> H#1 + H#2                                [  0.430  -0.7100 ] 
H#2 + #1 <-> H#1 + #2                                  [     0   -0.6400 ] 
CO#1 + H#1 <-> CHO#1 + #1                              [  0.620   0.5900 ] 
CO#1 + H#1 <-> COH#1 + #1                              [  2.220   1.9400 ] 
CHO#1 + H#1 <-> CH2O#1 + #1                            [  0.290  -0.3200 ] 
CHO#1 + H#1 <-> CHOH#1 + #1                            [  1.42    0.4200 ] 
CH2O#1 + H#1 <-> CH3O#1 + #1                           [  0.850  -0.0100 ] 
CH2O#1 + H#1 <-> CH2OH#1 + #1                          [  2.71    1.5800 ] 
COH#1 + H#1 <-> CHOH#1 + #1                            [  0.55   -0.9400 ] 
CHOH#1 + H#1 <-> CH2OH#1 + #1                          [  1.66    0.8400 ] 
CH2OH#1 + H#1 <-> CH3OH#1 + #1                         [  1.04    0.3100 ] 
CH3O#1 + H#1 <-> CH3OH#1 + #1                          [  1.91    1.9100 ] 
CH2O#1 <-> CH2O + #1                                   [  0.848   0.2200 ] 
CH3OH#1 <-> CH3OH + #1                                 [ -0.79   -1.4200 ] 
CHO#1 + #1 <-> CH#1 + O#1                              [  0.22   -0.9300 ] 
CH2O#1 + #1 <-> CH2#1 + O#1                            [  1.08   -1.070  ] 
CH2O#1 + #2 <-> CH2#2 + O#1                            [  2.68    0.2900 ] 
CHOH#1 + #1 <-> CH#1 + OH#1                            [  0.36   -1.1200 ] 
CH2OH#1 + #1 <-> CH2#1 + OH#1                          [  1.22   -2.4200 ] 
CH2OH#1 + #2 <-> CH2#2 + OH#1                          [  0.30   -1.0600 ] 
CH3O#1 + #1 <-> CH3#1 + O#1                            [  1.83   -0.5400 ] 
CH3O#1 + #2 <-> CH3#2 + O#1                            [  1.77   -0.5300 ] 
CH3OH#1 + #1 <-> CH3#1 + OH#1                          [  0.53   -2.2200 ] 
CH3OH#1 + #2 <-> CH3#2 + OH#1                          [  0.74   -2.2100 ] 
CH#1 + H#1 <-> CH2#1 + #1                              [  0.94   -0.460  ] 
CH#1 + OH#1 <-> CH2#1 + O#1                            [  0.53   -0.6900 ] 
CH2#1 + H#1 <-> CH3#1 + #1                             [  1.57    0.51   ] 
CH2#2 + H#1 <-> CH3#2 + #1                             [  0.55   -0.84   ] 
CH2#2 + OH#1 <->CH3#2 + O#1                            [  0.74   -1.07   ] 
CH2#1 + OH#1 <->CH3#1 + O#1                            [  0.87    0.28   ] 
CH3#2 + H#2 <-> CH4 + #2 + #2                          [  1.15   -1.710  ] 
CH3#2 + OH#1 <-> CH4 + #2 + O#1                        [  0.75   -1.30   ] 
CH3#2 + H#1 <-> CH4 + #1 + #2                          [  1.27   -1.07   ] 
CH3#1 + H#2 <-> CH4 + #1 + #2                          [  1.11   -1.70   ] 
CH3#1 + H#1 <-> CH4 + #1+ #1                           [  1.68   -1.06   ] 
CH3#1 + OH#1 <-> CH4 + #1+ O#1                         [  0.79   -1.29   ] 
CHO#1 + CO#1 <-> CH#1 + CO2#1                          [  1.14   -0.6300 ] 
CH2O#1 + CO#1 <-> CH2#1 + CO2#1                        [  1.930  -0.7800 ] 
CH3O#1 + CO#1 <-> CH3#1 + CO2#1                        [  2.59   -0.2500 ] 
CHO#1 + CO#2 <-> CH#1 + CO2#2                          [  2.35    0.3100 ] 
CH2O#1 + CO#2 <-> CH2#1 + CO2#2                        [  2.55    0.1700 ] 
CH3O#1 + CO#2 <-> CH3#1 + CO2#2                        [  3.32    0.7000 ] 
CH2#1 + CO#2 <-> CH2CO#1 + #2                          [  0.50    0.2800 ] 
CH2#1 + CO#1 <-> CH2CO#1 + #1                          [  0.88    0.4600 ] 
CH2#2 + CO#2 <-> CH2CO#2 +#2                           [  0.20   -0.6300 ] 
CH3#1 + CO#2 <-> CH3CO#1 + #2                          [  0.62    0.5000 ] 
CH3#1 + CO#1 <-> CH3CO#1 + #1                          [  0.97    0.6700 ] 
CH3#2 + CO#2 <-> CH3CO#2 +#2                           [  0.53   -0.1000 ] 
CH3CO#1 + #1 <-> CH2CO#1 + H#1                         [  0.21   -0.7300 ] 
CH3CO#1 + O#1 <-> CH2CO#1 + OH#1                       [  0.00   -0.4000 ] 
CH3CO#2 + #1 <-> CH2CO#2 + H#1                         [  1.89    0.3200 ] 
CH3CO#2 + O#1 <-> CH2CO#2 + OH#1                       [  2.31    0.5400 ] 
CH2CO#1 <-> CH2CO + #1                                 [ -0.093  -0.730  ] 
CH2CO#2 <-> CH2CO + #2                                 [ -0.575  -1.180  ] 
CO#2 + O#1 <-> CO2#1 + #2                              [  0.15    0.1200 ] 
CO#2 + OH#1 <-> COOH#2 +#1                             [  1.7800  1.4800 ] 
COOH#2 + O#1 <-> CO2#2 + OH#1                          [  0.280  -0.2400 ] 
COOH#2 + #1 <-> CO2#2 + H#1                            [  3.73   -0.4700 ] 
H#1 + OH#1 <-> H2O#1 + #1                              [  1.89    1.8900 ] 
OH#1 + OH#1 <-> H2O#1 + O#1                            [  1.600   1.6400 ] 
H#2 + OH#1 <-> H2O#1 + #2                              [  1.96    1.2300 ] 
CO2#1 <-> CO2 + #1                                     [  0.329  -0.3100 ] 
CO2#2 <-> CO2 + #2                                     [ -0.824  -1.4300 ] 
H2O#1 <-> H2O + #1                                     [ -0.227  -0.8400 ] 

Q0 = [1 1];
Ps = 25;       % total relative pressure of gas : 1 bar
Ng = 1; % total number of the gas molecular : 1 mol
Ns = 6E19/6.02E23; % number of active site : 1 mol
Vgt = 0.720/1000/3600; % rate of gas volume flow : m^3/s

Istart = 0;
MaxTime = 20;
CalcDRC = 1;
%BarrierMode = 2;
%As = [15.10  26.40];
%ThermoMode = 6;
%Mr_CH4 = 0.5;
%Mr_CO = 0.5;
%Mr_H2O = 0.5;
%Mr_H2 = 0.5;
%Mr_CO2 = 0.5;
SimMode = [0];
ProfMode = [4];
SimValue = 0;
P_H2_INIT = 17.857;
P_CO_INIT = 7.143;
P_CO2_INIT = 0;
P_H2O_INIT = 0;
P_CH3OH_INIT = 0
P_CH2CO_INIT = 0
P_CH2O_INIT = 0
P_CH4_INIT = 0
P_C3H6_INIT = 0
P_CH4_a_INIT = 0
% Q1_CO_EQUI = 1;
% Q2_CO_EQUI = 2;
% Q1_H2_EQUI = 3;
Q1_v_INIT = 1;
Q2_v_INIT = 1;
T = 673;
Thermo_CH2CO = [0 0 0 2.08 0.84 0];



