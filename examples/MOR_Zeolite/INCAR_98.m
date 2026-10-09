% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [0.515760100441182    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [1.4126712325194368    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [0.4855335706674109     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [1.4911582450519045    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [0.3680896540844188     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [1.0781380478847749     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [0.610621955703535     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [0.41680864710501253    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [1.3436587514418261     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [0.8290030118654468    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [0.9797351017377145     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [1.326926680244555    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [1.8018033429475258     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [0.6015600297299362    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [2.1530994112932036     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [0.6973826300879615     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [0.7223318916574532    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [1.9466747657567618     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [1.0174849515326765    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [1.4712210865962745     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [1.3193926286529345    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.14905236706312266    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [0.8609081476955354     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [1.966594354265046     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [1.1725572850913224    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [0.9640217451424211     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [1.7958657540757579     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [0.6760988631458222    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.96220794728959     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [0.8839754503924461    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.0847556997821113     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [0.9714364836723178    -1.07   ]
CO_in <-> CO                                                                                  [0.1234589308532591    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [0.2320197139739933    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [1.264308549686842    -0.58   ]
% #1i for site i, #11 for site 1
% forwards and reverse reaction separated by <-> ,
% and use it as an identifier for Reaction Mechanism Equations 
% CO O2 CO2 for gas phase 
% CO#1i CO2#1i O#1i for adsorbed species
% Species(p) : relative pressure item of Species in rate equation
% Species(c) : relative coverage item of Species in rate equation
% instead of 2A by repeative writing A + A
% format : variable = value; Or matlab command;
% Q0 : initial coverage state
% Ea : energy barriers of each reaction
% G0 : Gibbs free energy of each reaction
% T  : temperature of the reaction
% P_Species : relative pressure of the Species
% C_Species : relative concentration of the Species
% Qi_Species : relative coverage of the Species at site i, v for vacancy
% Q_Species : relative coverage of the Species for only one site case 
% X_Species_INIT : the initial Species parameters X: P C Qi
% X_Species_FROZ : freeze the Species parameters X: P C Qi
% X_Species_EQUI : deal with Reaction i in equilibrium to solve X_Species
% matlab codes are acceptable in this file

Q0 = 1;
Ps = 1;       % total relative pressure of gas : 1 bar
Ng = 1; % total number of the gas molecular : 1 mol
Ns = 3.5E17/6.02E23; % number of active site : 1 mol
Vgt = 2.4/1000/3600; % rate of gas volume flow : m^3/s

PlotType = [1, 2, 3, 4];  
PlotMode = [1, 2, 3, 4];  

SkipMode = 2;
% npar = 4;
% Istart = 0;
% MaxTime = 20;
% MaxOdeTime = 1000;
tspan = [0 1e6];
CalcDRC = 1;
P_CH2CO_in_INIT = 0.01;
P_CO_INIT = 0;
P_CO_in_INIT = 0;
P_H2_INIT = 1;
P_CH2CH2_INIT = 0;
P_CH2CH2_in_INIT = 0;
P_CH2CHCH3_INIT = 0;
P_CH2CHCH3_in_INIT = 0;
Q1_v_INIT = 1;
T = 650;
SimMode = [0];
ProfMode = [4];
SimValue = 0;
FigMode = 3
