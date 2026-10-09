% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [1.0007383489623964    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [0.4162779432652064    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [0.5844171915624643     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [1.2504955044724717    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [1.3434871725168178     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [0.7062843862865624     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [0.8858217621225595     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [1.2121855512210151    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [0.3638583311151545     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [1.2850156426096357    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [1.1906312011312492     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [1.4181256114829528    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [1.1579243573097264     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [0.8962025831624061    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [1.5139284896795195     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [1.0373024516478835     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [0.5184269598653416    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [2.090704025985853     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [1.1368100776252472    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [1.2422739208644875     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [0.6164781735198647    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.6076440589814626    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [1.056501612370702     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [1.5359540624496462     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [1.13728202458889    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [1.6759087601195854     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [0.9773814229630858     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [1.469854284699383    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.529531969318141     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [1.4909064072063591    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.7771909396839503     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [0.22101722522079154    -1.07   ]
CO_in <-> CO                                                                                  [0.9260811615604021    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [0.9916818307580181    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [0.6204124083668136    -0.58   ]
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
