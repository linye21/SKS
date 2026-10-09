% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [0.46424069688534875    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [0.3286090207005994    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [1.5608759532076633     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [0.969228272838807    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [0.7434491858770489     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [0.7841615122526313     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [0.708414104129007     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [0.7646724747825604    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [0.889806542798208     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [0.04606425542503808    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [0.8539854347307205     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [0.7690940195212617    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [1.4644394183095641     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [0.1297443158154078    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [2.1833204836081106     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [1.8393654201039589     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [0.9923734903020102    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [2.553175630574101     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [0.9492958944585521    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [0.9603563905961795     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [0.9368743344266698    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.2317502616990369    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [1.9357554793843272     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [1.9336568860124315     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [0.16619235756726994    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [1.2843979740553595     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [1.7662809675855131     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [0.5308940089949513    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.093726866147831     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [0.061703450121491824    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.9005547900006179     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [0.9464837162761746    -1.07   ]
CO_in <-> CO                                                                                  [1.4914329178598293    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [1.4443992347253634    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [0.2866446266994636    -0.58   ]
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
