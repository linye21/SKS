% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [0.7523230101885386    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [1.0650300584936807    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [1.3860470052800975     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [0.5381530058196443    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [0.34835414710664575     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [1.1625638259001618     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [1.100308025052124     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [1.3268122348182416    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [1.400457826004462     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [0.6360235698719721    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [1.6253863536030475     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [0.3845279536110614    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [0.6395008823996176     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [1.0701262269340437    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [1.3646955086844634     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [1.7107628095091416     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [0.020896309825436235    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [2.119227624604271     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [0.5193653499554979    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [0.42737782817335357     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [0.16739740618282334    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.8883523637528388    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [1.5031697222511333     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [0.908161615145537     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [0.20604572006235727    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [2.2122953227901125     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [0.9545058602002192     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [1.1390924368762074    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.9677337429010935     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [0.4040409346311264    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.5021580106547034     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [1.1513690798034666    -1.07   ]
CO_in <-> CO                                                                                  [1.0214814223733308    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [0.5880334206311484    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [0.4806903576150324    -0.58   ]
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
