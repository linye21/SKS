% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [0.07900904297419416    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [0.4304145942187356    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [1.4061996909561503     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [1.022718253307718    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [0.6812791074388127     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [1.940434679036485     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [1.0289670498648986     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [0.013420255003293392    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [0.4114257007942158     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [0.16670679991437964    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [1.4119784121361139     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [1.1899300194424711    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [1.7127664659180144     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [0.2161106934441851    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [1.718790973360747     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [0.7652906019717163     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [0.8802104774693786    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [2.879677149873933     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [0.3790059017372155    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [0.21762634932730596     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [1.3461849868082811    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.22284683754967743    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [1.2364356597305888     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [1.446834314190863     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [0.5639982533777135    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [1.751657456965986     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [1.8938232965328865     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [0.6136301557008913    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.4869063916986365     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [1.1535893687164784    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.985441041997591     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [1.0347628079796598    -1.07   ]
CO_in <-> CO                                                                                  [1.2665208817461793    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [0.09222071572522653    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [0.16386133572869976    -0.58   ]
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
