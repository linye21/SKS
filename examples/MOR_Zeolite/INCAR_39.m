% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1                                                                   [0.5818804806280583    -0.40   ]
CH2CO_H#1 <-> CH3CO#1                                                                         [1.0155490343739826    -0.31   ]
CH3CO#1 <-> CH3#1 + CO_in                                                                     [0.4342716307840348     0.42   ]
CH3#1 + CH2CO_in <-> CH3CH2CO#1                                                               [0.608084534514466    -0.64   ]
CH3CH2CO#1 <-> CH3CH2#1 + CO_in                                                               [0.5809066624059567     0.06   ]
CH3CH2#1 <-> CH2CH2_in + #1                                                                   [1.957594530409972     0.69   ]
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1                                                         [0.7714336146522056     0.00   ]
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in                                                         [1.380326146687191    -0.07   ]
CH3CH2CH2#1 <-> CH2CHCH3_in + #1                                                              [0.348513072316092     0.33   ]
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1                                                           [1.082916809690225    -0.15   ]
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in                                                           [1.3119579042401595     0.21   ]
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1                                                           [0.22024518080621192    -1.18   ]
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1                                                       [1.4670637664672033     0.43   ]
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in                                                       [0.4941990657667576    -0.72   ]
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1                                                            [2.099557024885474     0.88   ]
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1                                                 [1.8918467769040952     0.43   ]
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in                                                 [1.457138162031404    -0.74   ]
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1                                                           [2.3553820040911395     1.39   ]
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1                                                       [0.022765578233819035    -0.73   ]
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1                                                             [1.3554506973107916     0.17   ]
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in                                                            [1.0907860480302416    -1.04   ]
CH3COCH2CO#1 <-> CH3COCHCO_H#1                                                                [0.5209811249889799    -0.38   ]
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1                                                    [1.7304432548877897     0.63   ]
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1                                                       [1.2685609098842097     0.60   ]
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in                                                       [0.6174659649261829    -1.15   ]
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1                                                            [1.1796340996236443     0.74   ]
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1                                                 [1.652623900274285     0.67   ]
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in                                                 [1.063181208617205    -0.96   ]
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1                                                           [2.823194819383078     1.66   ]
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1                                                       [1.135185652163663    -1.06   ]
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1                                                             [1.8805046306978985     0.59   ]
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in                                                            [1.2892915561170608    -1.07   ]
CO_in <-> CO                                                                                  [0.2381601914926563    -0.50   ]
CH2CH2_in <-> CH2CH2                                                                          [0.9030665143673002    -0.56   ]
CH2CHCH3_in <-> CH2CHCH3                                                                      [0.2219089909299714    -0.58   ]
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
