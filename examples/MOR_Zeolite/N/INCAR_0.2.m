% '%' and '!' for comments
% From Reaction Mechanism Equations To Kinetic Equations
% Create Reaction Mechanism Equations in file ReactionEquation.m▫
% examples :

CH2CO_in + #1 <-> CH2CO_H#1              [ 0.7500      -0.4000    ] 
CH2CO_H#1 <-> CH3CO#1                    [ 0.7500      -0.3100    ] 
CH3CO#1 <-> CH3#1 + CO_in                [ 1.1700      0.4200     ] 
CH3#1 + CH2CO_in <-> CH3CH2CO#1          [ 0.7500      -0.6400    ] 
CH3CH2CO#1 <-> CH3CH2#1 + CO_in          [ 0.8100      0.0600     ] 
CH3CH2#1 <-> CH2CH2_in + #1              [ 1.4400      0.6900     ] 
CH3CH2#1 + CH2CO_in <-> CH3CH2CH2CO#1    [ 0.7500      0.0000     ] 
CH3CH2CH2CO#1 <-> CH3CH2CH2#1 + CO_in    [ 0.7500      -0.0700    ] 
CH3CH2CH2#1 <-> CH2CHCH3_in + #1         [ 1.0800      0.3300     ] 
CH3CO#1 + CH2CO_in <-> CH3COCH2CO#1      [ 0.7500      -0.1500    ] 
CH3COCH2CO#1 <-> CH3COCH2#1 + CO_in      [ 0.9600      0.2100     ] 
CH3COCH2#1 + CH2CO_in <-> FMR_CH3#1      [ 0.7500      -1.1800    ] 
FMR_CH3#1 + CH2CO_in <-> FMR_CH3CH2CO#1  [ 0.5600      0.4300     ] 
FMR_CH3CH2CO#1 <-> FMR_CH2CH3#1 + CO_in  [ 0.7500      -0.7200    ] 
FMR_CH2CH3#1 <-> CH2CH2_in + FMR#1       [ 2.0000      0.8800     ] 
FMR_CH2CH3#1 + CH2CO_in <-> FMR_CH2CH3CH2CO#1 [ 0.7900      0.4300     ] 
FMR_CH2CH3CH2CO#1 <-> FMR_CH2CH2CH3#1 + CO_in [ 1.4300      -0.7400    ] 
FMR_CH2CH2CH3#1 <-> FMR_HCH2CHCH3#1      [ 1.3900      1.3900     ] 
FMR_HCH2CHCH3#1 <-> CH2CHCH3_in + FMR#1  [ 0.0100      -0.7300    ] 
FMR#1 + CH2CO_in <-> FMR_HCH2CO#1        [ 0.9200      0.1700     ] 
FMR_HCH2CO#1 <-> FMR_CH3#1 + CO_in       [ 0.7500      -1.0400    ] 
CH3COCH2CO#1 <-> CH3COCHCO_H#1           [ 0.7500      -0.3800    ] 
CH3COCHCO_H#1 + CH3CO#1 <-> SMR_CH3#1 + #1 [ 1.3800      0.6300     ] 
SMR_CH3#1 + CH2CO_in <-> SMR_CH3CH2CO#1  [ 1.3500      0.6000     ] 
SMR_CH3CH2CO#1 <-> SMR_CH2CH3#1 + CO_in  [ 0.7500      -1.1500    ] 
SMR_CH2CH3#1 <-> CH2CH2_in + SMR#1       [ 1.4900      0.7400     ] 
SMR_CH2CH3#1 + CH2CO_in <-> SMR_CH2CH3CH2CO#1 [ 1.4200      0.6700     ] 
SMR_CH2CH3CH2CO#1 <-> SMR_CH2CH2CH3#1 + CO_in [ 0.7500      -0.9600    ] 
SMR_CH2CH2CH3#1 <-> SMR_HCH2CHCH3#1      [ 2.4100      1.6600     ] 
SMR_HCH2CHCH3#1 <-> CH2CHCH3_in + SMR#1  [ 0.7500      -1.0600    ] 
SMR#1 + CH2CO_in <-> SMR_HCH2CO#1        [ 1.3400      0.5900     ] 
SMR_HCH2CO#1 <-> SMR_CH3#1 + CO_in       [ 0.7500      -1.0700    ] 
CO_in <-> CO                             [ 0.0000      -0.5000    ] 
CH2CH2_in <-> CH2CH2                     [ 0.7500      -0.5600    ] 
CH2CHCH3_in <-> CH2CHCH3                 [ 0.7500      -0.5800    ] 


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
