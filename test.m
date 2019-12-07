%occupy=ones(n_e,1);   %denote whether the electron exists
%% Read table for elastic scattering 
% Find total cross section: 
% For energy > 50 eV
T_el = readtable("TotalElasticCS_Nitrogen.csv");
t_el = table2array(T_el);
% For 1 eV < energy < 50 eV
T_sel = readtable("TotalElasticCS_smalleV_Nitrogen.csv");
t_sel = table2array(T_sel);

% Find differential cross section: 
T_del = readtable("ElasticCS_Nitrogen.csv");
t_del = table2array(T_del);

%% Read table for total cross section of impact ionization
T_ion = readtable("TotalIonCS_Nitrogen.csv");
t_ion = table2array(T_ion);

%% Read table of parameters of impact excitation
% Fitting parameters for Nitrogen:
T_exc = readtable("ParameterExcCS_Nitrogen.csv");
t_exc = table2array(T_exc);

%% Read table of parameters of impact excitation
% Fitting parameters for Nitrogen:
T_exc = readtable("ParameterExcCS_Nitrogen.csv");
t_exc = table2array(T_exc);

%%  
n =10;
T = 0.01*rand(n,1)
%% Pre-set data
% constant
c = 3*10^8;
m = 9.109 * 10^(-31); % electron rest mass
E_e = m * c^2 / (1.6022 * 10 ^ (-19)); 

% Fitting parameter of Nitrogen
I = 15.6/1000; % Ionization threshold, unit : eV
T_a = 1000; % unit : eV
T_b = 2*I;
t_b = I;
T_s = 4.17; % unit : eV
t_s = 13.8; % unit : eV



% Calculated value
%% the energy higher than the threshold
GoToExcitation = find(T<I);
n = length(T);



tic
[sig_elastic,sig_ionization,sig_excitation]=find_CS0( t_el, t_sel, t_ion, t_exc, T );
toc
