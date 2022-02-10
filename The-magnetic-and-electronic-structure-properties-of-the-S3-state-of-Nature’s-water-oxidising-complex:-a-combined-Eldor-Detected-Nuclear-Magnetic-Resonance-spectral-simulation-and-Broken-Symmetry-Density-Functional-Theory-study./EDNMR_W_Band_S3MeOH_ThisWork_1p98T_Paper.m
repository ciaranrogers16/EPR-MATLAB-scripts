%% 1D EDNMR Simulation
%  Simulating 1D EDNMR data using horseradish S3 MeOH state 
%  Adapted from Nino Wili's code:
%  Jeschke, G. et. al., Phys. Chem. Chem. Phys., 2019, 21, 11676-11688
%  DOI: 10.1039/c9cp01760g
%
%  olivia.hardwick@student.manchester.ac.uk
%  ciaran.rogers@manchester.ac.uk
% =========================================================================
close all
clear all

% 1D EDNMR Simulation
MWFQ=94.017; % define microwave frequency in GHz (W Band)

%% Define EDNMR Experimental Parameters (Ref 24 Olivia's MChem report)
% EDNMR simulation at 1.98T
Exp.mwFreq = MWFQ;
Exp.Field = 1980; % Field position mT 
Exp.ExciteWidth = 100; %Excitation Bandwidth
Exp.nPoints = 1024; %Number of Points
Exp.Range = [0 300]; %Range
nu_eld = 4.7e7; % rad/s*conversion_mhz
Exp.nu1 = nu_eld*pi*1.59155e-7; % rad/s*conversion_mhz nu1 - HTA pulse power - in MHz (DOI 10.1007/s00723-017-0927-4)
Exp.Tm = 0.13; % decay time of EDNMR nutations in us
Exp.tHTA = 5; % HTA pulse length in us - long low power pulse
Exp.Q = 1; % Q0 of the cavity (set 1 for no frequency dependence)
Exp.Temperature = 4.8; % Kelvin
Exp.CrystalOrientation = []; % powder
Exp.Harmonic = 0; 

% Options
Opt.nKnots = 31;
Opt.Threshold.Probe = 1e-4;
Opt.Threshold.Pump = 1e-4;
Opt.Output = 'summed';

%% Spin system (Ref 24 Olivia's MChem report) Final HFC values - isotropic
Sys1.S = 3; 
Sys1.g = [1.99]; 
Sys1.Nucs = '55Mn, 1H'; 
Sys1.D = -3e4*[0.281 0.0461]; %S3MEOH values obtained from reference 24 of report
Sys1.A = [-104, -104, -97.5; 1 1 1];
Sys1.lwEndor = 5; 

Sys2.S = 3; 
Sys2.g = [1.99]; 
Sys2.Nucs = '55Mn, 1H'; 
Sys2.D = -3e4*[0.281 0.0461]; %S3MEOH values obtained from reference 24 of report
Sys2.A = [-96.5, -96.5, -91.3; 1 1 1];
Sys2.lwEndor = 5; 

Sys3.S = 3; 
Sys3.g = [1.99]; 
Sys3.Nucs = '55Mn, 1H'; 
Sys3.D = -3e4*[0.281 0.0461]; %S3MEOH values obtained from reference 24 of report
Sys3.A = [-2.4; 1];
Sys3.lwEndor = 5; 

Sys4.S = 3; 
Sys4.g = [1.99]; 
Sys4.Nucs = '55Mn, 1H'; 
Sys4.D = -3e4*[0.281 0.0461]; %S3MEOH values obtained from reference 24 of report
Sys4.A = [-7.2; 1];
Sys4.lwEndor = 5; 

%% Simulations of mononuclear systems
% Calculate EDNMR spectrum Mn 1
[x_a1,MnSim_a1] = horseradish_github(Sys1,Exp,Opt); 
MnSim_A1=MnSim_a1/max(MnSim_a1);

% Calculate EDNMR spectrum Mn 2
[x_a2,MnSim_a2] = horseradish_github(Sys2,Exp,Opt); 
MnSim_A2=MnSim_a2/max(MnSim_a2);

% Calculate EDNMR spectrum Mn 3
[x_a3,MnSim_a3] = horseradish_github(Sys3,Exp,Opt); 
MnSim_A3=MnSim_a3/max(MnSim_a3);

% Calculate EDNMR spectrum Mn 4
[x_a4,MnSim_a4] = horseradish_github(Sys4,Exp,Opt); 
MnSim_A4=MnSim_a4/max(MnSim_a4);

save('EDNMR_Wband_1p98T_OMalley_simulation_V9.mat')
