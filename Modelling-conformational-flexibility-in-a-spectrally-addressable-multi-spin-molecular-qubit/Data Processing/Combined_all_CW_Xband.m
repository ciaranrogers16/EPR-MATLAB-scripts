%% Plotting experimental and simulated CW X-band spectra
%
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

B_corr=24; % correction in Gauss via Strong Pitch calibration
%% 50 K
%
%
%% Load Raw Data
[B_50,spc_50] = eprload('.\DA340_50K_CW_9p3973GHz_5G_20dBAtt_6n');
B_50=B_50+B_corr; % field axis correction
B_mT_50=B_50/10; % Gauss2mT

%% Baseline Correction and Data Smoothing
spc_base_50=basecorr(spc_50,1,2); % (data to correct, dim, order of polynomial)
spc_base_50=real(spc_base_50);

%% Simulation

% Spin System
Cu_50.S = 1/2;
Cu_50.Nucs = 'Cu,14N,14N,14N';
Cu_50.g = [2.0628 2.2192];
Cu_A_50 =[64.3692 564.2601]; % MHz
N_A_50 = [40.0654 43.6985]; % MHz
Cu_50.A = [Cu_A_50; N_A_50; N_A_50; N_A_50];
Cu_50.HStrain=[27.2887 55.7070];
Cu_50.lwpp = [0.1491 0.0629]; % mT
Cu_50.weight=0.8594;

% Spin System - Nitroxide 
Nitroxide_50.g = [2.0057 2.0055 2.0011];
Nitroxide_50.Nucs = '14N';
Nitroxide_50.A = [21.4379 23.9962 101.3483]; % MHz
Nitroxide_50.lw = 1.2034; % mT
Nitroxide_50.weight = 1;


% Experimental Parameters
Exp_50.mwFreq = 9.3988; % in GHz
Exp_50.Range = [min(B_mT_50) max(B_mT_50)]; % in mT
Exp_50.nPoints = numel(B_mT_50); 
Exp_50.Temperature = 50; % Kelvin


% Optional Parameters
Opt_50.nKnots = 31;
Opt_50.Method = 'hybrid';

[B_sim_50,spc_sim_50]=pepper({Cu_50,Nitroxide_50},Exp_50,Opt_50);

%% 15 K 
%
%
%% Load Raw Data
[B_15,spc_15] = eprload('.\DA340_15K_CW_9p3989GHz_5G_20dBAtt_9n_tuned_correctly');
B_15=B_15+B_corr; % field axis correction
B_mT_15=B_15/10; % Gauss2mT

%% Baseline Correction and Data Smoothing
spc_base_15=basecorr(spc_15,1,2); % (data to correct, dim, order of polynomial)
spc_base_15=real(spc_base_15);

%% Simulation
clear Sys_50;
% Spin System
Cu_15.S = 1/2;
Cu_15.Nucs = 'Cu,14N,14N,14N';
Cu_15.g = [2.0628 2.2192];
Cu_A_15 =[64.3692 564.2601]; % MHz
N_A_15 = [40.0654 43.6985]; % MHz
Cu_15.A = [Cu_A_15; N_A_15; N_A_15; N_A_15];
Cu_15.HStrain=[27.2887 55.7070];
Cu_15.lwpp = [0.1491 0.0629]; % mT
Cu_15.weight=0.75;

% Spin System - Nitroxide 
Nitroxide_15.g = [2.0057 2.0055 2.0011];
Nitroxide_15.Nucs = '14N';
Nitroxide_15.A = [21.4379 23.9962 101.3483]; % MHz
Nitroxide_15.lw = 1.2034; % mT
Nitroxide_15.weight = 1;

 % Spin System - Cr7Ni Ring
Cr7Ni_15.S = 1/2;
Cr7Ni_15.g = [1.7826 1.8020 1.7463];
Cr7Ni_15.gStrain = [0.0050 0.0032 0.0054];
Cr7Ni_15.lwpp = [10]; % mT
Cr7Ni_15.weight = [0.25];

clear Exp_50;
% Experimental Parameters
Exp_15.mwFreq = 9.3988; % in GHz
Exp_15.Range = [min(B_mT_15) max(B_mT_15)]; % in mT
Exp_15.nPoints = numel(B_mT_15); 
Exp_15.Temperature = 15; % Kelvin

clear Opt_50;
% Optional Parameters
Opt_15.nKnots = 31;
Opt_15.Method = 'hybrid';

[B_sim_15,spc_sim_15]=pepper({Cu_15,Cr7Ni_15,Nitroxide_15},Exp_15,Opt_15);

%% 5.5 K
%
%
%% Load Raw Data
[B_5,spc_5] = eprload('.\DA340_5p5K_CW_9p3988GHz_5G_20dBAtt_10n');
B_5=B_5+B_corr; % field axis correction
B_mT_5=B_5/10; % Gauss2mT

%% Baseline Correction and Data Smoothing
spc_base_5 = basecorr(spc_5,1,2); % (data to correct, dim, order of polynomial)
spc_base_5 = real(spc_base_5);

clear Sys_15;
% Spin System
Cu_5.S = 1/2;
Cu_5.Nucs = 'Cu,14N,14N,14N';
Cu_5.g = [2.0628 2.2192];
Cu_A_5 =[64.3692 564.2601]; % MHz
N_A_5 = [40.0654 43.6985]; % MHz
Cu_5.A = [Cu_A_5; N_A_5; N_A_5; N_A_5];
Cu_5.HStrain=[27.2887 55.7070];
Cu_5.lwpp = [0.1491 0.0629]; % mT
Cu_5.weight=0.4;

% Spin System - Nitroxide 
Nitroxide_5.g = [2.0057 2.0055 2.0011];
Nitroxide_5.Nucs = '14N';
Nitroxide_5.A = [21.4379 23.9962 101.3483]; % MHz
Nitroxide_5.lw =  1.3028; % mT
Nitroxide_5.weight = 1;

 % Spin System - Cr7Ni Ring
Cr7Ni_5.S = 1/2;
Cr7Ni_5.g = [1.7875 1.8020 1.7482];
Cr7Ni_5.gStrain = [0.0050 0.0032 0.0054];
Cr7Ni_5.lwpp = [10]; % mT
Cr7Ni_5.weight = [2.5];

clear Exp_15;
% Experimental Parameters
Exp_5.mwFreq = 9.3988; % in GHz
Exp_5.Range = [min(B_mT_5) max(B_mT_5)]; % in mT
Exp_5.nPoints = numel(B_mT_5); 
Exp_5.Temperature = 5.5; % Kelvin

clear Opt_15;
% Optional Parameters
Opt_5.nKnots = 31;
Opt_5.Method = 'hybrid';

[B_sim_5,spc_sim_5]=pepper({Cu_5,Cr7Ni_5,Nitroxide_5},Exp_5,Opt_5);

%% Plotting and Graphical Rendering
h=figure(1); 
clf;
hold on; 

spc_sim_50=rescale(spc_sim_50,spc_base_50,'lsq1');
p1=plot(B_sim_50,spc_sim_50,'k.');
p2=plot(B_mT_50,spc_base_50-5,'color',[0.6157    0.0039    0.2588]);
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

spc_sim_15=rescale(spc_sim_15,spc_base_15,'lsq1');
p3=plot(B_sim_15,spc_sim_15-20,'k.');
p4=plot(B_mT_15,spc_base_15-25,'color',[0.9922    0.6863    0.3686]);
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

spc_sim_5=rescale(spc_sim_5,spc_base_5,'lsq1');
p5=plot(B_sim_5,spc_sim_5-40,'k.');
p6=plot(B_mT_5,spc_base_5-45,'color',[0.3686    0.3098    0.6353]);
set(get(get(p5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlabel('{\it B}_0 (mT)')
ylabel('Signal Intensity (a.u.)')
set(gca,'YTick',[])
set(gca,'TickDir','in')
axis([240 450 -300 200])
title ('X-band')
legend('50 K','15 K','5.5 K')
h=make_plot_nice_EDFS(h);
h.Renderer='opengl';

