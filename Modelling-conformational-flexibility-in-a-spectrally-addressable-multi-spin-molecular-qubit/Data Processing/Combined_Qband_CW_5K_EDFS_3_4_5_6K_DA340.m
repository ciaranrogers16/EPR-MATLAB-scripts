%% CW and EDFS Processing Template
% DA340; EDFS @ 34GHz ; 3K / 5K
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

%% Load Raw Data EDFS
% 3K
[B,spc,pars] = eprload('.\AB_CR024_DA340_EDFS_34GHz_3K_12132G_12_24_5n');
B_mT=B./10; % Gauss2mT
real_spc=real(spc);
real_spc_norm=real_spc./max(real_spc);

% 4 K
[B4,spc4,pars4] = eprload('.\DA340_EDFS_12_24_20pcAMP_HPA0dB_34GHz_4K_5n');
B_mT4=B4./10; % Gauss2mT
real_spc4=real(spc4);
real_spc_norm4=real_spc4./max(real_spc4);

% 5 K
[B5,spc5,pars5] = eprload('.\EDFS_DA340_3mm_5K_34GHz_10_20_20pcAWG_5n');
B_mT5=B5./10; % Gauss2mT
real_spc5=real(spc5);
real_spc_norm5=real_spc5./max(real_spc5);

% 6 K
[B6,spc6,pars6] = eprload('.\DA340_EDFS_6K_34GHz_12_24_75pcAmp_HPA_0dB');
B_mT6=B6./10; % Gauss2mT
real_spc6=real(spc6);
real_spc_norm6=real_spc6./max(real_spc6);

%% Load Raw Data CW
B_corr = 16;

[B_cw,spc]=eprload('.\AB_05_CRDA340_02_5K_13000Gcentre_5000G_100secs_20dB_2G_att18dB_nos5_bgrcorrect');
B_cw=B_cw+B_corr; % field axis correction
B_mT_cw=B_cw/10; % Gauss2mT

%% Baseline Correction
spc_base = basecorr(spc,1,3); % (data to correct, dim, order of polynomial)

%% Simulation
clear Sys
% Spin System
Cu.S = 1/2;
Cu.Nucs = 'Cu,14N';
Cu.g = [2.0496 2.2080];
Cu_A = [6.5109 -806.7639]; % MHz
N_A = [21.1392 77.4961]; % MHz
Cu.A = [Cu_A; N_A];
Cu.HStrain=[26.9524 366.6823];
Cu.lwpp = [5.8388 1.0586]; % mT
Cu.weight=2;

% Spin System - Cr7Ni Ring
Cr7Ni.S = 1/2;
Cr7Ni.g = [1.7826 1.8020 1.7463];
Cr7Ni.gStrain = [0.0050 0.0032 0.0054];
Cr7Ni.lwpp = [14.0759]; % mT
Cr7Ni.weight = [2.3992];

% Spin System - Nitroxide 
Nitroxide.g = [2.0059,2.0096,2.0019];
Nitroxide.Nucs = '14N';
Nitroxide.A = [21.5556,23.7915,100.0972]; % MHz
Nitroxide.lw = 1.6183; % mT
Nitroxide.weight = 1;

clear Exp;
% Experimental Parameters
Exp.mwFreq = 34; % in GHz
Exp.Range = [min(B_mT_cw) max(B_mT_cw)]; % in mT
Exp.nPoints = numel(B_mT_cw); 
Exp.Temperature = 5; % Kelvin

clear Opt;
% Optional Parameters
Opt.nKnots = 91;
Opt.Method = 'hybrid';

% Simulation
[B_sim,spc_sim]=pepper({Cu,Cr7Ni,Nitroxide},Exp,Opt);

%% Plotting and Graphical Rendering
h=figure(1); 
clf;
hold on; 

spc_base=spc_base./max(spc_base);
plot(B_mT_cw,spc_base,'k')
spc_sim=spc_sim./max(spc_sim);
spc_sim=rescale(spc_sim,spc_base,'lsq1');
plot(B_sim,spc_sim,'r')

plot(B_mT,real_spc_norm-2.3,'color',[0.3686    0.3098    0.6353])
plot(B_mT4,real_spc_norm4-2.1,'color',[0.6706    0.8706    0.6510])
plot(B_mT5,real_spc_norm5-1.9,'color',[0.9922    0.6863    0.3686])
plot(B_mT6,real_spc_norm6-1.7,'color',[0.6157    0.0039    0.2588])

xlabel('{\it B}_0 (mT)')
ylabel('Signal Intensity (a.u.)')
set(gca,'YTick',[])
set(gca,'TickDir','in')
axis([1050 1450 -10 10])
title ('Q-band')
legend('CW (5 K)','Simulation','EDFS (3 K)','EDFS (4 K)','EDFS (5 K)','EDFS (6 K)')

h=make_plot_nice_EDFS(h);
h.Renderer='opengl';