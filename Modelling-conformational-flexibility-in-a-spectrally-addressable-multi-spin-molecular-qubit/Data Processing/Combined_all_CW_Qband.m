%% Plotting experimental and simulated CW Q-band spectra
% 
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

B_corr=8; % correction in Gauss via Strong Pitch calibration

%% Load Raw Data
[B_20,spc_20] = eprload('.\AB_05_CRDA340_04_20K_13000Gcentre_5500G_110secs_20dB_2G_att18dB_nos5.dta');
B_20=B_20+B_corr; % field axis correction
B_mT_20=B_20/10; % Gauss2mT

[B_50,spc_50] = eprload('.\AB_05_CRDA340_05_50K_13000Gcentre_5500G_110secs_20dB_2G_att18dB_nos5.dta');
B_50=B_50+B_corr; % field axis correction
B_mT_50=B_50/10; % Gauss2mT

[B_100,spc_100] = eprload('.\AB_05_CRDA340_05_100K_13000Gcentre_5500G_110secs_20dB_2G_att18dB_nos5.dta');
B_100=B_100+B_corr; % field axis correction
B_mT_100=B_100/10; % Gauss2mT

%% Plotting and Graphical Rendering
h=figure(); 
clf;
hold on; 

spc_20=rescale(spc_20,spc_100,'maxabs');
plot(B_mT_20,spc_20,'r')
spc_50=rescale(spc_50,spc_100,'maxabs');
plot(B_mT_50,spc_50,'b')

plot(B_mT_100,spc_100,'m')

xlabel('{\it B}_0 (mT)')
ylabel('Signal Intensity (a.u.)')
set(gca,'YTick',[])
legend('5 K','20 K','50 K','100 K')

h=make_plot_nice(h);
h.Renderer='opengl';
