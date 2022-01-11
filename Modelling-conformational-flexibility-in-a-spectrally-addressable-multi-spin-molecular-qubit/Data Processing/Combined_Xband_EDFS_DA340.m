%% CW and EDFS Processing Template
% DA340; EDFS @ Xband ; 5 / 15 / 50 K
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

%% Load Raw Data EDFS
% 5.5 K
[B5,spc5,pars5] = eprload('.\AB_CR027_DA340_EDFS_5.5K_9_657GHz_20_40_5mm_tau_120');
B_mT5 = B5./10; % Gauss2mT
B_mT5=B_mT5-12.2;
real_spc5 = real(spc5);
real_spc_norm5 = real_spc5./max(real_spc5);

% 15 K
[B,spc,pars] = eprload('.\DA340_EDFS_15K_9p31GHz_20_40_HPA6dB');
B_mT=B./10; % Gauss2mT
real_spc=real(spc);
real_spc_norm=real_spc./max(real_spc);

% 50 K
[B50,spc50,pars50] = eprload('.\AB_031_DA340_EDFS_50K_9_637GHz_20_40_300tau_3n');
B_mT50=B50./10; % Gauss2mT
B_mT50=B_mT50-11.5;
real_spc50=real(spc50);
real_spc_norm50=real_spc50./max(real_spc50);

%% Plotting and Graphical Rendering
h=figure(1); 
clf;
hold on; 
plot(B_mT50,real_spc_norm50+1,'color',[0.6157    0.0039    0.2588])

plot(B_mT,real_spc_norm+0.5,'color',[0.9922    0.6863    0.3686])

plot(B_mT5,real_spc_norm5,'color',[0.3686    0.3098    0.6353])

xlabel('{\it B}_0 (mT)')
ylabel('Signal Intensity (a.u.)')
set(gca,'YTick',[])
set(gca,'TickDir','in')
axis([210 450 -10 10])
title ('X-band')
legend('50 K','15 K','5.5 K')

h=make_plot_nice_EDFS(h);
h.Renderer='opengl';