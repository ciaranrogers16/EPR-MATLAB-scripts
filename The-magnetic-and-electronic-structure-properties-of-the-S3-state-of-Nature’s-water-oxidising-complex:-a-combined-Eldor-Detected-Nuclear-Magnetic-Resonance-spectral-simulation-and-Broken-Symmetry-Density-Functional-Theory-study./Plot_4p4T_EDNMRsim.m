%% Plotting of 1D EDNMR Simulation and correcting for bandwidth
%  4.4T
%
% ciaran.roger@manchester.ac.uk
% =========================================================================
clear all
close all

h=figure(1);
clf
hold on

%% Simulation plot - O'Malley
load('EDNMR_Wband_4p4T_OMalley_simulation_V9')
Sum_abcd_gauss_Omalley = MnSim_a1 + MnSim_a2 + MnSim_a3 + MnSim_a4; 
Sum_abcd_gauss_Omalley=Sum_abcd_gauss_Omalley./max(Sum_abcd_gauss_Omalley);
% plot(x_a,Sum_abcd_gauss_Omalley);

% Apply a gaussian
LOR=lorentzian(x_a1,0,400);
LOR=LOR/max(LOR);
Sum_abcd_gauss_Omalley_gauss=1.75*LOR.*Sum_abcd_gauss_Omalley;
% plot(x_a1,LOR+2.25,'M--');
plot(x_a1,Sum_abcd_gauss_Omalley_gauss+2,'r');

%% Simulation plot - Cox
load('EDNMR_Wband_4p4T_Chyrsina_simulation_norm')
Sum_abcd_gauss_chyrsina = MnSim_a1 + MnSim_a2 + MnSim_a3 + MnSim_a4; 
Sum_abcd_gauss_chyrsina=Sum_abcd_gauss_chyrsina./max(Sum_abcd_gauss_chyrsina);
% plot(x_a,Sum_abcd_gauss_1);

% Apply a gaussian
LOR=lorentzian(x_a1,0,400);
LOR=LOR/max(LOR);
Sum_abcd_gauss_chyrsina_gauss=1.75*LOR.*Sum_abcd_gauss_chyrsina;
plot(x_a1,Sum_abcd_gauss_chyrsina_gauss+1,'b');

%% Experimental plot
% Load W Band grabit plot
load('EDNMR_S3MeOH_4_4T_experimental.mat')
plot(EDNMR_Grabit_experimental(:,1),EDNMR_Grabit_experimental(:,2)./max(EDNMR_Grabit_experimental(:,2)),'k');

%% Adjusting figure
ylabel('Hole Intenisty / a.u.')
xlabel('( \nu_{HTA} - \nu_{Obs} ) / MHz')
axis([0 300 -0.1 4.25]);
set(gca,'YTick',[])
title('4.4 T - S_{3}^{MeOH}')
%legend('This work','Chrysina et al. simulation', 'Experimental')
h=make_plot_nice_ENDOR(h);