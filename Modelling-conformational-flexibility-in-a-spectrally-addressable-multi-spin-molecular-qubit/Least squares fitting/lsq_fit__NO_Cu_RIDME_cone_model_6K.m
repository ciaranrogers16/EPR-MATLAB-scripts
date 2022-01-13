%% Least squares fitting of simulated data to experimental data
% NO - Cr7Ni ring RIDME simulations - DA340
%
% Please ensure the subfolders ./Simulations and ./Simulations_for_fitting
% exist prior to running
%
% alice.bowen@manchester.ac.uk
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

%% Input experimental data and select field positions
% note this code assumes 5 traces to fit - if there are more or less
% traces to fit you will need to edit each section to add or remove
% variables.

% please ensure the subfolders ./Simulations and ./Simulations_for_fitting
% exist prior to running

%% Load experimental data
[x1,y1,p1] = textread('.\Processed_DA340_6K_34GHz_5P_RIDME_12104G_MPFUs_20000000d4_4n_total_fit.dat','%f %f %f');

% Can set up variable field as this is the variable we will loop over
% edit for number of datsets to be fitted simultaneously
% make sure order matches up to the x1, x2, x3 order of files loaded above.

Experiment.Field = [1209.5];

% define length of traces to use for fitting:

fitindex = [50; ...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    50 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    100 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    175 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ;...
    245 ];

%% set zero times and do background correction %%

%plot results?
doplot = 1; % logical operator 0 or 1
background_correct = 0; % logical operator 0 or 1 - 0 if data is already background corrcted
modulation_depth_removal = 1; % logical operator 0 or 1 - 1 removes modulation dpeth for fitting to account for experimental differences accross bandwidth etc.

if background_correct == 1
    % fit polynomial background - note the index of x1 and y1 define the
    % background region to fit - remove this part if the data is already
    % background corrected
    coef_b1 = polyfit(x1(60:end),y1(60:end),5);
    b1 = (coef_b1(1)*x1.^5)+(coef_b1(2)*x1.^4)+(coef_b1(3)*x1.^3)+(coef_b1(4)*x1.^2)+(coef_b1(5)*x1.^1)+coef_b1(6);
    
    % Plot results
    if doplot == 1
        figure(1)
        hold on
        plot(x1,y1)
        plot(x1,b1)
        
        figure(2)
        plot(x1,y1./b1)
    end
    
    % Normalise traces
    y1c = (y1./b1);
    y1c = y1c./max(y1c);
    
    if doplot == 1
        figure(10)
        plot(x1,y1c)
    end
    
else
    y1c = y1./max(y1);
    
    if doplot == 1
        figure(15)
        plot(x1,y1c)
    end
end % background correction


% Set zero time and truncate traces for fitting - this assumes all traces
% have the same time step and same zero time - if this is not the case
% please edit as needed.
time_step = x1(2)-x1(1);
[M,I] = max(y1);
zero_time = x1(I);

% Correct time axis
x1 = x1 - zero_time;

% Truncate the traces to remove time before zero time.
[M,I] = min(abs(x1));
x1_fit = x1(I:end);
y1_fit = y1c(I:end);

% Normalise traces
y1_fit = y1_fit./max(y1_fit);

%% Calculate modulation depths incase you need this for later work %%
% Calculate modulation depths of experimental traces by averaging last
% n points - you need to vary this () depending on how long the trace is and
% where the oscialations go to zero
num_points_MD = 25; % number of points to average modulation depth over.

MD1 = 1 - sum(abs(y1_fit(end-(num_points_MD-1):end)))./num_points_MD;

if modulation_depth_removal == 1
    
    % As there are often inconsistent ratios we remove the modulation depth for fitting
    y1_fit_zero = y1_fit - (1-MD1)*ones(length(y1_fit),1);
    y1_fit_zero = y1_fit_zero./max(y1_fit_zero);
    
    if doplot ==1
        figure(20)
        hold on
        plot(x1_fit, y1_fit_zero)
    end
    
else
    % Fit without modulation depth removal
    y1_fit_zero= y1_fit;
end

fprintf('Press any key to continue'), pause % pause at this point to review plots

%% Do fitting of time traces

% This code assumes that you have saved the files using two loops
% J is a strucutral parameter
% I is the experimetnal field parameter
fittingN = 90; % number of itterative fitting loops
maxJ = 5000; % largest J
maxI = length(Experiment.Field); % largest I
D1 = maxJ; % number of simulations to be fitted

% Load first file to check it loads ok - here the files are in a subfolder
% called Simulations
I = Experiment.Field(1); % Field
J = 1;
load(['.\Simulations_NO_Cu_RIDME_cone_model_6K\Sim_J_' num2str(J) '.mat']);

% Fitting weighting factor - if some data sets are of better or worse
% quality you may wish to weight their contibution more or less to the
% overall fit.
w1 = 1;

%set zero axes for collecting data - one per datset that needs to be fitted
DIFF = zeros(D1,length(Experiment.Field));

% save best time domain traces for plotting
BEST1 = zeros(length(x1_fit),1);

%save best frequency domain traces for plotting
frequency_length = length(f_sum_bin);% This code assumes you have used the time-domain method.
% You may need to change which parameter you use depending on the type of calcualtion method you have used.

BEST_Freq_1 = zeros(1,frequency_length);

%define experimental data for fitting
E1 = y1_fit_zero;

bestmin = zeros(fittingN,3);
tic
for N = 1:fittingN %edit this to edit the number of fitting loops
    display(N)
    DIFF = zeros(D1,1);
    if N == 1 % this first loop is used to remove the modulation depth if you are fitting without the modulation depth and save datasets in a smaller number of files
        ss = 1; % simulation variable
        for J = 1:maxJ
            % load simfiles and save y_deer for each field position
            
            load(['.\Simulations_NO_Cu_RIDME_cone_model_6K\Sim_J_' num2str(J) '.mat']);
            Y1 = y_FTdeer_bin;
            F1 = f_sum_bin;
            V1 = v_b_bin;
            
            if modulation_depth_removal == 1
                %calculate and correct for modulation depths
                MDsimY1 = 1 - sum(abs(Y1(end-(num_points_MD-1):end)))./num_points_MD;
                
                %Remove modulation depth for fitting
                Y1_zero = Y1 - (1-MDsimY1)*ones(1, length(Y1));
                Y1_zero  = Y1_zero ./Y1_zero(1); %normalize
                
            elseif modulation_depth_removal == 0
                Y1_zero = Y1;
            end
            
            % Redefine time axis for fitting - make it the correct length and step size
            t_fit=[t];
            trace_Y1_fit = [Y1_zero'];
            
            % Truncate data and time domains to the correct length
            T1 = t_fit(1:length(x1_fit));
            
            trace_Y1_fit = trace_Y1_fit(1:length(x1_fit));
            
            % add in previous best fits - in the first loop this is zeros so here it is
            % a sanity chack to make sure variables BEST* are of correct size
            Y1 = trace_Y1_fit+BEST1;
            
            %normalise
            Y1 = Y1./max(Y1);
            
            %calculate the total difference between experiment and data
            %for fitting range
            DIFF(ss) = (sum(w1*sqrt((real(E1(1:fitindex(N,1))-Y1(1:fitindex(N,1)))).^2)));
            
            save(['./Simulations_NO_Cu_RIDME_cone_model_6K_for_fitting/Sim_J_' num2str(J) '_I_all.mat'],...
                'Y1','T1','F1','V1');
            Index_all(ss,:) = [J];
            display(J)
            ss = ss+1; % Increment variable
        end % J loop
    else
        ss = 1; % simulation variable
        for J = 1:maxJ
            display(J)
            %load simfile
            load(['./Simulations_NO_Cu_RIDME_cone_model_6K_for_fitting/Sim_J_' num2str(J) '_I_all.mat']);
            Y1_fit = Y1;
            
            Y1 = Y1_fit +BEST1;
            
            %normalise
            Y1 = Y1./max(Y1);
            
            %DIFF(ss) = [(sum(sqrt((E1-Y1).^2))];
            DIFF(ss) = (sum(w1*sqrt((real(E1(1:fitindex(N,1))-Y1(1:fitindex(N,1)))).^2)));
            
            ss = ss+1; %Increment variable
        end
    end
    
    % Find minimum in the difference vector
    [M(N),In(N)] = min(real(DIFF));
    M_norm(N) = M(N)./sum(fitindex(N,:));
    
    % load simfile corresponding to minimum in difference vector
    % K = ceil(In(N)/maxK);
    % J = ((In(N)/maxK)-(K-1))*maxJ;
    J_store(N) = Index_all(In(N),1);
    load(['./Simulations_NO_Cu_RIDME_cone_model_6K_for_fitting/Sim_J_' num2str(J_store(N)) '_I_all.mat']);
    
    % add new best fit into existing best fit data
    BEST1 = Y1+BEST1;
    
    %save frequency part too
    BEST_Freq_1 = BEST_Freq_1+F1;
    
    %calculate the total difference btwween data and best fit
    BEST1_norm = BEST1./max(BEST1);
    
    Total_DIFF(N)=(sum(w1*sqrt((real(E1)-BEST1_norm(1:length(E1))).^2)));
end
toc

% Normalise best fit traces
BEST1 = BEST1./max(BEST1);

% Normalise best fit dipolar spectra
BEST_Freq_1 = BEST_Freq_1./max(BEST_Freq_1);

%% Save workspace variables
save('Data_Simulations_NO_Cu_RIDME_cone_model_6K_full_index.mat');

%% Plot all fitting data
f=figure(106);
hold on
plot([1:fittingN],M,'o-')
xlabel('Fitting iteration')
ylabel('Total difference between Exp and Sim')
f=make_plot_nice_FormFactor(f);
f.Renderer='opengl';

f=figure(107);
hold on
plot([1:fittingN],M_norm,'o-')
xlabel('Fitting iteration')
ylabel('Total difference between Exp and Sim')
f=make_plot_nice_FormFactor(f);
f.Renderer='opengl';

f=figure(108);
hold on
plot([1:fittingN], Total_DIFF,'o-')
xlabel('Fitting iteration')
ylabel('Total difference between Exp and Sim')
f=make_plot_nice_FormFactor(f);
f.Renderer='opengl';

f=figure(109);
hold on
plot(T1,E1,'k')
plot(T1,BEST1)

xlabel('Time / µs')
ylabel('Intensity')
legend('Experimental Cu-NO RIDME, 6 K','Best fit')
f=make_plot_nice_FormFactor(f);
f.Renderer='opengl';