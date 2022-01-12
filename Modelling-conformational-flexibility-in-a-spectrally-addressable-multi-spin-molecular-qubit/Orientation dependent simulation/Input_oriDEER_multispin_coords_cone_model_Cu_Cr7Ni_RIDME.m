%% Orientation Dependent Dipolar Simulation Input loop
%  oriDEER/RIDME simulation input
%
%%  DA340 Cu(II)-Cr7Ni
%   Copper detection params. from DFT EPR-II calc.
%   Cr7Ni params. defined geometrically from single crystal XRD structure
%
%   ciaran.rogers@manchester.ac.uk
%==========================================================================
clear all
close all

%% Import files
NAtoms=textread('CuTPPI_input_file_oriDEER.xyz','%d',1);
Atoms=cell(NAtoms,1);
XYZ=zeros(NAtoms,3); % Angstrom, cartesian coordinate system.

[Atoms,XYZ(:,1),XYZ(:,2),XYZ(:,3)]=...
    textread('CuTPPI_input_file_oriDEER.xyz','%s %f %f %f','headerlines',2);

XYZ=XYZ/10; % Angstrom -> nm

%% Copper(II) - coordinates of Cu in g frame of Cu - scaled, moved to origin
cu_det = (([XYZ(5,3)-XYZ(5,3),XYZ(5,2)-XYZ(5,2),XYZ(5,1)-XYZ(5,1);
    XYZ(2,3)-XYZ(5,3),XYZ(2,2)-XYZ(5,2),XYZ(2,1)-XYZ(5,1);
    XYZ(3,3)-XYZ(5,3),XYZ(3,2)-XYZ(5,2),XYZ(3,1)-XYZ(5,1);
    XYZ(4,3)-XYZ(5,3),XYZ(4,2)-XYZ(5,2),XYZ(4,1)-XYZ(5,1);
    XYZ(1,3)-XYZ(5,3),XYZ(1,2)-XYZ(5,2),XYZ(1,1)-XYZ(5,1)]./0.144).*0.1996);

cu_centre = cu_det(1,:)';

%% Cr7Ni ring - Average NH2+...Metal distance 0.45 nm
Cr7Ni_centre = ([cu_centre(1,:)] + [2.34 0 0])';

% Rotate +z +45 degrees - Active rotations (counter-clockwise)
rot_spin_centre_1z = cosd(45).*(Cr7Ni_centre(3,:)+0.45) - sind(45).*0;
rot_spin_centre_1y = sind(45).*(Cr7Ni_centre(3,:)+0.45) + cosd(45).*0;

% Rotate +z +135 degrees
rot_spin_centre_2z = cosd(135).*(Cr7Ni_centre(3,:)+0.45) - sind(135).*0;
rot_spin_centre_2y = sind(135).*(Cr7Ni_centre(3,:)+0.45) + cosd(135).*0;

% Rotate -z +45 degrees
rot_spin_centre_3z = cosd(45).*(Cr7Ni_centre(3,:)-0.45) - sind(45).*0;
rot_spin_centre_3y = sind(45).*(Cr7Ni_centre(3,:)-0.45) + cosd(45).*0;

% Rotate -z +135 degrees
rot_spin_centre_4z = cosd(135).*(Cr7Ni_centre(3,:)-0.45) - sind(135).*0;
rot_spin_centre_4y = sind(135).*(Cr7Ni_centre(3,:)-0.45) + cosd(135).*0;

Cr7Ni_pump = [Cr7Ni_centre(3,:)+0.45,Cr7Ni_centre(2,:),Cr7Ni_centre(1,:);
    Cr7Ni_centre(3,:)-0.45,Cr7Ni_centre(2,:),Cr7Ni_centre(1,:);
    Cr7Ni_centre(3,:),Cr7Ni_centre(2,:)+0.45,Cr7Ni_centre(1,:);
    Cr7Ni_centre(3,:),Cr7Ni_centre(2,:)-0.45,Cr7Ni_centre(1,:);
    rot_spin_centre_1z,rot_spin_centre_1y,Cr7Ni_centre(1,:)';
    rot_spin_centre_2z,rot_spin_centre_2y,Cr7Ni_centre(1,:)';
    rot_spin_centre_3z,rot_spin_centre_3y,Cr7Ni_centre(1,:)';
    rot_spin_centre_4z,rot_spin_centre_4y,Cr7Ni_centre(1,:)']; % coords

%% Cr7Ni overlay at origin
% Rotate +x +45 degrees - Active rotations (counter-clockwise)
rot_cu_centre_1x = cosd(45).*(cu_centre(1,:)+0.45) - sind(45).*0;
rot_cu_centre_1y = sind(45).*(cu_centre(1,:)+0.45) + cosd(45).*0;

% Rotate +x +135 degrees
rot_cu_centre_2x = cosd(135).*(cu_centre(1,:)+0.45) - sind(135).*0;
rot_cu_centre_2y = sind(135).*(cu_centre(1,:)+0.45) + cosd(135).*0;

% Rotate -x +45 degrees
rot_cu_centre_3x = cosd(45).*(cu_centre(1,:)-0.45) - sind(45).*0;
rot_cu_centre_3y = sind(45).*(cu_centre(1,:)-0.45) + cosd(45).*0;

% Rotate -x +135 degrees
rot_cu_centre_4x = cosd(135).*(cu_centre(1,:)-0.45) - sind(135).*0;
rot_cu_centre_4y = sind(135).*(cu_centre(1,:)-0.45) + cosd(135).*0;

Cr7Ni_pump_overlay = [cu_centre(1,:)+0.45,cu_centre(2,:),cu_centre(3,:);
    cu_centre(1,:)-0.45,cu_centre(2,:),cu_centre(3,:);
    cu_centre(1,:),cu_centre(2,:)+0.45,cu_centre(3,:);
    cu_centre(1,:),cu_centre(2,:)-0.45,cu_centre(3,:);
    rot_cu_centre_1x,rot_cu_centre_1y,cu_centre(3,:)';
    rot_cu_centre_2x,rot_cu_centre_2y,cu_centre(3,:)';
    rot_cu_centre_3x,rot_cu_centre_3y,cu_centre(3,:)';
    rot_cu_centre_4x,rot_cu_centre_4y,cu_centre(3,:)'];  % coordinates of Cr7Ni overlayed on centre of Copper spin

%% Input
%------------P vectors relating the location of the centers --------------%
% Each pump center requires an individual p.Nx vector
load('.\p_vectors_only.mat')

% How the p vector is defined:
% p(1) - phi x-y plane (rads)
% p(2) - theta from z-axis (rads)
% p(3) - r distance displace in nm
% p(4) - Euler-angle alpha (rads) for entrie pump system - coords and Ham
% p(5) - Euler-angle beta
% p(6) - Euler-angle gamma

%--------------------------INPUT DATA ARRAY OPT---------------------------%

opt.check_coordinates = 1;  % calls plotting program to plot the spin centers and tensors
opt.calcualte_DEER = 1;     % calls the main OriDEERsim code to calcualte DEER trace

% Experimental PARAMETERS
opt.Det_mw=34;           % Detection pulse frequency in GHz
opt.Pump_mw=29.25;       % Pumping pulse frequency in GHz (in RIDME - doesn't matter as BW of pump pulse infinite)
opt.B_obs=1187.1;        % field of the experiment in mT + b_corr

% Detection pulse PARAMETERS
opt.Det_t180=24;          % time for pi detection pulses (ns)- edit to match experiment
opt.Det_cutoff=0.001;     % cut off value for detection pulse resonances (vs angles on the grid)  - don't edit
opt.Det_calc_ori=1;       % calculate detection/Pump orientations if ==1 - don't edit
opt.Det_pulse_profile='sinc5';  % 'sinc2', 'sinc5', 'gaussian' or 'user_defined' - selects the profile of the Detection pulse - don't edit unless you are using a shaped pulse
%opt.Det_Pulse_Profile = ;  % needs to be defined only if opt.Det_pulse_profile = 'user_defined' allows the inclusion of shaped pulses. Needs to be the same lenght as opt.B0
opt.Det_resfields_Perturb=0;    % pertubation only for Easyspin option - don't edit
opt.Det_bandwidth=1e9;          % bandwidth around detection pulses (MHz)- mimics resonator bandwidth - don't edit

% Detection system PARAMETERS - EDIT for your system
opt.Det_coordinates = cu_det'; % coords of the detection center
Cu_A = [-60.4642 -423.9036]; % MHz
opt.Det_Ham=struct('S',1/2,'g',[2.0425 2.1604],'gFrame',[0 0 0],...
    'HStrain',[55.7976 127.6371],'Nucs','65Cu','A',[Cu_A],...
    'lwpp',[1.0594 1.6886]);
opt.Det_spindensity=[0.1 0.1 0.1 0.1 0.6]; % spin density for each point/nuclei
opt.Det_magnetophotoselection = 0;

% Pump pulse PARAMETERS (for all pump centers)
opt.Pump_t180=0.00001;         % time for pi Pump pulse (ns) - edit to match you experiment - in RIDME 0.00001 ns is effectively infinite BW
opt.Pump_excitation_fraction_cutoff=0.01;  % cut off value for (Pump pulse*site 2 resonances) overlap integral - don't edit
opt.Pump_pulse_profile='sinc2';  % 'sinc2', 'gaussian' or 'user_defined' - selects the profile of the Pump Pulse
%opt.Pump_Pulse_Profile = ;  % needs to be defined only if opt.Det_pulse_profile = 'user_defined' allows the inclusion of shaped pulses. Needs to be the same lenght as opt.B0
opt.Pump_resfields_Perturb=0;   % pertubationonly for Easyspin option
opt.Pump_bandwidth=1e9;         % bandwidth around pump pulses (MHz) - mimics resonator bandwidth - don't edit

% Pump system PARAMETERS - Edit for your system
opt.Pump_number = 1; % define the number of pump centers, 1 center gives pairwise interaction only etc

% Pump system PARAMETERS (specific to center 1) - Edit for your system
opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
opt.Pump_spindensity.N1=[0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125]; % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coords for the center of the moeity in the starting frame used as center of rotation and translation
opt.Pump_coordinates.N1=Cr7Ni_pump_overlay'; % coords of the first pump center
opt.Pump_magnetophotoselection = 0; 

% Pump system PARAMETERS (specific to center 2) - Edit for your system
% opt.Pump_Ham.N2=struct('S',1/2,'g',[2.0101 2.00708 2.00344],'gFrame',[0 0 0],'HStrain',[15.795 13.886 19.954],'Nucs','14N','A',[12.815 19.6701 102.061],'AFrame',[0 0 0]);
% opt.Pump_spindensity.N2=[0.5 0.5]; % spin density for each point/nuclei
% opt.Pump_coordinate_centre.N2=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
% opt.Pump_coordinates.N2=[0 0 0; 0.5 0 0]';% coordinates of the second pump center - best to include these in the most anisotripc frame of the pump hamiltonian
%
% % Pump system PARAMETERS (specific to center 3) - Edit for your system
% opt.Pump_Ham.N3=struct('S',1/2,'g',[2.0101 2.00708 2.00344],'gFrame',[0 0 0],'HStrain',[15.795 13.886 19.954],'Nucs','14N','A',[12.815 19.6701 102.061],'AFrame',[0 0 0]);
% opt.Pump_spindensity.N3=[0.5 0.5]; % spin density for each point/nuclei
% opt.Pump_coordinate_centre.N3=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
% opt.Pump_coordinates.N3=[0 0 0; 0.5 0 0]';% coordinates of the third pump center - best to include these in the most anisotripc frame of the pump hamiltonian
%
% % Pump system PARAMETERS (specific to center 4) - Edit for your system
% opt.Pump_Ham.N4=struct('S',1/2,'g',[2.0101 2.00708 2.00344],'gFrame',[0 0 0],'HStrain',[15.795 13.886 19.954],'Nucs','14N','A',[12.815 19.6701 102.061],'AFrame',[0 0 0]);
% opt.Pump_spindensity.N4=[0.5 0.5]; % spin density for each point/nuclei
% opt.Pump_coordinate_centre.N4=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
% opt.Pump_coordinates.N4=[0 0 0; 0.5 0 0]';% coordinates of the fourth pump center - best to include these in the most anisotripc frame of the pump hamiltonian
%
% % Pump system PARAMETERS (specific to center 5) - Edit for your system
% opt.Pump_Ham.N5=struct('S',1/2,'g',[2.0101 2.00708 2.00344],'gFrame',[0 0 0],'HStrain',[15.795 13.886 19.954],'Nucs','14N','A',[12.815 19.6701 102.061],'AFrame',[0 0 0]);
% opt.Pump_spindensity.N5=[0.5 0.5]; % spin density for each point/nuclei
% opt.Pump_coordinate_centre.N5=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
% opt.Pump_coordinates.N5=[0 0 0; 0.5 0 0]';% coordinates of the fifth pump center - best to include these in the most anisotripc frame of the pump hamiltonian

% General settings for both spin systems and experiment
opt.lineshape='gaussian';       % line shape form making the spectrum can be 'gaussian' or 'lorentian' - don't edit
opt.nKnots=31;                  % DEER orientation selection defines the number of knots - make this ca. 5 if you are graphing the excitation fraction otherwise 90 degrees makes one point every 1 degree is usually a good starting point. If the time trace comes back into phase then make this larger - making this smaller makes the calculation quicker
opt.Sym='Ci';                   % symmetry for points on grid used for orientation calculation - don't edit
t=(0:0.008:2)';                 % DEER time trace in microseconds (us)it is good to match the time step to your experimental data if you want to fit one to the other
opt.b0=(1050:120/2047:1450)';   % field vector must cover CW EPR field range of both spin systems - make sure your spin system si resonant in this window.

% logical operators for calculation parts,
opt.calc_CWEPR=1;               % calculate EPR sepctrum on each cycle? should only be set to zero if the program is running one detection pulse type/postion within a folder - if it is set to zero the caluation for the CW spectrum is not performed and loaded from a file. If you are unsure leave this at 1.
opt.Det_calc_ori=1;             % calculate detection/Pump orientations? should only be set to zero if the program is running one detection pulse type/postion within a folder - if it is set to zero the caluation for the detection resonances is not performed and loaded from a file.  If you are unsure leave this at 1.
opt.Two_spin_frequency_bin = 1; % calculate the DEER trace using frequency binning. Works only if there is only one pump center - can also be used for checking
opt.Multi_spin_time_calc = 0;   % calculate the DEER trace using time dependent cos functions - works for all numbers of pump centers and calculated the multispin interactions

% options for processing of time domain calculated (multi spin) data to get frequency domain.
opt.window = 1;                 % apply a window function before FFT
opt.gauss_wid = 0;              % use a gaussian window
opt.hann_wid =1;                % use a hanning window
opt.hamm_wid =0;                % use a hamming window
opt.exp_wid =0;                 % use an exponential decay window
opt.zero_fill =  1;             % use zero-filling for FFT - you can perform hamming without zerofil but not zero fill without hamming
opt.zero_fill_number = 8192;    % number of points to zero fill to.

% logical operators for plotting and printing
opt.graph_CWEPR=0;                    % for graphing continuuos wave EPR spectrum - graphs the CW spectrum of detection center and pump center(s) with the profiles of the detection and pump pulses - turn off if you are looping over multiple inputs as it requires manual input to continu calculation
opt.Det_graph_ori=0;                  % plot orientation selection of detection spin (plots coloured spheres)- turn off if you are looping over input values.
opt.Det_graph_excitation_fraction=0;  % for graphing the excitation fraction of detection spins - only use for checking and with a small number of knots
opt.Pump_graph_excitation_fraction=0; % for graphing the excitation fraction of pump spins - only use for checking and with a small number of knots
opt.graph_DEERtrace_FFT=0;            % graph final DEER trace / dipolar spectrum - can be left turned off if you want to plot traces in your own style or after loop of input values.
opt.Print_coord=0;                    % print coordinates of system to the screen - these are not currently saved by the program as the depend on the grometry of the system.
opt.quiet = 0;                        % turn on or off all printing to screen - should be set to 1 to remove printing to screen

% resonant field calculation methods - Do not edit
opt.Det_resfields_method='easyspin';   % currently only 'easyspin' available
opt.Pump_resfields_method='easyspin';  % currently only 'easyspin' available
opt.f_dd_method='secular';             % currently only 'secular' available


%% ---------------------------RUN THE CALCULATION---------------------------

    for J = 1:length(p_pumpCr7Ni_detCu)
        p.N1 = p_pumpCr7Ni_detCu(:,J)';
        if opt.calcualte_DEER == 1
            
            [y_deer v_b f_sum ft t_wid y_deer_wid t_d_bin y_FTdeer_bin v_b_bin f_sum_bin f_mod_bin]= OriPDSSim_multispin_magnetophotoselection(p,t,opt);
            
            % OUTPUT
            %  y_deer   - time-domain trace
            %  v_b      - frequency axis
            %  f_sum    - frequency domain (fft of y_deer)
            %  ft       - individual time component for pairwise interactions
            %  t_wid    - time domain after window function is applied
            %  y_deer_wid - time-domain trace after window function is applied
            %  t_d_bin  - time domain for binning calculation
            %  y_FTdeer_bin - time trace for binning calculation
            %  v_b_bin  - frequency domain for binning calculation
            %  f_sum_bin - frequencies calcualted in binning calculation includigng
            %  zero frequencies
            %  f_mod_bin - frequencies calcualted in binning calculation excluding zero
            %  frequencies
            
            % Plotting Dipolar Simulations
            %figure, hold on, plot(t,y_deer./max(y_deer),'b'), plot(t,sum(ft.N1)./max(y_deer),'r'),
            %figure, hold on, plot(v_b,abs(f_sum)./max(abs(f_sum)),'b'),plot(v_b_bin,f_mod_bin./max(f_mod_bin),'g')
            %figure, hold on, plot(t,y_deer./max(y_deer),'b'), plot(t_d_bin,y_FTdeer_bin,'g')
            %figure, hold on, plot(t_wid,y_deer_wid./max(y_deer_wid),'b')
              save(['./Simulations_Cu_Cr7Ni_RIDME_Cone_model_FINAL_UPDATED/Sim_J_',num2str(J),'.mat'],...
                    'y_FTdeer_bin','f_sum_bin','v_b_bin','t','p_pumpCr7Ni_detCu');
        end
    end