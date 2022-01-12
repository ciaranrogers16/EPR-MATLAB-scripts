function [y_deer v_b f_sum ft t_wid y_deer_wid t_d_bin y_FTdeer_bin v_b_bin f_sum_bin f_mod_bin]= OriPDSSim_multispin_magnetophotoselection(p,t,opt)
% --- OUTPUT ---
%  y_deer   - time-domain trace
%  v_b      - frequency axis
%  f_sum    - frequency domain (fft of y_deer)
%  ft       - indivudual time component for pairwise interactions
%  t_wid    - time domain after window function is applied
%  y_deer_wid - time-domain trace after window function is applied
%  t_d_bin  - time domain for binning calculation
%  y_FTdeer_bin - time trace for binning calculation
%  v_b_bin  - frequency domain for binning calculation
%  f_sum_bin - frequencies calcualted in binning calculation includigng
%  zero frequencies
%  f_mod_bin - frequencies calcualted in binning calculation excluding zero
%  frequencies
%
% --- INPUT ---
% p(1) - phi x-y plane   (radians)
% p(2) - theta  from z-axis (radians)
% p(3) - r distance displace in nm
% p(4) - Euler-angle alpha (radians) for entrie pump system - coordinates and Hamiltonian
% p(5) - Euler-angle beta
% p(6) - Euler-angle gamma
% t    - time axis in us
% -------------------------------------------------------------------------
if opt.quiet == 0
    fprintf('\n--- Simulation with OriLITTERSim ---')
end
% ------------------------------INITIAL CHECKS---------------------------
[r c]=size(opt.b0); if c>r, opt.b0=opt.b0'; end   % ensure the field axis is a column vector
% ------------------------------DEFINE CONSTANTS---------------------------
bohr = 9.27*10^-24;      % J/T
planck = 6.626*10^-34;   % J/s   [h/(2pi)=1.055*10^-34 J/s]
perm = 4*pi*10^-7;       %
k1=(bohr^2*perm)/(4*pi*planck); % collect constants together

% ---------------SET UP PUMP CENTERS RELATIVE TO DETECTION CENTER-----------------
% --- rotate Pump Spin Systems (opt.Pump_Ham.Nx.gFrame,AFrame,DFrame) and coordinates(opt.Pump_coordinates.Nx) into molecular frame ---
opt=Pump_translate_rotate(p,opt);
% --- translate distances (R_12), unit vectors (N_12) and spin densties(spindensity_12) for each pump center realtive to the detection center---
[R_12 N_12 spindensity_12 opt]=Pump_det_N_spindensities(p,opt);
% --- linewidth switch---
switch(opt.lineshape), case 'gaussian', opt.lineshape=1; case 'lorentzian', opt.lineshape=0; end

%-----------------DEFINE THE PROFILE OF THE PUMP PULSE-----------------
opt.dv=(opt.b0-opt.B_obs)*28.02494;   % dv axis in MHz (B_obs axis in mT to MHz)
dw=2*pi*opt.dv; %  to radians
switch(opt.Pump_pulse_profile)
    case 'sinc2' % density matrix result
        w1=2*pi/(2*opt.Pump_t180/1e3);      % Pump B1 pi pulse field strength (MHz to radians)
        w_eff=sqrt(dw.^2+w1.^2);
        P_Pump=(w1./w_eff).^2.*sin(w_eff*opt.Pump_t180/2/1e3).^2;
    case 'gaussian' % gaussian
        BW=(1/opt.Pump_t180/1e3)/28.02;     % gaussian Pump pulse excitation in mT (FWHH)
        P_Pump=DEER_Sim_specmake(opt.b0,opt.B_obs,1,BW,0);
        P_Pump=P_Pump/max(P_Pump);
    case 'sinc5' % density matrix result - usually used for detection pulse profile
        w1=2*pi/(2*opt.Pump_t180/1e3);       % Pump B1 pi pulse field strength (MHz to radians)
        w_eff=sqrt(dw.^2+w1.^2);
        P_Pump=(w1./w_eff).^5.*sin(w_eff*opt.Det_t180/2/1e3).^5;
    case 'user_defined' % User defined pulse profile to enable the use of shaped pulses
        P_Pump=opt.Pump_Pulse_Profile;
end
%-----------------DEFINE THE PROFILE OF THE DETECTION PULSE-----------------
switch(lower(opt.Det_pulse_profile))
    case 'sinc5' % density matrix result correct Detection pulse profile
        w1=2*pi/(2*opt.Det_t180/1e3);       % Pump B1 pi pulse field strength (MHz to radians)
        w_eff=sqrt(dw.^2+w1.^2);
        P_det=(w1./w_eff).^5.*sin(w_eff*opt.Det_t180/2/1e3).^5;
    case 'sinc2' % density matrix result
        w1=2*pi/(2*opt.Det_t180/1e3);       % Pump B1 pi pulse field strength (MHz to radians)
        w_eff=sqrt(dw.^2+w1.^2);
        P_det=(w1./w_eff).^5.*sin(w_eff*opt.Det_t180/2/1e3).^2;
    case 'user_defined' % User defined pulse profile to enable the use of shaped pulses
        P_det=opt.Det_Pulse_Profile;
end
%-------------------------Gaussian RESONATOR BANDWIDTH--------------------------
Det_bandwidth=DEER_Sim_specmake(opt.dv,0,1,opt.Det_bandwidth,1);   % axis MHz
Det_bandwidth=Det_bandwidth/max(Det_bandwidth);                    %normalisation
Pump_bandwidth=DEER_Sim_specmake(opt.dv,0,1,opt.Pump_bandwidth,1); % axis MHz
Pump_bandwidth=Pump_bandwidth/max(Pump_bandwidth);                 %normalisation
%--------------------------TOTAL PULSE RESPONSE------------------------
P_Pump=P_Pump.*Pump_bandwidth;
P_det=P_det.*Det_bandwidth;
% CALCULATION OF CW SPECTRA FOR orientation selection /diagnostics
if opt.calc_CWEPR==1, CWEPR_simulate(opt,P_Pump,P_det), end
if opt.graph_CWEPR==1, CWEPR_diagnostic(opt,P_Pump,P_det), end
%-------- DETECTIONS ORIENTATIONS -----------------
if opt.quiet == 0
    fprintf('\n1. DETECTION ORIENTATION CALCULATIONS\n')
end
if opt.Det_calc_ori==1
    ori1 = det_resfields_orientations(opt,P_det);
else
    load det_orientations  % variables: ori1 %opt1 Exp1 Sys1
end
if opt.quiet == 0
    fprintf('...Finished.')
end

%-------- PUMP ORIENTATIONS -----------------------
if opt.quiet == 0
    fprintf('\n2. PUMP ORIENTATION CALCULATIONS')
end
for pumpN = 1:opt.Pump_number   % loop over pump centers
    fn = sprintf('N%d', pumpN);
    switch(lower(opt.Pump_resfields_method)) % choose RESFIELDS METHOD
        case 'easyspin' % Easy spin general method
            ExpPump = struct('mwFreq',opt.Pump_mw,'Range',[min(opt.b0) max(opt.b0)],'CrystalOrientation',ori1.angles');%,'Temperature',opt.Pump_Triplet_pop);
            Pump_Ham_Easyspin = opt.Pump_Ham.(fn);   
            if opt.Pump_Ham.(fn).S>=1
            ExpPump.Temperature=opt.Pump_Triplet_pop.(fn);
            end
            if opt.Pump_magnetophotoselection ==1
            ExpPump.Ordering=opt.Pump_Ordering;
            Pump_Ham_Easyspin.HStrain = mt2mhz(max(Pump_Ham_Easyspin.lw)); % convert largest linewidth to Hstrain for resfields
            end
            if isfield(Pump_Ham_Easyspin,'HStrain')==0 && isfield(Pump_Ham_Easyspin,'gStrain')==0 && isfield(Pump_Ham_Easyspin,'AStrain')==0 
            Pump_Ham_Easyspin.HStrain = mt2mhz(max(Pump_Ham_Easyspin.lwpp));
            end
            if isfield(opt,'Pump_resfields_Perturb')==0, opt.Pump_resfields_Perturb=1; end
            [En.(fn), Amp.(fn), Wid.(fn), g_Pump.(fn)]=resfields_method_EasySpin(Pump_Ham_Easyspin,ExpPump,opt.Pump_resfields_Perturb); % En in mT
            xaxis=opt.b0;
    end
end
if opt.quiet == 0
    fprintf('...Finished.')
end
%---SET VECTORS FOR THE EXCITATION FRACTION AND PAIRWISE/TOTAL DIPOLAR INTERACTIONS-----
for pumpN = 1:opt.Pump_number % loop over number of pump centers
    fn = sprintf('N%d', pumpN);
    excitation_fraction.(fn)=zeros(1,length(ori1.Vecs));
    f_dd_component1.(fn)=zeros(length(spindensity_12.(fn)),length(ori1.Vecs));
    f_dd_component2.(fn)=zeros(length(spindensity_12.(fn)),length(ori1.Vecs));
    f_dd.(fn)=zeros(length(spindensity_12.(fn)),length(ori1.Vecs));
    w_dd.(fn)=zeros(length(spindensity_12.(fn)),length(ori1.Vecs));
    f_zero=0;
    ft.(fn) = zeros(size(ori1.Vecs,2),length(t));
end
F_Nt = ones(size(ori1.Vecs,2),length(t));

%---IF MAGNETOPHOTOSELECTION ON PUMP CENTER CALCULATE ORDERING WEIGHTS IN PUMP FRAMES----
if opt.Pump_magnetophotoselection ==1
    [Vecs,GeometricWeights] = sphgrid(opt.Sym,opt.nKnots,'c');
    [phi,theta] = vec2ang(Vecs);
    if opt.Pump_foldoridist == 1
        orifun = foldoridist(opt.Pump_Ordering,opt.Sym);
    else
        orifun = opt.Pump_Ordering;
    end       
    OrderingWeights = orifun(phi,theta); %define ordering weights
    % rotate ordering grid to match with pump center(s)
    for pumpN = 1:opt.Pump_number % loop over number of pump centers
        fn = sprintf('N%d', pumpN);
        p_euler = p.(fn);
        R=erot(p_euler(4:6)); %note erot generates the passive rotation matrix - we need to do an active rotation
        PumpOri.fn.Vecs = R*ori1.Vecs;
        for nn = 1:length(PumpOri.fn.Vecs)
            similarity = dot(PumpOri.fn.Vecs(:,nn).*ones(3,length(Vecs)),Vecs);
            [M,I] = max(similarity);
            OriIndex.(fn)(nn) = I;
            PumpOrderingWeights.(fn)(nn) = OrderingWeights(I); %store ordering weights
            PumpGeometricWeights.(fn)(nn) = GeometricWeights(I); %store geometric weights
        end
    end
end
%---------------------CONSTRUCT DEER SPECTRUM--------------------------
if opt.quiet == 0
    fprintf('\n3. DEER SPECTRUM CONSTRUCTION\n')
end
if exist('opt.f_dd_method')==0; opt.f_dd_method='secular'; end
if strcmp(opt.f_dd_method,'secular')==1;
    if opt.Multi_spin_time_calc == 1
        for dd=1:size(ori1.Vecs,2);
            for pumpN = 1:opt.Pump_number % loop over number of pump centers
                fn = sprintf('N%d', pumpN);
                % spectrum needed to be excited by Pump pulse at EACH orientation
                %(phi,theta) of Det. Sys. with respect to the field where (phi,theta)
                % define B_obs vector in molecular frame (phi->xy-axis,theta->z-axis)
                %-------- CONSTRUCT THE pump SPECTRUM for ((phi,theta)) --------
                y0=DEER_Sim_specmake(xaxis,En.(fn)(:,dd),Amp.(fn)(:,dd),Wid.(fn)(:,dd),opt.lineshape);
                %-------- CALCULATE THE EXCITATION FRACTION FROM OVERLAP---------
                % Overlap = Integral of {(Pump pulse excitation profile)*(theoritical spectrum needed for Pump pulse)}
                excitation_fraction.(fn)(dd)=trapz(xaxis,y0.*P_Pump);%./trapz(xaxis,y0); %normalised to a maximum of 1
                if opt.Pump_magnetophotoselection ==1
                    excitation_fraction.(fn)(dd) = excitation_fraction.(fn)(dd).*PumpOrderingWeights.(fn)(dd); % include ordering weigths
                end                
                %-------- GRAPHING THE EXCITATION FRACTION-----------------------
                if opt.Pump_graph_excitation_fraction==1, figure(8), opt.Pump_graph_excitation_fraction=graph_excitation_fraction(xaxis,y0,P_Pump,ori1.Vecs(:,dd),dd,length(excitation_fraction),opt); end
                %--------CALCULATE THE DIPOLAR CONTRIBUTIONS FOR EXCITATIONS ABOVE THE CUT OFF -------
                if excitation_fraction.(fn)(dd) > opt.Pump_excitation_fraction_cutoff
                    % calculate cos(angles) between vectors e-e (=N_12) and n_B_obs (=ori1.Vecs(:,dd))
                    cos_theta.(fn)=N_12.(fn)'*ori1.Vecs(:,dd);
                    % Calculate the dipolar coupling between the spins in Hz
                    f_dd_component1.(fn)(:,dd)=ori1.g_det(dd)*g_Pump.(fn)(dd)*k1.*spindensity_12.(fn);         % The spin.combined term is the combined spin density of the two atoms considered it acts as a weighting factor% cofficient1
                    f_dd_component2.(fn)(:,dd)=3*(cos_theta.(fn).^2)-1;                      % cofficient2, contribution of angle
                    % Summed dipolar Frequency
                    f_dd.(fn)(dd)=sum(f_dd_component1.(fn)(:,dd).*f_dd_component2.(fn)(:,dd)./(R_12.(fn)*1e-9).^3)./10^6;           % dipolar coupling between the spins in MHz
                    w_dd.(fn)(dd) = f_dd.(fn)(dd)*2*pi;   % dipolar coupling between the spins in angular frequency units.
                    
                    %calculate time dependent pairwise dipolar interaction functions
                    wt = ori1.TotalWeights(dd)./max(ori1.TotalWeights);
                    ft.(fn)(dd,:) = wt.*(1-excitation_fraction.(fn)(dd).*(1-cos(w_dd.(fn)(dd).*t)));
                end  %excitation_fraction
                %calculate total time dependent dipolar interaction term
                F_Nt(dd,:) = F_Nt(dd,:).*ft.(fn)(dd,:);
                
            end %number of pump centers
            
        end %ori1.Vecs
        
        %sum time dependent functions to average over all angles of the
        %magnetic field
        y_deer = sum(F_Nt);
        
        %perform window function and zero filling
        if opt.window == 1 && opt.zero_fill ==0
            if opt.gauss_wid ==1
                wid = gausswin(length(y_deer)*2, 2)';
            elseif opt.hann_wid ==1
                wid = hann(length(y_deer)*2)';
            elseif opt.exp_wid == 1
                wid = [exp(-4*t)' exp(-4*t)'];
            elseif opt.hamm_wid == 1
                wid = hamming(length(y_deer)*2)';
            end
            y_deer_wid = (y_deer-(ones(1,length(y_deer))*mean(y_deer(length(y_deer)/2:end)))).*wid((length(y_deer)+1):end);
            
        elseif opt.window == 1  && opt.zero_fill ==1
            if opt.gauss_wid ==1
                wid = gausswin(length(y_deer)*2, 2)';
            elseif opt.hann_wid ==1
                wid = hann(length(y_deer)*2)';
            elseif opt.exp_wid == 1
                wid = [exp(-4*t)' exp(-4*t)'];
            elseif opt.hamm_wid == 1
                wid = hamming(length(y_deer)*2)';
            end
            if opt.zero_fill_number>=length(y_deer)
                y_deer_wid = zeros(1,opt.zero_fill_number);
                y_deer_wid(1:length(y_deer)) = (y_deer-(ones(1,length(y_deer))*mean(y_deer(length(y_deer)/2:end)))).*wid((length(y_deer)+1):end);
                
            else
                y_deer_wid = (y_deer-(ones(1,length(y_deer))*mean(y_deer(length(y_deer)/2:end)))).*wid((length(y_deer)+1):end);
                
            end
        else
            y_deer_wid = y_deer;
        end
        
        % caculate frequency domain and define frequency axis
        f_sum = fftshift(fft(y_deer_wid));
        nf=1/(t(2)-t(1));                         % nquist freq.
        N=length(y_deer_wid);                            % length
        v_b=nf*(-1/2:1/N:1/2-1/N);                % freq. axis for binning
        t_wid = (0:(t(2)-t(1)):((length(y_deer_wid)-1)*(t(2)-t(1))));
        if opt.graph_DEERtrace_FFT==1, graph_DEERtrace_FFT(t_wid,y_deer_wid,v_b,f_sum,0,opt), end
    else
        y_deer = 0;
        v_b = 0;
        f_sum = 0;
        ft.N1 = 0;
        t_wid = 0;
        y_deer_wid = 0;
    end
    
    
    if opt.Pump_number == 1 && opt.Two_spin_frequency_bin == 1;
        % ----------DEFINE FREQUENCY AXIS / STARTING INTENSITY FOR BINNING---------
        nf_bin=1/(t(2)-t(1));                         % nquist freq.
        N_bin=8*length(t);                            % length
        v_b_bin=nf_bin*(-1/2:1/N_bin:1/2-1/N_bin);                % freq. axis for binning
        f_mod_bin=0*v_b_bin;                              % freq. intensities (modulated part)
        f_zero_bin=0;                                 % freq. intensities (unmodulated part)
        v_min_bin=min(v_b_bin);                           % need for freq. binning index
        y_FTdeerelta_bin=v_b_bin(2)-v_b_bin(1);               % need for freq. binning index
        %---SET VECTOR FOR THE EXCITATION FRACTION AND DIPOLAR INTERACTION-----
        excitation_fraction_bin=zeros(1,length(ori1.Vecs));
        f_dd_component1_bin=zeros(length(spindensity_12.N1),length(ori1.Vecs));
        f_dd_component2_bin=zeros(length(spindensity_12.N1),length(ori1.Vecs));
        for dd=1:size(ori1.Vecs,2);
            % spectrum needed to be excited by Pump pulse at EACH orientation
            %(phi,theta) of Det. Sys. with respect to the field where (phi,theta)
            % define B_obs vector in molecular frame (phi->xy-axis,theta->z-axis)
            %-------- CONSTRUCT THE pump SPECTRUM for ((phi,theta)) --------
            y0=DEER_Sim_specmake(xaxis,En.N1(:,dd),Amp.N1(:,dd),Wid.N1(:,dd),opt.lineshape);
            %-------- CALCULATE THE EXCITATION FRACTION FROM OVERLAP---------
            % Overlap = Integral of {(Pump pulse excitation profile)*(theoritical spectrum needed for Pump pulse)}
            excitation_fraction_bin(dd)=trapz(xaxis,y0.*P_Pump);
            if opt.Pump_magnetophotoselection ==1
                excitation_fraction_bin(dd) = excitation_fraction_bin(dd).*PumpOrderingWeights.N1(dd); % include ordering weigths
            end 
            %-------- GRAPHING THE EXCITATION FRACTION-----------------------
            if opt.Pump_graph_excitation_fraction==1, figure(8), opt.Pump_graph_excitation_fraction=graph_excitation_fraction(xaxis,y0,P_Pump,ori1.Vecs(:,dd),dd,length(excitation_fraction)); end
            %--------CALCULATE THE DIPOLAR CONTRIBUTION TO SPECTRUM FOR EXCITATIONS ABOVE THE CUT OFF --------
            if excitation_fraction_bin(dd) > opt.Pump_excitation_fraction_cutoff
                % calculate cos(angles) between vectors e-e (=N_12) and n_B_obs (=ori1.Vecs(:,dd))
                cos_theta_bin=N_12.N1'*ori1.Vecs(:,dd);
                % Calculate the dipolar coupling between the spins in Hz
                f_dd_component1_bin(:,dd)=ori1.g_det(dd).*g_Pump.N1(dd).*k1.*spindensity_12.N1;         % The spin.combined term is the combined spin density of the two atoms considered it acts as a weighting factor% cofficient1
                f_dd_component2_bin(:,dd)=3*(cos_theta_bin.^2)-1;                      % cofficient2, contribution of angle
                % Summed dipolar Frequency
                f_dd_bin=sum(f_dd_component1_bin(:,dd).*f_dd_component2_bin(:,dd)./(R_12.N1*1e-9).^3);           % dipolar coupling between the spins in Hz
                %f_dd_s(dd)=f_dd; - removed by alice jan 2011 - speed up factor of
                %10 if removed - not needed elsewhere in the program!
                % FREQUENCY BINNING
                wt=ori1.TotalWeights(dd);
                nn=round(( f_dd_bin/1e6-v_min_bin) / y_FTdeerelta_bin) +1 ;   % index for + freq. axis
                f_mod_bin(nn)=f_mod_bin(nn) + excitation_fraction_bin(dd)*wt/2;   % 2 for +ve/-ve freq.
                nn=round((-f_dd_bin/1e6-v_min_bin) / y_FTdeerelta_bin) +1 ;  % index for - freq. axis
                f_mod_bin(nn)=f_mod_bin(nn) + excitation_fraction_bin(dd)*wt/2;   % 2 for +ve/-ve freq.
                f_zero_bin=f_zero_bin+(1-excitation_fraction_bin(dd))*wt;
            else % NO excitation_fraction so add constant FREQUENCY BINNING
                f_zero_bin=f_zero_bin+ori1.TotalWeights(dd);
            end  %excitation_fraction
        end %ori1.Vecs
        f_mod_bin=f_mod_bin/sum(ori1.TotalWeights);      % normalize mod. part
        f_zero_bin=f_zero_bin/sum(ori1.TotalWeights);    % normalise unmod. part
        % form time axis from freq. bin spectrum
        f_sum_bin=f_mod_bin;                        % modulated freq. part
        f_sum_bin(N_bin/2+1)=f_sum_bin(N_bin/2+1)+f_zero_bin;   % add unmodulated (zero freq.) part
        y_FTdeer_bin=fft(fftshift((f_sum_bin)));    % fourier transform
        dt=t(2)-t(1);                       % time increments
        t_d_bin=dt*(0:length(y_FTdeer_bin)-1);      % time axis
        t_d_bin=t_d_bin(1:N_bin/8);                     % take on positive values
        y_FTdeer_bin=real(y_FTdeer_bin(1:N_bin/8));     % take on real positive values
        if opt.quiet == 0
            fprintf('...Finished. \n--- BINNING SIMULATION FINISHED! ---')
        end
        if opt.graph_DEERtrace_FFT==1, graph_DEERtrace_FFT(t_d_bin,y_FTdeer_bin,v_b_bin,f_mod_bin,f_zero_bin,opt), end
    else
        y_FTdeer_bin =0;
        t_d_bin = 0;
        v_b_bin =0;
        f_sum_bin =0;
        f_mod_bin=0;
    end
    
    if opt.quiet == 0
        fprintf('...Finished. \n--- SIMULATION FINISHED! ---')
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ---------------------------- FUNCTIONS ----------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function CWEPR_simulate(opt,P_Pump,P_det)
% Powder simulation for Detect (opt.Det_Ham) and Pump system (opt.Pump_Ham)
if opt.quiet == 0
    fprintf('\n  - POWDER SIMULATION FOR Det. and Pump Systems ...')
end
exp_powder_det=struct('Range',[min(opt.b0) max(opt.b0)], 'nPoints',length(opt.b0),'Harmonic',0, 'mwFreq',opt.Det_mw);
exp_powder_Pump=struct('Range',[min(opt.b0) max(opt.b0)], 'nPoints',length(opt.b0),'Harmonic',0, 'mwFreq',opt.Pump_mw);
if opt.Det_Ham.S>=1
exp_powder_det.Temperature=opt.Det_Triplet_pop;
end
if opt.Det_magnetophotoselection ==1
    exp_powder_det.Ordering=opt.Det_Ordering;
end
for pumpN = 1:opt.Pump_number   % loop over pump centers
    fn = sprintf('N%d', pumpN);
    opt.Pump_Ham_CW = opt.Pump_Ham.(fn);
    opt.Det_Ham_CW = opt.Det_Ham;
    if opt.Pump_Ham_CW.S>=1
    exp_powder_Pump.Temperature=opt.Pump_Triplet_pop.(fn);
    end
    if opt.Pump_magnetophotoselection ==1
    exp_powder_Pump.Ordering=opt.Pump_Ordering;
    end
    [B0,yp1] = pepper(opt.Det_Ham_CW,exp_powder_det);   % Detect System
    [B0,yp2.(fn)] = pepper(opt.Pump_Ham_CW,exp_powder_Pump); % Pump System
end
save CW_powder B0 yp1 yp2

% Energy level system
% calculate orientations
v=opt.dv+opt.Det_mw*1e3;   % unit mHz
[phi,theta] = sphgrid('Ci',10);
for ii=1:length(phi)
    En = levels(opt.Det_Ham_CW,phi(ii),theta(ii),opt.B_obs);
    E_det(ii)=En(2)-En(1);
end
w=ones(size(E_det));
a=ones(size(E_det));
fp1=DEER_Sim_specmake(v,E_det,a,w,1);   % axis MHz

for ii=1:length(phi)
    En = levels(opt.Pump_Ham_CW,phi(ii),theta(ii),opt.B_obs);
    E_pump(ii)=En(2)-En(1);
end
w=ones(size(E_pump));
a=ones(size(E_pump));
fp2=DEER_Sim_specmake(v,E_pump,a,w,1);   % axis MHz

plot(v,fp2)


function CWEPR_diagnostic(opt,P_Pump,P_det)
load CW_powder
figure(1), clf, hold on
plot(opt.B_obs*[1 1],[0.0 2],'y','linewidth',1.25),
plot(B0,P_det,'m',B0,P_Pump+1.05,'c',B0,yp1/max(yp1),'r')
for pumpN = 1:opt.Pump_number   % loop over pump centers
    fn = sprintf('N%d', pumpN);
    plot(B0,yp2.(fn)/max(yp2.(fn))+1.05,'b')
end
legend('Obs. position','Det pulse','Pump pulse',['Det Ham at ' num2str(opt.Det_mw) ' GHz'],['Pump Ham at ' num2str(opt.Pump_mw) ' GHz'])
%fprintf('... Press any key to continue'), pause

function [Pos, Amp, Wid, g]=resfields_method_EasySpin(Ham,ExpParms,Perturb)
% EasySpin to calculate resonance psostions
% Pos in mT, Amp, Wid in mT, g-values for each position
% pertubation for resonance fields
if Perturb==1, OptRes.Perturb=1; else OptRes.Perturb=0; end
[Pos,Amp,Wid]=resfields(Ham,ExpParms);
% g_Pump needed for DEER dipolar frequency
B_g=mean(Pos,1);     % centre field value corresponding to g-value position
g=71.44775*ExpParms.mwFreq./B_g; % g-value for DEER formulae


function [ori1 g_det]=det_resfields_orientations(opt,P_det)
% --- output variables ---
% ori1.TotalWeights;               % total resonance weights
% ori1.angles=[phi;theta;0*phi];   % orientations
% ori1.Vecs                        % vectors for corresponding angles
% g_det                            % g-value for each orientation
% --- Input ---
% opt.Det_Ham, opt.B_obs, ,opt.Det_mw, opt.nKnots, opt.Det_graph_ori, opt.Sym
% opt.Det_resfields_method, opt.Det_t180
% -------------------------------------------------------------
Det_Ham=opt.Det_Ham;
% calculate orientations
[Vecs,GeometricWeights] = sphgrid(opt.Sym,opt.nKnots,'c');  % triangular distribution of knots
tri = sphtri(opt.Sym,opt.nKnots);                           % triangulation info for plotting
[phi,theta] = vec2ang(Vecs);
% input parameters for freq. axis calculations
ExpDet=struct('mwFreq',opt.Det_mw,'Range',[min(opt.b0) max(opt.b0)],'Vecs',Vecs,'B_obs',opt.B_obs);
if opt.Det_magnetophotoselection ==1
    ExpDet.Ordering=opt.Det_Ordering;
end
if opt.Det_Ham.S>=1
ExpDet.Temperature=opt.Det_Triplet_pop;
end
xaxis=opt.dv+ExpDet.mwFreq*1e3;    % frequency axis [MHz]
switch(lower(opt.Det_resfields_method))
    case 'easyspin'      % resonance field positions
        ExpDet=struct('mwFreq',opt.Det_mw,'Range',[min(opt.b0) max(opt.b0)],'CrystalOrientation',[phi;theta]','Vecs',Vecs, 'Harmonic',0);
        Det_Ham_Easyspin = opt.Det_Ham;
        if opt.Det_Ham.S>=1
        ExpDet.Temperature=opt.Det_Triplet_pop;
        end
        if opt.Det_magnetophotoselection ==1
        ExpDet.Ordering=opt.Det_Ordering;
        Det_Ham_Easyspin.HStrain = mt2mhz(max(Det_Ham_Easyspin.lw)); % convert largest linewidth to Hstrain for resfields
        end
        if isfield(Det_Ham_Easyspin,'HStrain')==0 && isfield(Det_Ham_Easyspin,'gStrain')==0 && isfield(Det_Ham_Easyspin,'AStrain')==0 
             Det_Ham_Easyspin.HStrain = mt2mhz(max(Det_Ham_Easyspin.lwpp));
        end
        [En,Amp,Wid]=resfields(Det_Ham_Easyspin,ExpDet);    % all resonat positions, amplitudes, wid: as column
        xaxis=opt.b0;
        B_g=mean(En,1); g_det=1e4*0.07144775*ExpDet.mwFreq./(10*B_g); % g-value for DEER formulae
end
ResonaceWeights=zeros(length(GeometricWeights),1);
for dd=1:size(En,2)
    y0=DEER_Sim_specmake(xaxis,En(:,dd),Amp(:,dd),Wid(:,dd),opt.lineshape);
    %----------CALCULATE THE EXCITATION FRACTION FROM OVERLAP----------
    ResonaceWeights(dd)=trapz(xaxis,y0.*P_det);
    % graph overlapp integral
    if opt.Det_graph_excitation_fraction==1
        figure(2)
        opt.Det_graph_excitation_fraction=graph_excitation_fraction(xaxis,y0,P_det,ExpDet.Vecs(:,dd),dd,length(ResonaceWeights),opt);
    end
end
if opt.Det_magnetophotoselection ==1
    if opt.Det_foldoridist == 1
        orifun = foldoridist(ExpDet.Ordering,opt.Sym);
    else
        orifun = ExpDet.Ordering;
    end
       
OrderingWeights = orifun(phi,theta);
TotalWeights = ResonaceWeights(:).*GeometricWeights(:).*OrderingWeights(:); %% add in magnetoselection on detection here?

else
TotalWeights = ResonaceWeights(:).*GeometricWeights(:); 
end
TotalWeights = TotalWeights';
% Display orientational selection
if isfield(opt,'Det_graph_ori')==1
    if opt.Det_graph_ori==1
        figure(4)
        x = Vecs(1,:); y = Vecs(2,:); z = Vecs(3,:);
        trisurf(tri,x,y,z,ResonaceWeights); title('Resonance Weights')
        %colorbar
        shading interp, axis off, axis('equal')
        if isfield(opt,'gaxis')~=0
            hold on
            plot3(opt.gaxis*[0 -1],[0 0],[0 0],'k','linewidth',2)
            plot3([0 0],opt.gaxis*[0 -1],[0 0],'k','linewidth',2)
            plot3([0 0],[0 0],opt.gaxis*[0 1],'k','linewidth',2)
            hold off
        end
        
        figure(5)
        x = Vecs(1,:); y = Vecs(2,:); z = Vecs(3,:);
        trisurf(tri,x,y,z,TotalWeights); title('TotalWeights')
        %colorbar
        shading interp, axis off, axis('equal')
        if isfield(opt,'gaxis')~=0
            hold on
            plot3(opt.gaxis*[0 -1],[0 0],[0 0],'k','linewidth',2)
            plot3([0 0],opt.gaxis*[0 -1],[0 0],'k','linewidth',2)
            plot3([0 0],[0 0],opt.gaxis*[0 1],'k','linewidth',2)
            hold off
        end
        
    end
end

if isfield(opt,'Det_cutoff')==1
    % Display orientational selection cuttoff
    % this part is only for graphing
    if isfield(opt,'Det_graph_ori')==1
        if opt.Det_graph_ori==1
            figure(6)
            plot(1:length(TotalWeights),TotalWeights/max(TotalWeights),'k',1:length(ResonaceWeights),ResonaceWeights/max(ResonaceWeights),'r')
            title('totalweights - black, ResonaceWeights - red')
            
            TotalWeightsCutoff=zeros(1,length(TotalWeights));
            ResonaceWeights_normalized=ResonaceWeights/max(ResonaceWeights);
            for ii=1:length(TotalWeightsCutoff)
                if ResonaceWeights_normalized(ii)<opt.Det_cutoff
                    TotalWeightsCutoff(ii)=0;
                else
                    TotalWeightsCutoff(ii)=1;
                end
            end
            figure(7)
            x = Vecs(1,:); y = Vecs(2,:); z = Vecs(3,:);
            trisurf(tri,x,y,z,TotalWeightsCutoff); title('TotalWeightsCutoff')
            %colorbar
            shading flat
        end
        % remove all values less than Det_cutoff%
        % ResonaceWeights is not weighted by geometric weighting
        Det_cutoff=max(ResonaceWeights)*opt.Det_cutoff;
        cc=gt(ResonaceWeights,Det_cutoff);
        TotalWeights=TotalWeights(cc);
        phi=phi(cc);
        theta=theta(cc);
        Vecs=Vecs(:,cc);
        g_det=g_det(cc);
    end
end
% output variables
ori1.angles=[phi;theta];           % order for pepper
ori1.TotalWeights=TotalWeights;    % total resonance weights
ori1.Vecs=Vecs;                    % vector for corresponding angles
ori1.g_det=g_det;                  % g for corresponding angles

% save orientational related variables to harddrive
save det_orientations ori1

function graph_DEERtrace_FFT(t,y_FTdeer,v_b,f_mod,f_zero,opt)
%--------DISPLAY THE DEER TRACE AND FOURIER TRANSFORM GRAPHICALLY--------
figure(9), clf, plot(t,y_FTdeer,'g')
title(['Deer Time Domain Trace: v-det=' num2str(opt.Det_mw) ' ,v-Pump=' num2str(opt.Pump_mw) ])
xlabel('Time (micro-seconds)')
ylabel('Normalized Intensity')
figure(10), clf, plot(v_b,f_mod)
title(['Modulated Part of Frequency Domain. Zero freq. part=' num2str(f_zero)])
xlabel('frequency (MHz)')
ylabel('Normalized Intensity (unit area)')

function opt=Pump_translate_rotate(p,opt)
% rotates Pump system COORDINATES and SPIN HAMILTONIAN
% Input
% p_euler                         %Euler angles to rotate Pump coordinate system
% opt.Pump_coordinates
% opt.Pump_coordinate_centre
% opt.Pump_Ham.gFrame,opt.Pump_Ham.AFrame,opt.Pump_Ham.DFrame % rotate if they exist

% Rotation matrix from Euler angles in radians

for pumpN = 1:opt.Pump_number % loop over number of pump centers
    fn = sprintf('N%d', pumpN);
    p_euler = p.(fn);
    R=erot(p_euler(4:6)); %note erot generates the passive rotation matrix - we need to do an active rotation
    % --- Rotate around a given mean 'org' ---
    C=opt.Pump_coordinates.(fn); % coordinates
    org=opt.Pump_coordinate_centre.(fn); % centre
    for ii=1:size(C,2), C(:,ii)=C(:,ii)-org; end % translate to centre
    for ii=1:size(C,2), C(:,ii)=R*C(:,ii);  end % Rotate coordinates
    for ii=1:size(C,2), C(:,ii)=C(:,ii)+org; end % translate back
    opt.Pump_coordinates.(fn)=C; % save coordinates
    
    % --- Rotate Spin Hamiltonian
    % D-matrix axes via Euler angles ---
    if isfield(opt.Pump_Ham.(fn),'DFrame')==0, opt.Pump_Ham.(fn).DFrame=[0 0 0]; end %set to zero if not present
    DFrame=opt.Pump_Ham.(fn).DFrame;  % g-matirx Euler angles
    RD=erot(DFrame);  % g-matix axes [gx gy gz] in columns
    RD=R*RD'; % rotate axes
    RD=RD';
    [alpha,beta,gamma] = eulang(RD); % Euler angles for rotated g-axes
    opt.Pump_Ham.(fn).DFrame=[alpha beta gamma];
    
    % g-matrix axes via Euler angles ---
    if isfield(opt.Pump_Ham.(fn),'gFrame')==0, opt.Pump_Ham.(fn).gFrame=[0 0 0]; end %set to zero if not present
    gFrame=opt.Pump_Ham.(fn).gFrame;  % g-matirx Euler angles
    Rg=erot(gFrame);  % g-matix axes [gx gy gz] in columns
    Rg=R*Rg'; % rotate axes
    Rg=Rg';
    [alpha,beta,gamma] = eulang(Rg); % Euler angles for rotated g-axes
    opt.Pump_Ham.(fn).gFrame=[alpha beta gamma];
    
    % ---  A-matrix axes via Euler angles ---
    if isfield(opt.Pump_Ham.(fn),'A')==1
        if isfield(opt.Pump_Ham.(fn),'AFrame')==0; opt.Pump_Ham.(fn).AFrame=0*opt.Pump_Ham.(fn).A; end
        AFrame=opt.Pump_Ham.(fn).AFrame;
        for ii=1:size(opt.Pump_Ham.(fn).AFrame,1)
            RA=erot(AFrame(ii,:));   % hyperfine axes
            RA=R*RA';              % rotate axes
            RA=RA';
            [alpha,beta,gamma] = eulang(RA); % Euler angles for rotated A-axes
            AFrame(ii,:)=[alpha beta gamma]; % save new A-matrix orientation
        end
        opt.Pump_Ham.(fn).AFrame=AFrame; % save back into opt variable
    end
end

function output=graph_excitation_fraction(xaxis,y0,P_profile,Vecs,dd,dd_total,opt)
total_excitation_fraction=cumtrapz(xaxis,y0.*P_profile);
load CW_powder   % B0 yp1 yp2
clf, hold on

plot(xaxis,y0/max(y0)+.01,'r')  % Pump resonance for a orientation with det. resonances
plot(xaxis,P_profile+.02,'k')                    % pulse profile
plot(xaxis,total_excitation_fraction+.03,'m') %,opt.b0,total_excitation_fraction/max(total_excitation_fraction),'--m') % This value changes with the angle of the applied field changes so the position of the resonance changes

% note xaxis could be an energy scale, it has same dimension as B0 and
% for graphical purposes xaxis == B0

plot(B0,y0/max(y0)+.001,'r')  % Pump resonance for a orientation with det. resonances
plot(B0,P_profile+.002,'k')                    % pulse profile
plot(B0,total_excitation_fraction+.003,'m') %,opt.b0,total_excitation_fraction/max(total_excitation_fraction),'--m') % This value changes with the angle of the applied field changes so the position of the resonance changes
plot(B0,yp1/(max(yp1))+.004,'c')           % CW Spectrum for detection pulse
for pumpN = 1:opt.Pump_number   % loop over pump centers
    fn = sprintf('N%d', pumpN);
    plot(B0,yp2.(fn)/(max(yp2.(fn)))+.005,'b')           % CW Spectrum for Pumping pulse
end
hold off
xlabel(['B_0 [mT]'])
title([num2str(dd) ' of ' num2str(dd_total) ', B_0=[' num2str(Vecs(1)) ',' num2str(Vecs(2)), num2str(Vecs(3)) '], Excitation Fraction=' num2str(max(total_excitation_fraction))])
legend('Single-Orientation','Pulse','Excitation fraction','Det','Pump')
output=input('Graph_excitation_fraction (any key to quit)','s');
if isempty(output), output=1; else output=0; end

function [R_12 N_12 spindensity_12 opt]=Pump_det_N_spindensities(p,opt)
% output
% R_12    % distances between all spin density centre on det/Pump system
% N_12    % unit vectors between all spin density centre on det/Pump system
% spindensity_12  % spin density product p(1)*p(2) between det/Pump system spin densities
% opt.Pump_coordinates % translated corrdinates according to p_translate(phi,theta,r)

for pumpN = 1:opt.Pump_number   % loop over pump centers
    fn = sprintf('N%d', pumpN);
    p_translate = p.(fn);
    
    if opt.Print_coord==1,
        fprintf(['\nDistances & Spin Density Products'])
    end
    
    % displacement vector joining midpoint det/Pump systems,
    % phi, rotation in xy-plane in radians, theta, rotation from z in radians
    r_12  = p_translate(3);   % distance between midpoints in nm
    dR=r_12*ang2vec(p_translate(1),p_translate(2));  % [phi,theta]=[p_translate(1) p_translate(2)]
    % add vector dr=dr(phi,theta,r)
    II=ones(1,size(opt.Pump_coordinates.(fn),2));
    opt.Pump_coordinates.(fn)=opt.Pump_coordinates.(fn)+dR*II;
    % matrix of distances, unit vectors, and spindensities
    cc=1;
    for ii=1:size(opt.Det_coordinates,2)
        for jj=1:size(opt.Pump_coordinates.(fn),2)
            n=opt.Pump_coordinates.(fn)(:,jj)-opt.Det_coordinates(:,ii); % vector joint two centres
            R_12.(fn)(cc,1)=norm(n);    % distances (column)
            N_12.(fn)(:,cc)=n/norm(n);  % unit vector (in columns)
            spindensity_12.(fn)(cc,1)=opt.Det_spindensity(:,ii)*opt.Pump_spindensity.(fn)(:,jj); % spin density (in column)
            if opt.Print_coord==1,
                fprintf(['\n R=' num2str(R_12.(fn)(cc)) 'nm, p=' num2str(spindensity_12.(fn)(cc))  ])
            end
            cc=cc+1;
        end
    end
end

function y0=DEER_Sim_specmake(xaxis,En,a,w,alpha)
% y0=DEER_Sim_specmake_frequency(v,En,a,w,alpha)
% --- INPUT ---
% v      - freq. axis - must be column
% Energy - line position (MUST BE A COLUMN FOR EACH SPECTRA)
% a      - amplitude
% w      - widths, Lorentzian + Gaussian  (same units as pos)
% alpha  - 0 = loren., 1=Gaussian.
% --- OUTPUT ---
% y0    - absorption (unit area =1)
y0=zeros(length(xaxis),size(En,2));
for j=1:size(En,2)
    for i=1:size(En,1)
        % Absorption
        fwhm=w(i,j);
        xaxis0=En(i,j);
        if alpha==0
            % lorentzian
            gamma = 0.57735026918963*fwhm; % distance from x0 to inflexion point
            pre = 0.36755259694786; % 2/pi/sqrt(3)
            k = (xaxis-xaxis0)/gamma;
            y = pre/gamma./(1+4/3*k.^2);
        else
            % gaussian
            % prefactor is 1/sqrt(2*log(2))
            gamma = 0.849321800288*fwhm;   % distance from x0 to inflexion point
            pre = 0.797884560803;  % sqrt(2/pi)
            k = (xaxis-xaxis0)/gamma;
            y = pre/gamma*exp(-2*k.^2);
        end
        y0(:,j)=y0(:,j)+a(i,j)*y;
    end % i
    y0(:,j)=y0(:,j)/trapz(xaxis,y0(:,j)); % give the spectrum unit area
    %y0(:,j)=y0(:,j)/sum(a(:,j)); % give the spectrum unit area
end % j


