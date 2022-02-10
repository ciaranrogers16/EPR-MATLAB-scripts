% horseradish_github:   Calculate EDNMR spectrum for arbitrary spin system
% calling the github version of EasySpin
%   Written for EasySpin 6.0-alpha1, probably works for older version, too.
%
%   [nu,spec] = horseradish(Sys,Exp)
%   ... = horseradish(Sys,Exp,Opt)
%
%   Input:
%     Sys      spin system structure
%     Exp      experiment structure
%       Exp.mwFreq       mw observer Frequency
%       Exp.Field        magnetic field, in mT
%       Exp.Range        mw offset range [min max], in MHz
%       Exp.nPoints      number of points
%       Exp.ExciteWidth  excitation bandwidth, in MHz
%       Exp.Q            Q of the cavity
%       Exp.tHTA         HTA pulse length, in us
%       Exp.nu1          nu1, in MHz, normalized for spin1/2 with g=gfree
%       Exp.Tm           decay time of electron coherences, in us
%       Exp.Temperature  Temperature in K for inital density matrix
%       Exp.CrystalOrientation Orientation of the Crystallite, can be a
%                        vector
%       Exp.CrystalWeight Weights of the CrystalOrientation vectors
%     Opt      options structure
%       Opt.Symmetry     symmetry of spin system
%       Opt.nKnots       number of knots for the orientation grid
%       Opt.Output       for Powders: 'summed': sum of isotopologues, else
%                        'separate'
%       Opt.Threshold.Probe    cutoff for transition selection
%       Opt.Threshold.Pump     cutoff if HTA has an effect on
%       Opt.Threshold.Iso      cutoff for inclusion of isotopologues
%
%
%  The simulation code here is heavily based on
%  "N.Cox, A. Nalepa, W. Lubitz, A. Savitsky, Journal of Magnetic
%  Resonance, 2017, 280, 63-78" and the deposited code on the EasySpin
%  Forum. http://easyspin.org/forum/viewtopic.php?f=8&t=404
%
%  The generation of Isotopologues and the recursive call of the function
%  was written by peeking into the code for the salt() routine of EasySpin
%
%
% Nino Wili, ETH Zurich Nov 2018, nino.wili@phys.chem.ethz.ch

function [nu,spec] = horseradish_github(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Householding and Bookkeeping                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Display help if no input is given
if (nargin==0), help(mfilename); return; end

%check input and put in default values if needed
[Sys,Exp,Opt] = ParseInputs(varargin{:});

%check whether isotopologues subspectra or all should be given in output
summedOutput = 0;
if strcmp(Opt.Output,'summed')
    summedOutput = 1;
end

%check whether the function was not called by a single isotopologues
if ~isfield(Sys,'singleiso') || ~Sys.singleiso
    
    %generate isotopologues list
    SysList = isotopologues(Sys,Opt.Threshold.Iso);
    nIsotopologues = numel(SysList);
    
    %check whether a PowderSimulation was requested
    PowderSimulation = ~isfield(Exp,'CrystalOrientation')  || isempty(Exp.CrystalOrientation);
    
    if ~PowderSimulation && ~summedOutput && nIsotopologues<2
        warning('There is no support for separate output of several crystal orientations without isotopologues. Ignored.')
    end
    
    appendSpectra = (PowderSimulation && ~summedOutput) || (~PowderSimulation && ~summedOutput && nIsotopologues>1);
    if appendSpectra
        spec = [];
    else
        spec = 0;
    end
    
    % Loop over all components and isotopologues
    for iIsotopologue = 1:nIsotopologues
        
        % Simulate single-isotopologue spectrum
        Sys_ = SysList(iIsotopologue);
        Sys_.singleiso = true;
        [nu,spec_] = horseradish_github(Sys_,Exp,Opt);
        
        % Accumulate or append spectra
        if appendSpectra
            spec = [spec; spec_*Sys_.weight];
        else
            spec = spec + spec_*Sys_.weight;
        end
        
    end
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Simulation of single isotopologue                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Orientational Grid
if ~isfield(Exp,'CrystalOrientation') || isempty(Exp.CrystalOrientation)
    %get symmetry group and symmetry frame (otherwise, the weighting will be incorrect)
    [Opt.Symmetry,Opt.SymmFrame] = symm(Sys);
    % Get orientations for the knots, molecular frame.
    [grid,tri] = sphgrid(Opt.Symmetry,Opt.nKnots(1),'cf');
    Vecs = grid.vecs;
    PowderWeight = grid.weights;
    % Transform vector to reference frame representation and convert to polar angles.
    [phi,theta] = vec2ang(Opt.SymmFrame*Vecs);
    chi = zeros(size(phi));
    Exp.AverageOverChi = 1; %integration over the third euler angle is achieved by calculating x and y only.
else
    phi = Exp.CrystalOrientation(:,1);
    theta = Exp.CrystalOrientation(:,2);
    chi = Exp.CrystalOrientation(:,3);
    Exp.AverageOverChi = 0; %the third euler average is accounted for explicitly
    if isfield(Exp,'CrystalWeight')
        PowderWeight = Exp.CrystalWeight;
    else
        PowderWeight = ones(size(Exp.CrystalOrientation,1),1);
    end
end
nOrientations = numel(PowderWeight);

% Pre-calculate spin Hamiltonian components in molecular frame
[Hzf,GxM,GyM,GzM] = sham(Sys);
zeemanop = @(v) v(1)*GxM + v(2)*GyM + v(3)*GzM;

% initialize output
%%% hey nino %%% - think about possible outputs as per orientation
[nu,spec] = makespec(Exp.Range,Exp.nPoints,Exp.Range+[eps -eps],[0 0]);

% loop over orientations
parfor iOrient = 1:nOrientations
    
    Exp_  = Exp; %local copy for the parallel loop
    Exp_ .mwFreq = Exp_ .mwFreq*1e3; %GHz to MHz
    Exp_ .CrystalOrientation = [phi(iOrient) theta(iOrient) chi(iOrient)];
    
    % calculate energy levels and probabilites between them
    [EnergyLevels, Probabilities] = SolveEigenProblem(Hzf,zeemanop,Exp_ ,Opt);
    
    % get populations
    %%% hey nino %%% - make user-supplied starting states possible
    Populations = diag(  sigeq(diag(EnergyLevels),Exp_ .Temperature)  );
    
    % transition frequencies between levels in MHz
    TransFreqs = abs(EnergyLevels-EnergyLevels');
    TransOffsets = (TransFreqs-Exp_ .mwFreq);
    
    %calculate probe weight of each transition -  ensures orientation
    %selection
    ProbeWeights = Probabilities.*exp(-2*((TransFreqs-Exp_ .mwFreq)/Exp_ .ExciteWidth).^2);
    ProbeWeights(ProbeWeights<Opt.Threshold.Probe) = 0; %apply cutoff
    
    
    if ~any(ProbeWeights) %abort this orientation if now transition is observed
        continue;
    end
    
    % Calculate W1(HTA) at EDNMR position using experimental loaded Q-value
    w1 = 2*pi*Exp_ .nu1*sqrt(1./(1+(TransOffsets*4*Exp_ .Q/(Exp_ .mwFreq)).^2/4));
    
    
    % calculate the Pump Weights using the bloch equations
    PumpWeights = (1-Bloch(w1,Exp_ .tHTA,Probabilities,Exp_ .Tm))/2; %between 0, and 1 (inversion)
    PumpWeights(PumpWeights<Opt.Threshold.Pump) = 0;  % apply cutoff
    PumpWeights(TransOffsets<Exp_ .Range(1) | TransOffsets>Exp_ .Range(2)) = 0; % only hit when within range
    
    [Amp, Pos] = EDNMR_PeakList(TransFreqs,ProbeWeights,PumpWeights,Populations);
    
    %add up the spectrum from this powder orientation
    if ~isempty(Pos)
        %get rid of peaks slightly outside range
        %cannot really be done earlier, because it needs the knowledge of connected transitions)
        Amp(Pos<Exp_ .Range(1) | Pos>Exp_ .Range(2)) = [];
        Pos(Pos<Exp_ .Range(1) | Pos>Exp_ .Range(2)) = [];
        spec = spec+PowderWeight(iOrient)*makespec(Exp_ .Range,Exp_ .nPoints,Pos,Amp);
    end
    
end %end of powder loop


%build final spectrum
dnu = nu(2)-nu(1);
spec = convspec(spec,dnu,Sys.lwEndor);

end % end of main function

function [Sys,Exp,Opt] = ParseInputs(varargin)

% Guard against wrong number of input or output arguments.
if (nargin<1), error('Please supply a spin system as first parameter.'); end
if (nargin<2), error('Please supply experimental parameters as second input argument.'); end
if (nargin>3), error('Too many input arguments, the maximum is three.'); end

if (nargout>4), error('Too many output arguments.'); end

% Put vargin into structures, check, and add default values
Sys = varargin{1};
Exp = varargin{2};
%%% hey nino %%% - you should add some more Exp-checks here
if ~isfield(Exp,'Temperature') Exp.Temperature = 300; end

if nargin==3
    Opt = varargin{3};
else
    Opt = struct;
end
if ~isfield(Opt,'nKnots') Opt.nKnots = 31; end
if ~isfield(Opt,'Threshold')
    Opt.Threshold = struct;
end
if ~isfield(Opt.Threshold,'Probe')  Opt.Threshold.Probe = 1e-4; end
if ~isfield(Opt.Threshold,'Pump')  Opt.Threshold.Pump = 1e-4; end
if ~isfield(Opt.Threshold,'Iso')  Opt.Threshold.Iso = 1e-3; end
if ~isfield(Opt,'Output')
    Opt.Output = 'summed';
else
    if ~strcmp(Opt.Output,'summed') && ~strcmp(Opt.Output,'separate')
        error('Please specify Opt.Output as summed or separate ')
    end
end

end

function [EnergyLevels, Probabilities] = SolveEigenProblem(Hzf,zeemanop,Exp,Opt)

% Set up lab axes for this orientation
[xL,yL,zL] = erot(Exp.CrystalOrientation,'rows');
% zL = z laboratoy axis: external static field
% xL = x laboratory axis: B1 excitation field
% yL = y laboratory vector: needed for intensity integration over third Euler angle

% Set up spin Hamiltonian and calculate energy levels and transition
% probabilities (normalized to 1 for a spin-1/2 with g=gfree)
H0 = Hzf + Exp.Field*zeemanop(zL);
[U,Hd] = eig(H0);
EnergyLevels = diag(Hd); %Energy levels

gamma = gfree*bmagn/planck*1e-9; %MHz/mT
if Exp.AverageOverChi %average is done ad hoc over x and y only
    Probabilities = 2/gamma^2*(abs(U'*zeemanop(xL)*U).^2+abs(U'*zeemanop(yL)*U).^2);
else
    Probabilities = 2*2/gamma^2*(abs(U'*zeemanop(xL)*U).^2);
end

end

function [Mz] = Bloch(w1,t,P,Tm)

%for now assuming  [0 0 1] as startingvector
% analytical solution of Mz in the presence of transverse relaxation

fac = exp(-t/(2*Tm));
Ch = cosh(t*sqrt(1-4*P.*Tm.^2.*w1.^2)/(2*Tm));
Sh = sinh(t*sqrt(1-4*P.*Tm.^2.*w1.^2)/(2*Tm));
Mz = fac.*(Ch + Sh./(sqrt(1-4*P.*Tm.^2.*w1.^2)));

end

function [Amp, Pos] = EDNMR_PeakList(TransFreqs,ProbeWeights,PumpWeights,Populations)

%initalize output
Pos = [];
Amp = [];

%abbreviations
Pops = Populations;

%loop over all probe transition
for i_probe = 1:size(TransFreqs,1)
    for j_probe = (i_probe+1):size(TransFreqs,2)
        
        if ProbeWeights(i_probe,j_probe)
            
            %loop over all pump transition
            for i_pump = 1:size(TransFreqs,1)
                for j_pump = (i_pump+1):size(TransFreqs,2)
                    
                    %only transitions connected to the probe and affected
                    %by HTA are considered
                    % the first condition is fulfilled if and only if
                    % exactly one level is shared
                    if numel(unique([i_probe j_probe i_pump j_pump]))==3 && PumpWeights(i_pump,j_pump)
                        
                        Pos = [Pos TransFreqs(i_pump,j_pump)-TransFreqs(i_probe,j_probe)]; %postion of pumped transition wrt to observed one(not! the resonator.)
                        
                        %calculate populations after the pump pulse
                        PumpPops = Pops;
                        PumpPops(i_pump) = (1-PumpWeights(i_pump,j_pump))*Pops(i_pump)+PumpWeights(i_pump,j_pump)*Pops(j_pump);
                        PumpPops(j_pump) = (1-PumpWeights(i_pump,j_pump))*Pops(j_pump)+PumpWeights(i_pump,j_pump)*Pops(i_pump);
                        
                        %calculate polarization difference of the probed
                        %transition
                        PolDif = (PumpPops(j_probe)-PumpPops(i_probe))-(Pops(j_probe)-Pops(i_probe));
                        
                        %calculate peak amplitude and append to list
                        Amp = [Amp PolDif*ProbeWeights(i_probe,j_probe)];
                    end
                    
                end
            end
            
        end
        
    end
end

end %end of peak function




% Denä wos guät geit, 
% giängs besser, 
% giängs denä besser, 
% wos weniger guät geit, 
% was aber nid geit 
% ohni dases denä 
% weniger guät geit, 
% wos guät get.

% Drum geit weni, 
% für dases denä 
% besser geit, 
% wos weniger guät geit. 
% Und drum geits o denä nid besser, 
% wos guät geit.
% 
% Mani Matter - 1970