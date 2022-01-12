function [opt_out]= OriDEER_multispin_coordinates(p,opt);
% --- rotate Pump Spin System (opt.Pump_Ham.gFrame,Apa) and coordinates(opt.Pump_coordinates) ---
opt=Pump_translate_rotate(p,opt);
% --- translate distances (R_12), unit vectors (N_12) and spin densties(spindensity_12) ---
[R_12 N_12 spindensity_12 opt]=Pump_det_N_spindensities(p,opt);
%define outputs
opt_out = opt;
end


%%% Functions %%% taken from OriLITTERSim
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
end

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
end





