%% Check best fitting p-vectors are self consistent in a three spin system
%
%
% ciaran.rogers@manchester.ac.uk
% =========================================================================
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

%% Nitroxide spin - N-O distance = + 0.136
NO_pump_overlay = [cu_centre(1,:),cu_centre(2,:),cu_centre(3,:);
    cu_centre(1,:)+0.136,cu_centre(2,:),cu_centre(3,:)];  % coordinates of NO overlayed on centre of Copper spin

%% Load in p-vectors
% Load in lsq fit data - NO det frame (relax Cr7Ni)
load('Y:\Personal Folders\CR\Matlab\02_Theoretical\oriDEER\DA340\RIDME\UPDATED_FINAL\NO_secondmax\lsq_fit_NO_Cr7Ni_RIDME_UPDATED_FINAL_for_fitting_all_max_and_3rd.mat');
J_store_sorted = sort(J_store)';
[J_store_count,Pval] = groupcounts(J_store_sorted);
J_matrix = [J_store_count,Pval];
display(J_matrix)
p_values_NO_3K = J_matrix(:,2);

fprintf('Press any key to continue'), pause

%% Rotate and translate between pump and detection frames
counter = 1;
for m=1:numel(p_values_NO_3K)
    p.N1 = p_pumpCr7Ni_detNO(:,m)';
    
    % Set up opt array
    N_A = [21.5556 23.7915 100.0972];
    opt.Det_Ham=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Det_spindensity=[0.5 0.5]; % spin density for each point/nuclei
    opt.Det_coordinates=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    % Pump system PARAMETERS
    opt.Pump_number = 1; %define the number of pump centers
    
    opt.Pump_coordinates.N1 = Cr7Ni_pump_overlay';
    opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
        'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    
    opt.Print_coord = 0;
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Center system at center of coordinates for Cr7Ni
    Pump_coordinates_centre_OPTOUT=[sum(opt_out.Pump_coordinates.N1(1,:))/8;...
        sum(opt_out.Pump_coordinates.N1(2,:))/8;...
        sum(opt_out.Pump_coordinates.N1(3,:))/8]';
        
    Cr7Ni_center_Cr7Ni(:,:,m) = opt_out.Pump_coordinates.N1(:,:)'-Pump_coordinates_centre_OPTOUT(1,:);
        
    NO_center_Cr7Ni(:,:,m) = opt.Det_coordinates(:,:)'-Pump_coordinates_centre_OPTOUT(1,:);
    
    center_Cr7Ni_center_Cr7Ni(:,m) = Pump_coordinates_centre_OPTOUT(1,:) -Pump_coordinates_centre_OPTOUT(1,:);
    center_NO_center_Cr7Ni(:,m) = [0 0 0] - Pump_coordinates_centre_OPTOUT(1,:);
    
    % Calcuate rotation matrix based on p-vector
    g_rot = erot(p_pumpCr7Ni_detNO(4:6,m)); % this is passive rotation
    
    % Rotate coordinates and centres from g frame of nitroxide
    % into g frame of ring - using an active rotation (g_rot transpose)
    
    % Coords
    Cr7Ni_center_Cr7Ni_gCr7Ni(:,:,m) = g_rot'*Cr7Ni_center_Cr7Ni(:,:,m)';
    NO_center_NO_gNO(:,:,m) = g_rot'*NO_center_Cr7Ni(:,:,m)';
    
    % Centres
    center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m) = g_rot'*center_Cr7Ni_center_Cr7Ni(:,m);
    center_NO_center_Cr7Ni_gCr7Ni(:,m) = g_rot'*center_NO_center_Cr7Ni(:,m);

    % Translate g matrices to correct position, rotate g matrices into g frame of ring and re-center
    gCr7Ni_gCr7Ni(:,:,m) = (g_rot'*(erot(p_pumpCr7Ni_detNO(4:6,m))+Pump_coordinates_centre_OPTOUT(1,:) - Pump_coordinates_centre_OPTOUT(1,:)))-center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m); % Should be [1 0 0; 0 1 0; 0 0 1] as we are now in the g-frame of the ring
    
    gNO_gCr7Ni(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]'));

    %%
    % Generate new p_vectors (back in the g-frame of ring)
    [phiCr7Ni_NO(m),thetaCr7Ni_NO(m)] = vec2ang(center_NO_center_Cr7Ni_gCr7Ni(:,m));
    g_angles_NO_Cr7Ni(:,m) = eulang(gNO_gCr7Ni(:,:,m));
    p_pumpNO_detCr7Ni_fitted(:,m) = [phiCr7Ni_NO(m),thetaCr7Ni_NO(m), norm(center_NO_center_Cr7Ni_gCr7Ni(:,m)), g_angles_NO_Cr7Ni(:,m)'];
    
    counter = counter+1;
    display(counter)
end

%% Test the fitted p_vector created (these p-vectors are for the fitted detection sequence on Cr7Ni)

% Set up opt array
opt.Det_coordinates = Cr7Ni_pump_overlay';
opt.Det_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
opt.Det_spindensity=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei

% Pump system PARAMETERS
opt.Pump_number = 1; % define the number of pump centers

% Pump system PARAMETERS (specific to center 1) - Edit for your system
N_A = [21.5556 23.7915 100.0972];
opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
    'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[0.5 0.5]; % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
opt.Pump_coordinates.N1=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

opt.Print_coord = 0;

f=figure(1);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ko');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plot=1:length(p_pumpNO_detCr7Ni_fitted)
    
    p.N1 = p_pumpNO_detCr7Ni_fitted(:,plot)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    NO_center_NO_gNO_for_plotting(:,:,plot) = NO_center_NO_gNO(:,:,plot)';
    
    % Plot coordinates
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'bo');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(NO_center_NO_gNO_for_plotting(:,1,plot),NO_center_NO_gNO_for_plotting(:,2,plot),NO_center_NO_gNO_for_plotting(:,3,plot),'bx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3([0 NO_center_NO_gNO_for_plotting(1,1,plot)], [0 NO_center_NO_gNO_for_plotting(1,2,plot)],[0 NO_center_NO_gNO_for_plotting(1,3,plot)],'b-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % Plot g axes for one center as a check
    pump_g1=erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(NO_center_NO_gNO_for_plotting(1,:,plot),pump_g1(1,:)+NO_center_NO_gNO_for_plotting(1,:,plot),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_center_NO_gNO_for_plotting(1,:,plot),pump_g1(2,:)+NO_center_NO_gNO_for_plotting(1,:,plot),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_center_NO_gNO_for_plotting(1,:,plot),pump_g1(3,:)+NO_center_NO_gNO_for_plotting(1,:,plot),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    
    %% Make .bild files for plotting in Cr7Ni detection frame
    
    FN = p_values_NO_3K(plot);
    
    Pump_coordinates_centre = (opt_out.Pump_coordinates.N1(:,2)'+opt_out.Pump_coordinates.N1(:,1)')/2;
    
    gz_FN = ((pump_g1(3,:)+Pump_coordinates_centre))*10; % nm --> Angstrom
    gz_centre = (Pump_coordinates_centre.*10); % nm --> Angstrom
    GZ = cat(2,gz_centre,gz_FN);
    GZ_Jmatrix = cat(2,gz_centre,J_matrix(plot,1)./7);
    
    % Write a BILD file for best fitting pump centre orientations - Chimera
    filename = sprintf('BILD_Cr7Nidet_NO_3K_spheres_%d.bild',FN);
    fileID = fopen(filename,'w+');
    fprintf(fileID,'.comment -- This file shows X,Y,Z axes as red, green, blue arrows --');
    fprintf(fileID,'\n.comment -- Edit "scale" value to adjust size --');
    fprintf(fileID,'\n.comment -- p.N1 = %d',FN);
    fprintf(fileID,'\n.scale 1');
    fprintf(fileID,'\n.transparency 0.5');
    fprintf(fileID,'\n.color chartreuse');
%     fprintf(fileID,'\n.color 0.0 0.0 1.0');
    fprintf(fileID,'\n.sphere %2.1f %2.1f %2.1f %2.1f',GZ_Jmatrix);
%     fprintf(fileID,'\n.arrow %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f 0.25',GZ);
end

axis equal
grid on
box on
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%%
%%
%%
%%
%% Load in lsq fit data - Cu det frame  (relax Cr7Ni)
load('.\Simulations_Cu_Cr7Ni_RIDME_UPDATED_FINAL_data_small2big.mat');
J_store_sorted = sort(J_store)';
[J_store_count,Pval] = groupcounts(J_store_sorted);
J_matrix = [J_store_count,Pval];
display(J_matrix)
p_values_Cu_3K = J_matrix(:,2);

fprintf('Press any key to continue'), pause


%% Rotate and translate between pump and detection frames
counter = 1;
for m=1:numel(p_values_Cu_3K)
    p.N1 = p_pumpCr7Ni_detCu(:,m)';
    
    % Detection system PARAMETERS - EDIT for your system
    opt.Det_coordinates = cu_det'; % coords of the detection center
    Cu_A = [-60.4642 -423.9036]; % MHz
    opt.Det_Ham=struct('S',1/2,'g',[2.0425 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]);
    opt.Det_spindensity=[0.1 0.1 0.1 0.1 0.6]; % spin density for each point/nuclei
    
    % Pump system PARAMETERS
    opt.Pump_number = 1; %define the number of pump centers
    
    opt.Pump_coordinates.N1 = Cr7Ni_pump_overlay';
    opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
        'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    
    opt.Print_coord = 0;
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Center system at center of coordinates for Cr7Ni
    Pump_coordinates_centre_OPTOUT=[sum(opt_out.Pump_coordinates.N1(1,:))/8;...
        sum(opt_out.Pump_coordinates.N1(2,:))/8;...
        sum(opt_out.Pump_coordinates.N1(3,:))/8]';
    
    Cr7Ni_center_Cr7Ni(:,:,m) = opt_out.Pump_coordinates.N1(:,:)'-Pump_coordinates_centre_OPTOUT(1,:);
    Cu_center_Cr7Ni(:,:,m) = opt.Det_coordinates(:,:)'-Pump_coordinates_centre_OPTOUT(1,:);
    
    center_Cr7Ni_center_Cr7Ni(:,m) = Pump_coordinates_centre_OPTOUT(1,:) - Pump_coordinates_centre_OPTOUT(1,:);
    center_Cu_center_Cr7Ni(:,m) = [0 0 0] - Pump_coordinates_centre_OPTOUT(1,:);
    
    % Calcuate rotation matrix based on p-vector
    g_rot = erot(p_pumpCr7Ni_detCu(4:6,m)); % this is passive rotation
    
    % Rotate coordinates and centres from g frame of nitroxide
    % into g frame of ring - using an active rotation (g_rot transpose)
    
    % Coords
    Cr7Ni_center_Cr7Ni_gCr7Ni(:,:,m) = g_rot'*Cr7Ni_center_Cr7Ni(:,:,m)';
    Cu_center_Cu_gCu(:,:,m) = g_rot'*Cu_center_Cr7Ni(:,:,m)'; %% in check plot %%
    
    % Centres
    center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m) = g_rot'*center_Cr7Ni_center_Cr7Ni(:,m);
    center_Cu_center_Cr7Ni_gCr7Ni(:,m) = g_rot'*center_Cu_center_Cr7Ni(:,m);

    % Translate g matrices to correct position, rotate g matrices into g frame of ring and re-center
    gCr7Ni_gCr7Ni(:,:,m) = (g_rot'*(erot(p_pumpCr7Ni_detCu(4:6,m))+Pump_coordinates_centre_OPTOUT(1,:) - Pump_coordinates_centre_OPTOUT(1,:)))-center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m); % Should be [1 0 0; 0 1 0; 0 0 1] as we are now in the g-frame of the ring
    
    gCu_gCr7Ni(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]'));%-center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m));

    %%
    % Generate new p_vectors (back in the g-frame of ring)
    [phiCr7Ni_Cu(m),thetaCr7Ni_Cu(m)] = vec2ang(center_Cu_center_Cr7Ni_gCr7Ni(:,m));
    g_angles_Cu_Cr7Ni(:,m) = eulang(gCu_gCr7Ni(:,:,m));
    p_pumpCu_detCr7Ni_fitted(:,m) = [phiCr7Ni_Cu(m),thetaCr7Ni_Cu(m), norm(center_Cu_center_Cr7Ni_gCr7Ni(:,m)), g_angles_Cu_Cr7Ni(:,m)'];
    
    counter = counter+1;
    display(counter)
end

%% Test the fitted p_vector created (these p-vectors are for the fitted detection sequence on Cr7Ni)

% Set up opt array
opt.Det_coordinates = Cr7Ni_pump_overlay';
opt.Det_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
opt.Det_spindensity=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei

% Pump system PARAMETERS
opt.Pump_number = 1; % define the number of pump centers

Cu_A = [-60.4642 -423.9036 -423.9036];
opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]); % if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
opt.Pump_coordinates.N1=cu_det';  % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

opt.Print_coord = 0;

f=figure(2);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ko');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plot=1:length(p_pumpCu_detCr7Ni_fitted)
    
    p.N1 = p_pumpCu_detCr7Ni_fitted(:,plot)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    Cu_center_Cu_gCu_for_plotting(:,:,plot) = Cu_center_Cu_gCu(:,:,plot)';
    
    % Plot coordinates
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ro');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(Cu_center_Cu_gCu_for_plotting(:,1,plot),Cu_center_Cu_gCu_for_plotting(:,2,plot),Cu_center_Cu_gCu_for_plotting(:,3,plot),'rx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3([0 Cu_center_Cu_gCu_for_plotting(1,1,plot)], [0 Cu_center_Cu_gCu_for_plotting(1,2,plot)],[0 Cu_center_Cu_gCu_for_plotting(1,3,plot)],'r-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % Plot g axes for one center as a check
    pump_g1=erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(Cu_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(1,:)+Cu_center_Cu_gCu_for_plotting(1,:,plot),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Cu_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(2,:)+Cu_center_Cu_gCu_for_plotting(1,:,plot),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Cu_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(3,:)+Cu_center_Cu_gCu_for_plotting(1,:,plot),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    
    %% Make .bild files for plotting in Cr7Ni detection frame
    
    FN = p_values_Cu_3K(plot);
    
    Pump_coordinates_centre = opt_out.Pump_coordinates.N1(:,1)';
    
    gz_FN = ((pump_g1(3,:)+Pump_coordinates_centre))*10; % nm --> Angstrom
    gz_centre = (Pump_coordinates_centre.*10); % nm --> Angstrom
    GZ = cat(2,gz_centre,gz_FN);
    GZ_Jmatrix = cat(2,gz_centre,J_matrix(plot,1)./8);
    
    % Write a BILD file for best fitting pump centre orientations - Chimera
    filename = sprintf('BILD_Cr7Nidet_Cu_3K_spheres_%d.bild',FN);
    fileID = fopen(filename,'w+');
    fprintf(fileID,'.comment -- This file shows X,Y,Z axes as red, green, blue arrows --');
    fprintf(fileID,'\n.comment -- Edit "scale" value to adjust size --');
    fprintf(fileID,'\n.comment -- p.N1 = %d',FN);
    fprintf(fileID,'\n.scale 1');
    fprintf(fileID,'\n.transparency 0.5');
    fprintf(fileID,'\n.color orchid');
%     fprintf(fileID,'\n.color 0.0 0.0 1.0');
    fprintf(fileID,'\n.sphere %2.1f %2.1f %2.1f %2.1f',GZ_Jmatrix);
%     fprintf(fileID,'\n.arrow %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f 0.25',GZ);
end

axis equal
grid on
box on
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';


%% Plot both Cr7Ni RIDME interactions on the same graph
for plot=1:length(p_pumpCu_detCr7Ni_fitted)
    % Display p value
    fprintf(['\nP-vector number ' num2str(plot) ', Distances & Spin Density Products'])
    
    p.N1 = p_pumpNO_detCr7Ni_fitted(:,plot)';
    p.N2 = p_pumpCu_detCr7Ni_fitted(:,plot)';
    
    % Set up opt array
    opt.Det_coordinates = Cr7Ni_pump_overlay';
    opt.Det_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
        'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
    opt.Det_spindensity=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
    
    % Pump system PARAMETERS
    opt.Pump_number = 2; % define the number of pump centers
    
    Cu_A = [-60.4642 -423.9036 -423.9036];
    opt.Pump_Ham.N2=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]); % if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N2=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N2=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N2=cu_det';  % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    
    N_A = [21.5556 23.7915 100.0972];
    opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[0.5 0.5]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N1=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    opt.Print_coord= 1;
    
    [opt_out]= OriDEER_multispin_coordinates_distances_spindensity(p,opt);
    
    %% Plot of detection and pump centers from oriDEER/RIDME simulation
    % p-vector should give same as manual input after rotation and translation
    
    f=figure(3);
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    m = plot3(opt_out.Det_coordinates(1,:), opt_out.Det_coordinates(2,:),...
        opt_out.Det_coordinates(3,:),'ko','MarkerSize',10);
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    det_cart_coord = opt_out.Det_coordinates;
    det_g= erot(opt_out.Det_Ham.gFrame);
    h = mArrow3(cu_centre(:,1),det_g(1,:),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(cu_centre(:,1),det_g(2,:),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(cu_centre(:,1),det_g(3,:),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    
    % plot pump center from transformation
    for pumpN = 1:opt.Pump_number
        fn = sprintf('N%d', pumpN);
        m = plot3(opt_out.Pump_coordinates.N1(1,:), opt_out.Pump_coordinates.N1(2,:),...
            opt_out.Pump_coordinates.N1(3,:),'bo','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        m = plot3(opt_out.Pump_coordinates.N2(1,:), opt_out.Pump_coordinates.N2(2,:),...
            opt_out.Pump_coordinates.N2(3,:),'ro','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Pump_coordinates_centre=[opt_out.Pump_coordinates.(fn)(:,1)'];
        pump_g = erot(opt_out.Pump_Ham.(fn).gFrame);
        h = mArrow3(Pump_coordinates_centre,pump_g(1,:)+Pump_coordinates_centre,...
            'color','red','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(Pump_coordinates_centre,pump_g(2,:)+Pump_coordinates_centre,...
            'color','green','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(Pump_coordinates_centre,pump_g(3,:)+Pump_coordinates_centre,...
             'color','blue','stemWidth',0.02,...
             'tipWidth',0.05,'facealpha',0.2);
    end
    
    axis equal
    box on
    grid on
    legend('g_{x}','g_{y}','g_{z}')
    
    f=make_plot_nice_p_vector(f);
    f.Renderer='opengl';

end

%% Load in lsq fit data - NO det frame (relax Cu) - 6 K
load('.\Data_Simulations_NO_Cu_RIDME_cone_model_6K.mat');
J_store_sorted = sort(J_store)';
[J_store_count,Pval] = groupcounts(J_store_sorted);
J_matrix = [J_store_count,Pval];
display(J_matrix)
p_values_NO_6K = J_matrix(:,2);

fprintf('Press any key to continue'), pause

%% Rotate and translate between pump and detection frames
counter = 1;
for m=1:numel(p_values_NO_6K)
    p.N1 = p_pumpCu_detNO(:,m)';
    
    % Set up opt array
    N_A = [21.5556 23.7915 100.0972];
    opt.Det_Ham=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Det_spindensity=[0.5 0.5]; % spin density for each point/nuclei
    opt.Det_coordinates=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    
    % Pump system PARAMETERS
    opt.Pump_number = 1; %define the number of pump centers
    
    Cu_A = [-60.4642 -423.9036 -423.9036];
    opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]); % if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N1=cu_det';  % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    opt.Print_coord = 0;
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Center system at center of coordinates for Cr7Ni
    Cu_center_Cu(:,:,m) = opt_out.Pump_coordinates.N1(:,:)'-opt_out.Pump_coordinates.N1(:,1)';
    NO_center_Cu(:,:,m) = opt.Det_coordinates(:,:)'-opt_out.Pump_coordinates.N1(:,1)';
    
    center_Cu_center_Cu(:,m) = opt_out.Pump_coordinates.N1(:,1)' - opt_out.Pump_coordinates.N1(:,1)';
    center_NO_center_Cu(:,m) = [0 0 0] - opt_out.Pump_coordinates.N1(:,1)';
    
    % Calcuate rotation matrix based on p-vector
    g_rot = erot(p_pumpCu_detNO(4:6,m)); % this is passive rotation
    
    % Rotate coordinates and centres from g frame of nitroxide
    % into g frame of ring - using an active rotation (g_rot transpose)
    
    % Coords
    Cu_center_Cu_gCu(:,:,m) = g_rot'*Cu_center_Cu(:,:,m)';
    NO_center_Cu_gCu(:,:,m) = g_rot'*NO_center_Cu(:,:,m)'; %% in check plot %%
    
    % Centres
    center_Cu_center_Cu_gCu(:,m) = g_rot'*center_Cu_center_Cu(:,m);
    center_NO_center_Cu_gCu(:,m) = g_rot'*center_NO_center_Cu(:,m);

    % Translate g matrices to correct position, rotate g matrices into g frame of ring and re-center
    gCu_gCu(:,:,m) = (g_rot'*(erot(p_pumpCu_detNO(4:6,m))+opt_out.Pump_coordinates.N1(:,1)' - opt_out.Pump_coordinates.N1(:,1)'))-center_Cu_center_Cu_gCu(:,m); % Should be [1 0 0; 0 1 0; 0 0 1] as we are now in the g-frame of the ring
    
    gNO_gCu(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]'));%-center_Cr7Ni_center_Cr7Ni_gCr7Ni(:,m));

    %%
    % Generate new p_vectors (back in the g-frame of ring)
    [phiNO_Cu(m),thetaNO_Cu(m)] = vec2ang(center_NO_center_Cu_gCu(:,m));
    g_angles_Cu_NO(:,m) = eulang(gNO_gCu(:,:,m));
    p_pumpNO_detCu_fitted(:,m) = [phiNO_Cu(m),thetaNO_Cu(m), norm(center_NO_center_Cu_gCu(:,m)), g_angles_Cu_NO(:,m)'];
    
    counter = counter+1;
    display(counter)
end

%% Test the fitted p_vector created (these p-vectors are for the fitted detection sequence on Cr7Ni)

% Detection system PARAMETERS - EDIT for your system
opt.Det_coordinates = cu_det'; % coords of the detection center
Cu_A = [-60.4642 -423.9036]; % MHz
opt.Det_Ham=struct('S',1/2,'g',[2.0425 2.1604],'gFrame',[0 0 0],...
    'HStrain',[55.7976 127.6371],'Nucs','65Cu','A',[Cu_A],...
    'lwpp',[1.0594 1.6886]);
opt.Det_spindensity=[0.1 0.1 0.1 0.1 0.6]; % spin density for each point/nuclei

% Pump system PARAMETERS
opt.Pump_number = 1; % define the number of pump centers

% Pump system PARAMETERS (specific to center 1) - Edit for your system
N_A = [21.5556 23.7915 100.0972];
opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
    'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[0.5 0.5]; % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
opt.Pump_coordinates.N1=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

opt.Print_coord = 0;

f=figure(4);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ko');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plot=1:length(p_pumpNO_detCu_fitted)
    
    p.N1 = p_pumpNO_detCu_fitted(:,plot)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    NO_center_Cu_gCu_for_plotting(:,:,plot) = NO_center_Cu_gCu(:,:,plot)';
    
    % Plot coordinates
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ro');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(NO_center_Cu_gCu_for_plotting(:,1,plot),NO_center_Cu_gCu_for_plotting(:,2,plot),NO_center_Cu_gCu_for_plotting(:,3,plot),'rx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3([0 NO_center_Cu_gCu_for_plotting(1,1,plot)], [0 NO_center_Cu_gCu_for_plotting(1,2,plot)],[0 NO_center_Cu_gCu_for_plotting(1,3,plot)],'r-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % Plot g axes for one center as a check
    pump_g1=erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(NO_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(1,:)+NO_center_Cu_gCu_for_plotting(1,:,plot),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(2,:)+NO_center_Cu_gCu_for_plotting(1,:,plot),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_center_Cu_gCu_for_plotting(1,:,plot),pump_g1(3,:)+NO_center_Cu_gCu_for_plotting(1,:,plot),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    
    %% Make .bild files for plotting in Cr7Ni detection frame
    
    FN = p_values_NO_6K(plot);
    
    Pump_coordinates_centre = opt_out.Pump_coordinates.N1(:,1)';
    
    gz_FN = ((pump_g1(3,:)+Pump_coordinates_centre))*10; % nm --> Angstrom
    gz_centre = (Pump_coordinates_centre.*10); % nm --> Angstrom
    GZ = cat(2,gz_centre,gz_FN);
    GZ_Jmatrix = cat(2,gz_centre,J_matrix(plot,1)./7);
    
    % Write a BILD file for best fitting pump centre orientations - Chimera
    filename = sprintf('BILD_Cr7Nidet_NO_6K_spheres_%d.bild',FN);
    fileID = fopen(filename,'w+');
    fprintf(fileID,'.comment -- This file shows X,Y,Z axes as red, green, blue arrows --');
    fprintf(fileID,'\n.comment -- Edit "scale" value to adjust size --');
    fprintf(fileID,'\n.comment -- p.N1 = %d',FN);
    fprintf(fileID,'\n.scale 1');
    fprintf(fileID,'\n.transparency 0.5');
    fprintf(fileID,'\n.color coral');
%     fprintf(fileID,'\n.color 0.0 0.0 1.0');
    fprintf(fileID,'\n.sphere %2.1f %2.1f %2.1f %2.1f',GZ_Jmatrix);
%     fprintf(fileID,'\n.arrow %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f 0.25',GZ);
end

axis equal
grid on
box on
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Plot both NO-ring // cu RIDME interactions on the same graph
for plot=1:length(p_pumpNO_detCu_fitted)
    % Display p value
    fprintf(['\nP-vector number ' num2str(plot) ', Distances & Spin Density Products'])
    
    p.N1 = p_pumpNO_detCu_fitted(:,plot)';
    p.N2 = p_pumpNO_detCr7Ni_fitted(:,plot)';
    
    % Set up opt array
    opt.Det_coordinates = Cr7Ni_pump_overlay';
    opt.Det_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
        'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
    opt.Det_spindensity=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
    
    % Pump system PARAMETERS
    opt.Pump_number = 2; % define the number of pump centers
    
    Cu_A = [-60.4642 -423.9036 -423.9036];
    opt.Pump_Ham.N2=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]); % if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N2=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N2=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N2=cu_det';  % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    
    N_A = [21.5556 23.7915 100.0972];
    opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[0.5 0.5]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N1=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    
    opt.Print_coord= 1;
    
    [opt_out]= OriDEER_multispin_coordinates_distances_spindensity(p,opt);
    
    %% Plot of detection and pump centers from oriDEER/RIDME simulation
    % p-vector should give same as manual input after rotation and translation
    
    f=figure(3);
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    m = plot3(opt_out.Det_coordinates(1,:), opt_out.Det_coordinates(2,:),...
        opt_out.Det_coordinates(3,:),'ko','MarkerSize',10);
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    det_cart_coord = opt_out.Det_coordinates;
    det_g= erot(opt_out.Det_Ham.gFrame);
    h = mArrow3(cu_centre(:,1),det_g(1,:),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(cu_centre(:,1),det_g(2,:),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(cu_centre(:,1),det_g(3,:),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    
    % plot pump center from transformation
    for pumpN = 1:opt.Pump_number
        fn = sprintf('N%d', pumpN);
        m = plot3(opt_out.Pump_coordinates.N1(1,:), opt_out.Pump_coordinates.N1(2,:),...
            opt_out.Pump_coordinates.N1(3,:),'bo','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        m = plot3(opt_out.Pump_coordinates.N2(1,:), opt_out.Pump_coordinates.N2(2,:),...
            opt_out.Pump_coordinates.N2(3,:),'ro','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Pump_coordinates_centre=[opt_out.Pump_coordinates.(fn)(:,1)'];
        pump_g = erot(opt_out.Pump_Ham.(fn).gFrame);
        h = mArrow3(Pump_coordinates_centre,pump_g(1,:)+Pump_coordinates_centre,...
            'color','red','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(Pump_coordinates_centre,pump_g(2,:)+Pump_coordinates_centre,...
            'color','green','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(Pump_coordinates_centre,pump_g(3,:)+Pump_coordinates_centre,...
             'color','blue','stemWidth',0.02,...
             'tipWidth',0.05,'facealpha',0.2);
    end
    
    axis equal
    box on
    grid on
    legend('g_{x}','g_{y}','g_{z}')
    
    f=make_plot_nice_p_vector(f);
    f.Renderer='opengl';

end










