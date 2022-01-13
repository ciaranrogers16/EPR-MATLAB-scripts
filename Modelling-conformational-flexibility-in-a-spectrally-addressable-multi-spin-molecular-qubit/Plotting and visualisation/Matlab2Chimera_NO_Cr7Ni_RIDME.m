%% Script to update a pdb file for visualising in Chimera
%  Including XYZ arrows for plotting g-tensor orientations
%
%
% alice.bowen@manchester.ac.uk
% ciaran.rogers@manchester.ac.uk
% =========================================================================
clear all
close all

%% Read in DFT optimised single crystal XRD structure .pdb file
PDBdata = pdb2mat('Final_DFT_XRD_DA340_aligned_NO_Cr7Ni_g_aligned.pdb'); % read in data from PDB file 
%% Atomic coordinates in Gaussian frame [x,y,z], column vectors
% Geometry Optimised Coordinates
% A: Detection spin (NO)
Na=[PDBdata.X(284)    PDBdata.Y(284)   PDBdata.Z(284)];
Oa=[PDBdata.X(274)    PDBdata.Y(274)   PDBdata.Z(274)];
A=cat(1,Na,Oa);

% B: Pump spin (Cr7Ni ring)
Cr1=[PDBdata.X(1)  PDBdata.Y(1)  PDBdata.Z(1)];
Cr2=[PDBdata.X(2)  PDBdata.Y(2)  PDBdata.Z(2)];
Cr3=[PDBdata.X(3)  PDBdata.Y(3)  PDBdata.Z(3)];
Cr4=[PDBdata.X(4)  PDBdata.Y(4)  PDBdata.Z(4)];
Cr5=[PDBdata.X(5)  PDBdata.Y(5)  PDBdata.Z(5)];
Cr6=[PDBdata.X(6)  PDBdata.Y(6)  PDBdata.Z(6)];
Cr7=[PDBdata.X(7)  PDBdata.Y(7)  PDBdata.Z(7)];
Ni1=[PDBdata.X(8)  PDBdata.Y(8)  PDBdata.Z(8)];
B=cat(1,Cr1,Cr2,Cr3,Cr4,Cr5,Cr6,Cr7,Ni1);

% Read in data for whole molecule
All_atom = [PDBdata.X; PDBdata.Y; PDBdata.Z]';

% Translation of origin to the centre of the copper porphyrin - spin centres:
centre_A=Na;
for i=1:2
    A_centred(i,:)=A(i,:)-centre_A;
end
for i=1:8
    B_centred(i,:)=B(i,:)-centre_A;
end

% Translation of origin to the centre of the copper porphyrin - molecule:
for i=1:425
    All_atom_centred(i,:)=All_atom(i,:)-centre_A;
end

%% Plot centred atoms in the detection spin g-frame
figure(1)
hold on
grid on
title('Detection Spin g-frame - Spin Centres')

% Plot atoms of the detection spin A
for i=1:2
    plot3(A_centred(i,1),A_centred(i,2),A_centred(i,3),'or')
end

% Plot atoms of the pump spin B
for i=1:8
    plot3(B_centred(i,1),B_centred(i,2),B_centred(i,3),'ok')
end

axis equal
xlabel('x')
ylabel('y')
zlabel('z')

figure(2)
hold on
grid on
title('Detection Spin g-frame - Molecule')

% Plot all atoms in molecule
for i=1:425
    plot3(All_atom_centred(i,1),All_atom_centred(i,2),All_atom_centred(i,3),'ob')
end

for i=278
    plot3(All_atom_centred(i,1),All_atom_centred(i,2),All_atom_centred(i,3),'og')
end

for i=278
    plot3(All_atom_centred(i,1),All_atom_centred(i,2),All_atom_centred(i,3),'og')
end

% Plot spin centres as different colours
for i=1:2
    plot3(A_centred(i,1),A_centred(i,2),A_centred(i,3),'or')
end

for i=1:8
    plot3(B_centred(i,1),B_centred(i,2),B_centred(i,3),'ok')
end

axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%% Write updated pdb file
NAtom=pdbread('Final_DFT_XRD_DA340_aligned.pdb'); % read in pdb file

Xvals = All_atom_centred(:,1);
Yvals = All_atom_centred(:,2);
Zvals = All_atom_centred(:,3);

cellXvals = num2cell(Xvals);
cellYvals = num2cell(Yvals);
cellZvals = num2cell(Zvals);

[NAtom(:,1).Model.Atom.X] = cellXvals{:};
[NAtom(:,1).Model.Atom.Y] = cellYvals{:};
[NAtom(:,1).Model.Atom.Z] = cellZvals{:};

%pdbwrite('Final_DFT_XRD_DA340_aligned_centred_orig_130degree.pdb',NAtom);

%% Load in p-vectors
% Load in lsq fit data
load('.\lsq_fit_NO_Cr7Ni_RIDME_UPDATED_FINAL_for_fitting_all.mat');
J_store_sorted = sort(J_store)';
[J_store_count,Pval] = groupcounts(J_store_sorted);
J_matrix = [J_store_count,Pval];
display(J_matrix)
p_values = J_matrix(:,2);

fprintf('Press any key to continue'), pause

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

% %% Read in Pvals most contributing to the fit from J_store in lsq fit code
for mm=1:numel(p_values)
    % Display p value
    fprintf(['\nP-vector number ' num2str(p_values(mm)) ', Distances & Spin Density Products'])
    % Input for orientation simulation
    p.N1 = p_pumpCr7Ni_detNO(:,mm)';
    
    % Set up opt array
    N_A = [21.5556 23.7915 100.0972];
    opt.Det_Ham=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Det_spindensity=[0.5 0.5]; % spin density for each point/nuclei
    opt.Det_coordinates=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian
    opt.Det_magnetophotoselection = 0;
    
    % Pump system PARAMETERS
    opt.Pump_number = 1; %define the number of pump centers
    
    opt.Pump_coordinates.N1 = Cr7Ni_pump_overlay';
    opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
        'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
    opt.Pump_spindensity.N1=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_magnetophotoselection = 0;
    
    opt.Print_coord = 1;
    
    [opt_out]= OriDEER_multispin_coordinates_distances_spindensity(p,opt);
    
    %% Plot of detection and pump centers from oriDEER/RIDME simulation
    % p-vector should give same as manual input after rotation and translation
    
    f=figure(101);
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    m = plot3(opt_out.Det_coordinates(1,:), opt_out.Det_coordinates(2,:),...
        opt_out.Det_coordinates(3,:),'ro','MarkerSize',10);
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    det_cart_coord = opt_out.Det_coordinates;
    det_g= erot(opt_out.Det_Ham.gFrame);
    h = mArrow3(NO_pump_overlay(1,:),det_g(1,:),'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(NO_pump_overlay(1,:),det_g(2,:),'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    h = mArrow3(NO_pump_overlay(1,:),det_g(3,:),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    
    % plot pump center from transformation
    for pumpN = 1:opt.Pump_number
        fn = sprintf('N%d', pumpN);
        m = plot3(opt_out.Pump_coordinates.(fn)(1,:), opt_out.Pump_coordinates.(fn)(2,:),...
            opt_out.Pump_coordinates.(fn)(3,:),'gx','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        Pump_coordinates_centre=[sum(opt_out.Pump_coordinates.(fn)(1,:))/8;...
            sum(opt_out.Pump_coordinates.(fn)(2,:))/8;...
            sum(opt_out.Pump_coordinates.(fn)(3,:))/8]';
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
    
    %% Define variables for BILD file showing pump g-tensor
    FN = p_values(mm);
    
    gz_FN = ((pump_g(3,:)+Pump_coordinates_centre))*10; % nm --> Angstrom
    gz_centre = (Pump_coordinates_centre.*10); % nm --> Angstrom
    GZ = cat(2,gz_centre,gz_FN);
    GZ_Jmatrix = cat(2,gz_centre,J_matrix(mm,1)./5);
    
%     % Write a BILD file for best fitting pump centre. orientations - Chimera
      filename = sprintf('BILD_Cr7Ni_spheres_%d.bild',FN);
      fileID = fopen(filename,'w+');
      fprintf(fileID,'.comment -- This file shows X,Y,Z axes as red, green, blue arrows --');
      fprintf(fileID,'\n.comment -- Edit "scale" value to adjust size --');
      fprintf(fileID,'\n.comment -- p.N1 = %d',FN);
      fprintf(fileID,'\n.scale 1');
      fprintf(fileID,'\n.transparency 0.5');
      fprintf(fileID,'\n.color dark red');
      %fprintf(fileID,'\n.color 0.0 0.0 1.0');
      fprintf(fileID,'\n.sphere %2.1f %2.1f %2.1f %2.1f',GZ_Jmatrix);
      %fprintf(fileID,'\n.arrow %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f 0.25',GZ);
end

    
    
    
    
    
    
    
    
    
    
    




