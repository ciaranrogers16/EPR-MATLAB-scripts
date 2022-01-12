%% Script to build a displacement cone model and output P-vector
%
%
%%  DA340 Cu(II)-NO DEER
%   Copper detection params. from DFT EPR-II calc.
%   NO params. defined geometrically from modified single crystal
%   XRD structure
%
%   ciaran.rogers@manchester.ac.uk
%   alice.bowen@manchester.ac.uk
%==========================================================================
clear all
close all
clc

% Note this code uses the easyspin version 5.2.31 (latest version June 2021)
% - it will not work with earlier versions!

% Randomize variables?
random = 1; % 1 = randomize angles, 0 = set angles to be ca. equally spaced in a regular grid

% This model assumes that the Cu postion stays constant and the NO spin
% can be modelled as a cone comming out form the center of this spin. We
% calculate in the g frame of the copper centered at [0,0,0]

%% Overlay all spin centres in g-frame of Cr7Ni spin centred at [0, 0, 0]
% cu_det
% NO_pump_overlay

%% Import files
NAtoms=textread('CuTPPI_input_file_oriDEER.xyz','%d',1);
Atoms=cell(NAtoms,1);
XYZ=zeros(NAtoms,3); % Angstrom, cartesian coordinate system.

[Atoms,XYZ(:,1),XYZ(:,2),XYZ(:,3)]=...
    textread('CuTPPI_input_file_oriDEER.xyz','%s %f %f %f','headerlines',2);

XYZ=XYZ/10; % Angstrom -> nm

%% Copper(II) - coordinates of Cu in g frame of Cu - scaled, moved to origin
cu_det = (([XYZ(5,3),XYZ(5,2)-0.0012,XYZ(5,1)-0.0413;
    XYZ(2,3),XYZ(2,2)-0.0012,XYZ(2,1)-0.0413;
    XYZ(3,3),XYZ(3,2)-0.0012,XYZ(3,1)-0.0413;
    XYZ(4,3),XYZ(4,2)-0.0012,XYZ(4,1)-0.0413;
    XYZ(1,3),XYZ(1,2)-0.0012,XYZ(1,1)-0.0413]./0.144).*0.1996);

cu_centre = (([XYZ(5,3),XYZ(5,2)-0.0012,XYZ(5,1)-0.0413]./0.144).*0.1996)';

%% Nitroxide spin - N-O distance = +0.136
NO_pump_overlay = [cu_centre(1,:),cu_centre(2,:),cu_centre(3,:);
    cu_centre(1,:)+0.136,cu_centre(2,:),cu_centre(3,:)];  % coordinates of NO overlayed on centre of Copper spin

%% Define centre of conical distribution

% Vector between center of the Cu spin and NO- this will be the center of
% the conical distibution - this needs to be in the g frame of the copper
cu_NO_vec = [3.05 0 0];

% Set total number of postions to sample for each mobile center - currently
% unused but you may want to keep only a certain number of centers?
% restrict_size = 1; %logic operator to restric size of model.
% sample_no = 200; % number of orientations to keep.

%% Calcualte cones
% First cone - NO
cone_angle = 90; % set the angle of the cone to sample in degrees
grid_size = 31; %set the gride size (no of knots)
r_vary = 0.4; % set variation in r in nm
vary_percone_r = 5; % number of length variations per cone positon
angle_vary = 5; % set variation in angle to the cone vector in degrees
vary_percone_angle = 5; % number of angle variaons per cone positon

% Model cone - NO
if random == 1
    [phi,theta] = sphrand(1000);
    counter = 1;
    for x = 1:length(theta)
        if theta(x) < (cone_angle*(pi/180))
            vecs(counter,:) = [(ang2vec(phi(x),theta(x)))]';
            counter = counter+1;
        end
    end
else
    grid_all = sphgrid('Ci',grid_size); % note this works with the newst version of easyspin -you can edit the grid size to increase/reduce the number of outputs
    counter = 1;
    for x = 1:length(grid_all.theta)
        if grid_all.theta(x) < (cone_angle*(pi/180))
            vecs(counter,:) = [(ang2vec(grid_all.phi(x),grid_all.theta(x)))]';
            counter = counter+1;
        end
    end
end

% set variation in r
if random == 1
    counter = 1;
    for x = 1:length(vecs)
        for xx = 1:vary_percone_r
            linker_vary_all(x,xx) =(rand-0.5)*2*r_vary;
            counter = counter+1;
        end
    end
else
    linker_vary_all = (-r_vary:(r_vary*2)./(vary_percone_r-1):r_vary).*ones(length(vecs),vary_percone_r);
end

% set variation in angle
if random == 1
    counter = 1;
    for x = 1:(length(vecs)*vary_percone_r)
        for xx = 1:vary_percone_angle
            angle_vary_all(x,xx) =(rand-0.5)*2*angle_vary;
            counter = counter+1;
        end
    end
else
    angle_vary_all = (-angle_vary:(angle_vary*2)./(vary_percone_angle-1):angle_vary).*ones(length(vecs),vary_percone_r,vary_percone_angle);
end

% Store vectors and variation in r
linker_vary_all_NO = linker_vary_all;
angle_vary_all_NO = angle_vary_all;
unitvecs_NO = vecs; % these are unit vectors

%% Orientate cone and set lengths

% NO cone
[phi_NO,theta_NO] = vec2ang(cu_NO_vec); % angles of center of cone
rotation_matrix_NO = erot(phi_NO,theta_NO,0); % Rotation matrix to rotate the center of the cone
r_NO = norm(cu_NO_vec); % average distance r in nm

% Rotate center of cone to equilibrium vector.
rotate_unitvecs_NO = (rotation_matrix_NO*unitvecs_NO')';

% Set length of vectors
rotate_vecs_NO = [];
for x = 1:size(linker_vary_all_NO,2)
    rotate_vecs_NO = [rotate_vecs_NO; rotate_unitvecs_NO.*((r_NO.*ones(size(linker_vary_all_NO,1),1))+linker_vary_all_NO(:,x))];
end

% Define zero matrices to store P_vector variables
store_phi_NO = zeros(length(x),1);
store_theta_NO = zeros(length(x),1);

% Plot to check vectors
f=figure(1);
clf;
hold on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')

n = plot3(cu_det(:,1),cu_det(:,2),cu_det(:,3),'mo','MarkerSize',10); % plot ring coordinates
for x = 1:20 % Plot 2 sets of vectors for each
    % NO vector
    v1 = ceil(rand.*length(rotate_vecs_NO));
    plot3([0 rotate_vecs_NO(v1,1)],[0 rotate_vecs_NO(v1,2)],[0 rotate_vecs_NO(v1,3)],'b-')
end

% Define g-tensor of Cu ring in g-frame of Cu ring
g_cu = ([XYZ(7,3),XYZ(7,2)-0.0012,XYZ(7,1)-0.0413;
    XYZ(8,3),XYZ(8,2)-0.0012,XYZ(8,1)-0.0413;
    XYZ(9,3),XYZ(9,2)-0.0012,XYZ(9,1)-0.0413]./0.144).*0.1996;

unit_g_cu_x = g_cu(1,:)./norm(g_cu(1,:));
unit_g_cu_y = g_cu(2,:)./norm(g_cu(2,:));
unit_g_cu_z = g_cu(3,:)./norm(g_cu(3,:));

Cu_g_tensor_det_frame = [unit_g_cu_x;unit_g_cu_y;unit_g_cu_z];

% Add g-z to plot using mArrow3 function
h = mArrow3(0,unit_g_cu_z(1,:),'color','blue','stemWidth',0.05,...
    'tipWidth',0.1,'facealpha',0.2);

title('Plot of vectors relative to copper porphyrin')
legend('Cu(II)','Nitroxide')
axis equal

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Position and orientate the coordinates at the end of the vectors
% Define zero matrices to store P_vector variables
% P = [phi, theta, r, -gamma, -beta, -alpha];

store_alpha_NO = zeros(length(x),1);
store_beta_NO = zeros(length(x),1);
store_gamma_NO = zeros(length(x),1);

store_phi_NO = zeros(length(x),1);
store_theta_NO = zeros(length(x),1);

store_rNO = zeros(length(x),1);

% NO
NO_angle_to_vector = 17.5; % angle of NO bond (gx) to cone vector in degrees

% NO
counter = 1;
for x = 1:length(rotate_vecs_NO)
    for xx = 1:size(angle_vary_all_NO,2)
        random_vector = [rand rand rand];
        random_unitvector = random_vector./norm(random_vector);
        % Find vector perpendicular to the cone vector and the random vector - this
        % will be an axis of rotation to rotate the NO bond away from the cone axis
        rotation_axis = cross(rotate_vecs_NO(ceil(rand.*length(rotate_vecs_NO)),:),random_unitvector);
        unit_rotation_axis = rotation_axis./norm(rotation_axis);% axis to rotate about
        unit_rotate_vecs_NO = rotate_vecs_NO(x,:)./norm(rotate_vecs_NO(x,:)); %vector to rotate
        theta = NO_angle_to_vector+angle_vary_all_NO(x,xx); % angle to rotate
        [rotatedUnitVector,R] = rotVecAroundArbAxis(unit_rotate_vecs_NO,unit_rotation_axis,theta);
        gx = rotatedUnitVector;
        gy = cross(gx,rotate_vecs_NO(x,:));
        unit_gy = gy./norm(gy);
        gz = cross(gx,gy);
        unit_gz = gz./norm(gz);
        NO_g_tensors(:,:,counter) = [gx;unit_gy;unit_gz];
        [P_alpha_NO,P_beta_NO,P_gamma_NO]=eulang(NO_g_tensors(:,:,counter));
        % Store alpha, beta and gamma each loop
        store_alpha_NO(x) = [P_alpha_NO];
        store_beta_NO(x) = [P_beta_NO];
        store_gamma_NO(x) = [P_gamma_NO];
        NO_rotated_g_tensor(:,:,counter) = Cu_g_tensor_det_frame*NO_g_tensors(:,:,counter)';
        NO_rotated = (NO_g_tensors(:,:,counter)'*NO_pump_overlay')';
        NO_final_coord(:,:,counter) = NO_rotated+(rotate_vecs_NO(x,:).*ones(size(NO_rotated,1),size(NO_rotated,2)));
        r_NO_norm = norm((NO_final_coord(1,:,counter)));
        % Store distance each loop
        store_rNO(x) = [r_NO_norm];
        [P_phi_Cu_NO,P_theta_Cu_NO]=vec2ang(NO_final_coord(1,:,counter)); % vec2ang of centre of rotated NO coordinates
        % Store phi and theta of NO vector each loop
        store_phi_NO(x) = [P_phi_Cu_NO];
        store_theta_NO(x) = [P_theta_Cu_NO];
        counter = counter+1;
        display(counter)
    end
end

%% Plotting and graphical rendering

% Plot to check coordinates
f=figure(2);
clf;
hold on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
plot3(cu_det(:,1),cu_det(:,2),cu_det(:,3),'mo','MarkerSize',10); % Plot ring coordinates
for x = 1:20 
    v2 = ceil(rand.*length(angle_vary_all_NO));
    m=plot3([0 NO_final_coord(1,1,v2)],[0 NO_final_coord(1,2,v2)],[0 NO_final_coord(1,3,v2)],'b-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    plot3(NO_final_coord(:,1,v2),NO_final_coord(:,2,v2),NO_final_coord(:,3,v2),'bo')
    % Plot Cu g-tensor
    h = mArrow3(0,unit_g_cu_z(1,:),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
    % Plot rotated NO g-tensor
    h = mArrow3(NO_final_coord(1,:,v2),NO_rotated_g_tensor(3,:,v2)+NO_final_coord(1,:,v2),'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
end

title('Plot of rotated vectors relative to copper porphyrin')
legend('Cu(II)','Nitroxide')
axis equal

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Final P_vectors
Pvals_cell_Cu_NO_DEER=cell(length(store_phi_NO),6);
for nn=1:length(store_phi_NO)
    P_NO = [store_phi_NO(nn),store_theta_NO(nn),store_rNO(nn),store_gamma_NO(nn),store_beta_NO(nn),store_alpha_NO(nn)];
    Pvals_cell_Cu_NO_DEER{nn} = P_NO; % save P vector into 1 x nn cell array
    % Save output for main simulation code
    Pvals_DA340_Cu_NO = cell2mat(Pvals_cell_Cu_NO_DEER);
    %save('Cu_NO_DEER_DA340_Pvector_cone_model','Pvals_DA340_Cu_NO')
end

%% Test a few of the P vectors
for nn = 1:5
    % Initial input for orientation simulation
    P = [store_phi_NO(nn),store_theta_NO(nn),store_rNO(nn),store_gamma_NO(nn),store_beta_NO(nn),store_alpha_NO(nn)];
    p.N1 = P; 
    
    % Detection system PARAMETERS - EDIT for your system
    opt.Det_coordinates = cu_det'; % coords of the detection center
    Cu_A = [-60.4642 -423.9036]; % MHz
    opt.Det_Ham=struct('S',1/2,'g',[2.0425 2.1604],'gFrame',[22.5*pi/180 0 0],...
        'HStrain',[55.7976 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]);
    opt.Det_spindensity=[0.1 0.1 0.1 0.1 0.6]; % spin density for each point/nuclei
    
    % Pump system PARAMETERS (specific to center 1) - Edit for your system
    N_A =[22,19,94];
    opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
        'Nucs','14N','A',N_A,'lw',1.904);
    opt.Pump_spindensity.N1=[0.5 0.5]; % spin density for each point/nuclei
    opt.Pump_coordinate_centre.N1=[NO_pump_overlay(1,:)]'; % coords for the center of the moeity in the starting frame used as center of rotation and translation
    opt.Pump_coordinates.N1=NO_pump_overlay'; % coords of the first pump center in the starting frame used as center of rotation and translation
    
    opt.Pump_number = 1; % Number of pump centres
    opt.Print_coord=1;
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    %% Plot of detection and pump centers from oriDEER/RIDME simulation
    % p-vector should give same as manual input after rotation and translation
    
    f=figure(3);
    hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    m = plot3(opt_out.Det_coordinates(1,:), opt_out.Det_coordinates(2,:),...
        opt_out.Det_coordinates(3,:),'mo','MarkerSize',10);
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
        m = plot3(opt_out.Pump_coordinates.(fn)(1,:), opt_out.Pump_coordinates.(fn)(2,:),...
            opt_out.Pump_coordinates.(fn)(3,:),'bo','MarkerSize',10);
        m.Annotation.LegendInformation.IconDisplayStyle = 'off';
        pump_g = erot(opt_out.Pump_Ham.(fn).gFrame);
        h = mArrow3(opt_out.Pump_coordinates.(fn)(:,1)',pump_g(1,:)+opt_out.Pump_coordinates.(fn)(:,1)',...
            'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(opt_out.Pump_coordinates.(fn)(:,1)',pump_g(2,:)+opt_out.Pump_coordinates.(fn)(:,1)',...
            'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);
        h = mArrow3(opt_out.Pump_coordinates.(fn)(:,1)',pump_g(3,:)+opt_out.Pump_coordinates.(fn)(:,1)',...
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