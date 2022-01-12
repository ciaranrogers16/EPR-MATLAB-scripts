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
input.random = random;
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
cu_det = (([XYZ(5,3)-XYZ(5,3),XYZ(5,2)-XYZ(5,2),XYZ(5,1)-XYZ(5,1);
    XYZ(2,3)-XYZ(5,3),XYZ(2,2)-XYZ(5,2),XYZ(2,1)-XYZ(5,1);
    XYZ(3,3)-XYZ(5,3),XYZ(3,2)-XYZ(5,2),XYZ(3,1)-XYZ(5,1);
    XYZ(4,3)-XYZ(5,3),XYZ(4,2)-XYZ(5,2),XYZ(4,1)-XYZ(5,1);
    XYZ(1,3)-XYZ(5,3),XYZ(1,2)-XYZ(5,2),XYZ(1,1)-XYZ(5,1)]./0.144).*0.1996);

cu_centre = cu_det(1,:)';

%% Nitroxide spin - N-O distance = + 0.136
NO_pump_overlay = [cu_centre(1,:),cu_centre(2,:),cu_centre(3,:);
    cu_centre(1,:)+0.136,cu_centre(2,:),cu_centre(3,:)];  % coordinates of NO overlayed on centre of Copper spin

%% Define centre of conical distribution

% Vector between center of the Cu spin and NO- this will be the center of
% the conical distibution - this needs to be in the g frame of the copper
cu_NO_vec = [3.01 0 0];
input.CuNO.vec = cu_NO_vec;
% Set total number of postions to sample for each mobile center - currently
% unused but you may want to keep only a certain number of centers?
% restrict_size = 1; %logic operator to restric size of model.
% sample_no = 200; % number of orientations to keep.

%% Calcualte cones
% First cone - NO
cone_angle = 70; % set the angle of the cone to sample in degrees
grid_size = 31; %set the gride size (no of knots)
r_vary = 0.5; % set variation in r in nm
vary_percone_r = 5; % number of length variations per cone positon
angle_vary = 45; % set variation in angle to the cone vector in degrees
vary_percone_angle = 1; % number of angle variaons per cone positon

input.CuNO.cone_angle = cone_angle; % set the angle of the cone to sample in degrees
input.CuNO.grid_size = grid_size; %set the gride size (no of knots)
input.CuNO.r_vary = r_vary; % set variation in r in nm
input.CuNO.vary_percone_r = vary_percone_r; % number of length variations per cone positon
input.CuNO.angle_vary = angle_vary; % set variation in gx angle to the cone vector in degrees
input.CuNO.vary_percone_angle = vary_percone_angle; % number of angle variations per cone positon


% Model cone - NO
if random == 1
    [phi,theta] = sphrand(5000);
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

%% Orientate cone and set length

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

%% Fig 1: Plot to check cone vectors
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

Cu_g_tensor_det_frame = [1 0 0; 0 1 0 ; 0 0 1];

% Add g-z to plot using mArrow3 function
h = mArrow3(0,Cu_g_tensor_det_frame(:,3),'color','blue','stemWidth',0.05,...
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

store_alpha_rot_NO = zeros(length(x),1);
store_beta_rot_NO = zeros(length(x),1);
store_gamma_rot_NO = zeros(length(x),1);

store_phi_NO = zeros(length(x),1);
store_theta_NO = zeros(length(x),1);

store_rNO = zeros(length(x),1);

% NO
NO_angle_to_vector = 17.5; % angle of NO bond (gx) to cone vector in degrees
input.NO.angle_to_vector  = NO_angle_to_vector;

% NO
counter = 1;
vary_yz_NO = 30; % you need to edit this for the size of vatiation in the yz plane -  it is +/- this calue so total flexibility is 2* value 
for x = 1:length(rotate_vecs_NO)
    for xx = 1:size(angle_vary_all_NO,2)
        random_vector = [rand-0.5 rand-0.5 rand-0.5];
        random_unitvector = random_vector./norm(random_vector);
        % Find vector perpendicular to the cone vector and the random vector - this
        % will be an axis of rotation to rotate the NO bond away from the cone axis
        rotation_axis = cross(rotate_vecs_NO(x,:),random_unitvector);
        unit_rotation_axis = rotation_axis./norm(rotation_axis);% axis to rotate about
        unit_rotate_vecs_NO = rotate_vecs_NO(x,:)./norm(rotate_vecs_NO(x,:)); % vector to rotate
        theta = NO_angle_to_vector+angle_vary_all_NO(x,xx); % angle to rotate
        [rotatedUnitVector] = rotVecAroundArbAxis(unit_rotate_vecs_NO,unit_rotation_axis,theta);
        gx = rotatedUnitVector;
        % add in extra rotation for x-y plane
        theta_vary = (rand-0.5)*vary_yz_NO;
        [rotatedUnitVector] = rotVecAroundArbAxis(unit_rotate_vecs_NO,gx,theta_vary);
        gy = cross(gx,rotatedUnitVector);
        unit_gy = gy./norm(gy); % need each g-vector to be a unit vector in order to be used as a rotation matrix - det +-1
        gz = cross(gx,gy);
        unit_gz = gz./norm(gz);
        NO_g_tensors(:,:,counter) = [gx;unit_gy;unit_gz];
        [P_alpha_NO,P_beta_NO,P_gamma_NO]=eulang(NO_g_tensors(:,:,counter));
        % Store alpha, beta and gamma each loop
        store_alpha_NO(x) = [P_alpha_NO];
        store_beta_NO(x) = [P_beta_NO];
        store_gamma_NO(x) = [P_gamma_NO];
        % You need to translate the g tensor for NO to the postion of
        % the nitroxide then rotate then translate back if you need to
        % apply rotations in the frame of another center but in this
        % case you are already in the frame of Cu so you don't need
        % to do any thing BUT I would record the coordinate postions for
        % testing
        % Transposes of rotation matrix are active rotations into the g
        % frames - easyspin does passive rotations
        NO_rotated = (NO_g_tensors(:,:,counter)'*NO_pump_overlay')';
        NO_final_coord(:,:,counter) = NO_rotated+(rotate_vecs_NO(x,:).*ones(size(NO_rotated,1),size(NO_rotated,2)));
        r_NO_norm = norm((NO_final_coord(1,:,counter)));
        [phi_NO, theta_NO] = vec2ang(NO_final_coord(1,:,counter));
        % Store distance and anges in each loop
        store_rNO(counter) = [r_NO_norm];   
        store_phi_NO(counter) = phi_NO;
        store_theta_NO(counter) = theta_NO;
        % Generate P_vector in frame of Cr7Ni ring
        p_pumpNO_detCu(counter,:) = [store_phi_NO(counter) store_theta_NO(counter) store_rNO(counter) -store_gamma_NO(counter) -store_beta_NO(counter) -store_alpha_NO(counter)]; % this way around because active rotation - easyspin language is passive needs to be this way around (see Easyspin euler angle page)
        counter = counter+1;
        display(counter)
    end
end

%% Test the p_vector created (these p-vectors are for a detection sequence on Cu)

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

f=figure(2);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ko');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plotting = 1:10  % do 5 test plots
    plot = ceil(rand*(counter-1));
    
    p.N1(1,:) = [store_phi_NO(plot) store_theta_NO(plot) store_rNO(plot) -store_gamma_NO(plot) -store_beta_NO(plot) -store_alpha_NO(plot)];
 
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Plot coordinates
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'bo');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(NO_final_coord(:,1,plot),NO_final_coord(:,2,plot),NO_final_coord(:,3,plot),'bx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3([0 NO_final_coord(1,1,plot)], [0 NO_final_coord(1,2,plot)],[0 NO_final_coord(1,3,plot)],'b-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % Plot g axes for one center as a check
    pump_g1=erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(NO_final_coord(1,:,plot),pump_g1(1,:)+NO_final_coord(1,:,plot),'color','red','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_final_coord(1,:,plot),pump_g1(2,:)+NO_final_coord(1,:,plot),'color','green','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(NO_final_coord(1,:,plot),pump_g1(3,:)+NO_final_coord(1,:,plot),'color','blue','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector   
end
axis equal
grid on
box on
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Generate P vectors for the rotations into the frames of the NO (50 K DEER - detecting on NO)

%% Rotate into g-frame of NO
counter = 1;  
for m = 1:length(p_pumpNO_detCu)
    % Center system at center of coordinates for NO
    NO_center_NO(:,:,m) = NO_final_coord(:,:,m)-NO_final_coord(1,:,m);
    Cu_center_NO(:,:,m) = cu_det-NO_final_coord(1,:,m); 

    center_NO_center_NO(:,m) = NO_final_coord(1,:,m) - NO_final_coord(1,:,m); 
    center_Cu_center_NO(:,m) = [0 0 0] - NO_final_coord(1,:,m); 

    % Calcuate rotation matrix
    g_rot = erot(p_pumpNO_detCu(m,4:6)); % this is passive rotation
    
    % Rotate coordinates and ceters into g frame of ring 1 - using an active
    % rotation (g_rot transpose)
    NO_center_NO_gNO(:,:,m) = g_rot'*NO_center_NO(:,:,m)';
    Cu_center_NO_gNO(:,:,m) = g_rot'*Cu_center_NO(:,:,m)';

    center_NO_center_NO_gNO(:,m) = g_rot'*center_NO_center_NO(:,m);
    center_Cu_center_NO_gNO(:,m) = g_rot'*center_Cu_center_NO(:,m);
    
    % Translate g matrices to correct position, rotate g matrices into g frame of Cu and re-center
    gNO_gNO(:,:,m) = (g_rot'*(erot(p_pumpNO_detCu(m,4:6))+NO_final_coord(1,:,m)' - NO_final_coord(1,:,m)'))-center_NO_center_NO_gNO(:,m);
    gCu_gNO(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]' - NO_final_coord(1,:,m)'))-center_Cu_center_NO_gNO(:,m);
    
    % Generate p_vectors
    [phiCu_NO(m),thetaCu_NO(m)] = vec2ang(center_Cu_center_NO_gNO(:,m));
    g_angles_Cu_NO(:,m) = eulang(gCu_gNO(:,:,m));
    p_pumpCu_detNO(:,m) = [phiCu_NO(m),thetaCu_NO(m), norm(center_Cu_center_NO_gNO(:,m)), g_angles_Cu_NO(:,m)'];
    
    counter = counter+1;
    display(counter)
end

%% Test the p_vectors created.

% Set up opt array
N_A = [21.5556 23.7915 100.0972];
opt.Det_Ham=struct('S',1/2,'g',[2.0059,2.0096,2.0019],'gFrame',[0 0 0],...
    'Nucs','14N','A',N_A,'lw',1.904); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Det_spindensity=[0.5 0.5]; % spin density for each point/nuclei
opt.Det_coordinates=NO_pump_overlay'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

% Pump system PARAMETERS 
opt.Pump_number = 1; %define the number of pump centers

% Pump system PARAMETERS (specific to center 1) - Edit for your system
Cu_A = [-60.4642 -423.9036 -423.9036];
opt.Pump_Ham.N1=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
        'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
        'lwpp',[1.0594 1.6886]); % if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation
opt.Pump_coordinates.N1=cu_det';  % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

opt.Print_coord = 0;

f=figure(3);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'bo');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plotting = 1:10 % do 5 test plots
    plot1 = ceil(rand*length(p_pumpNO_detCu));

    p.N1(1,:) = p_pumpCu_detNO(:,plot1)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Plot coordinates

    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ko');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(Cu_center_NO_gNO(1,:,plot1),Cu_center_NO_gNO(2,:,plot1),Cu_center_NO_gNO(3,:,plot1),'kx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    % Plot g axes for one center as a check
    pump_g1= erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(opt_out.Pump_coordinates.N1(:,1)',pump_g1(1,:)+opt_out.Pump_coordinates.N1(:,1)','color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(opt_out.Pump_coordinates.N1(:,1)',pump_g1(2,:)+opt_out.Pump_coordinates.N1(:,1)','color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(opt_out.Pump_coordinates.N1(:,1)',pump_g1(3,:)+opt_out.Pump_coordinates.N1(:,1)','color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);%+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
end

grid on
box on
axis equal
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Save the output - three sizes of variables - with information about the cone sizes etc. in input strucutre

%save('p_vectors_only','input','p_pumpCu_detNO','p_pumpNO_detCu')


