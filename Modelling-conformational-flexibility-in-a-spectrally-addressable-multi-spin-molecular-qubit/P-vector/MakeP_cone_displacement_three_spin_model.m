%% Script to build a displacement cone model and P-vector libraries
%
%
%%  DA340 RIDME
%   Copper detection params. from DFT EPR-II calc.
%   Cr7Ni and NO params. defined geometrically from modified single crystal
%   XRD structure DA333
%
%   ciaran.rogers@manchester.ac.uk
%   alice.bowen@manchester.ac.uk
%==========================================================================
clear all
close all
clc

%% Notes on P vectors generated

% As the P vectors in the g frame of the ring are generated independently
% there are two variables associated with these vectors x and xx
% These vary in the following way: 
% p_pumpNO_detCr7Ni(x,:)
% p_pumpCu_detCr7Ni(xx,:)

% In the g frame of the NO these vary: 
% p_pumpCr7Ni_detNO(:,x) 
% p_pumpCu_detNO(:,x,xx) - this has a double loop as I've made it such that
% each NO-ring vector also has all possible ring-Cu vectors. 

% In the g frame of the Cu these vary: 
% p_pumpCr7Ni_detCu(:,xx) 
% p_pumpNO_detCu(:,xx,x) - this has a double loop as I've made if such that
% each NO-ring also has all possible ring Cu vectors. 

%% 
% Note this code uses the easyspin version 5.2.31 (latest version June 2021)
% - it will not work with earlier versions! (Only true for grid size option 
% in non-random angles....)

% Randomize variables?
random = 1; % 1 = randomize angles, 0 = set angles to be ca. equally spaced in a regular grid
input.random = random;
% This model assumes that the ring postion stays constant and the Cu and NO
% can be modelled as cones comming out form the center of this ring. We
% calculate in the g frame of the ring with the ring centered at 000 - the
% nitroxide is above the plane of the ring and the Cu is below

%% Overlay all spin centres in g-frame of Cr7Ni spin centred at [0, 0, 0]
% cu_det
% NO_pump_overlay
% Cr7Ni_pump_overlay

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

%% Define centre of conical distribution
% Nitroxide is above the plane of the ring and the Cu is below

% Vector between center of the ring and NO- this will be the center of the
% conical distibution - this needs to be in the g frame of the ring
ring_NO_vec = [0 .4 1.50];
input.NO.vec = ring_NO_vec;
% Vector between center of the ring and Cu- this will be the center of the
% conical distibution - this needs to be in the g frame of the ring
ring_Cu_vec = [0 2 -1.75];
input.Cu.vec = ring_Cu_vec;

% Set total number of postions to sample for each mobile center - currently
% unused but you may want to keep only a certain number of centers?
% restrict_size = 1; %logic operator to restrict size of model.
% sample_no = 200; %number of orientations to keep.

%% Calcualte cones
% First cone - NO
cone_angle = 34; % set the angle of the cone to sample in degrees
grid_size = 31; %set the gride size (no of knots)
r_vary = 0.5; % set variation in r in nm
vary_percone_r = 1; % number of length variations per cone positon
angle_vary = 1; % set variation in gx angle to the cone vector in degrees
vary_percone_angle = 1; % number of angle variations per cone positon

input.NO.cone_angle = cone_angle; % set the angle of the cone to sample in degrees
input.NO.grid_size = grid_size; %set the gride size (no of knots)
input.NO.r_vary = r_vary; % set variation in r in nm
input.NO.vary_percone_r = vary_percone_r; % number of length variations per cone positon
input.NO.angle_vary = angle_vary; % set variation in gx angle to the cone vector in degrees
input.NO.vary_percone_angle = vary_percone_angle; % number of angle variations per cone positon

% Model cone - NO
if random == 1
    [phi,theta] = sphrand(5000);
    counter = 1;
    for x = 1:length(theta)
        if theta(x) <= (cone_angle*(pi/180))
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

% Second cone - Cu
cone_angle = 34; % set the angle of the cone to sample
grid_size = 31; % set the gride size (no of knots)
r_vary = 0.5;% set variation in r in nm
vary_percone_r = 1; % number of length variaons per cone positon
angle_vary = 1; % set variation in gz angle to the cone vector in degrees
vary_percone_angle = 1; % number of angle variations per cone positon

input.Cu.cone_angle = cone_angle; % set the angle of the cone to sample in degrees
input.Cu.grid_size = grid_size; %set the gride size (no of knots)
input.Cu.r_vary = r_vary; % set variation in r in nm
input.Cu.vary_percone_r = vary_percone_r; % number of length variations per cone positon
input.Cu.angle_vary = angle_vary; % set variation in gx angle to the cone vector in degrees
input.Cu.vary_percone_angle = vary_percone_angle; % number of angle variations per cone positon

% Model cone - Cu
if random == 1
    [phi,theta] = sphrand(5000);
    counter = 1;
    for x = 1:length(theta)
        if theta(x) <= (cone_angle*(pi/180))
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


% store vectors and variation in r
linker_vary_all_Cu = linker_vary_all;
angle_vary_all_Cu = angle_vary_all;
unitvecs_Cu = vecs; % these are unit vectors


%% Orientate cones and set lengths

% NO cone
[phi_NO,theta_NO] = vec2ang(ring_NO_vec); % angles of center of cone
rotation_matrix_NO = erot(phi_NO,theta_NO,0); % Rotation matrix to rotate the center of the cone
r_NO = norm(ring_NO_vec); % average distance r in nm

% Cu cone
[phi_Cu,theta_Cu] = vec2ang(ring_Cu_vec); % angles of center of cone
rotation_matrix_Cu = erot(phi_Cu,theta_Cu,0); % rotation matrix to rotate the center of the cone
r_Cu = norm(ring_Cu_vec); % average distance r in nm

% Rotate center of cone to equilibrium vector.
rotate_unitvecs_NO = (rotation_matrix_NO*unitvecs_NO')';
rotate_unitvecs_Cu = (rotation_matrix_Cu*unitvecs_Cu')';

% Set length of vectors
rotate_vecs_Cu = [];
for x = 1:size(linker_vary_all_Cu,2)
    rotate_vecs_Cu = [rotate_vecs_Cu; rotate_unitvecs_Cu.*((r_Cu.*ones(size(linker_vary_all_Cu,1),1))+linker_vary_all_Cu(:,x))];
end

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

n = plot3(Cr7Ni_pump_overlay(:,1),Cr7Ni_pump_overlay(:,2),Cr7Ni_pump_overlay(:,3),'ko','MarkerSize',10); % plot ring coordinates

for x = 1:1000
    % Cu vector
    v1 = ceil(rand.*length(rotate_vecs_Cu));
    plot3([0 rotate_vecs_Cu(v1,1)],[0 rotate_vecs_Cu(v1,2)],[0 rotate_vecs_Cu(v1,3)],'r-') 
    % NO vector
    v2 = ceil(rand.*length(rotate_vecs_NO));
    plot3([0 rotate_vecs_NO(v2,1)],[0 rotate_vecs_NO(v2,2)],[0 rotate_vecs_NO(v2,3)],'b-')
end

%% Define g-tensor of Cr7Ni ring in g-frame of Cr7Ni ring
Cr7Ni_pump_overlay_norm = Cr7Ni_pump_overlay(1,:)./norm(Cr7Ni_pump_overlay(1,:));

ring_NO_vec_cr7ni_gframe = [0 0 1.54];
Initial_Cr7Ni_gz = ring_NO_vec_cr7ni_gframe./norm(ring_NO_vec_cr7ni_gframe);
Initial_Cr7Ni_gy = cross(Initial_Cr7Ni_gz,Cr7Ni_pump_overlay_norm(1,:));
Initial_Cr7Ni_gx = cross(Initial_Cr7Ni_gz,Initial_Cr7Ni_gy);
%%%% editied to make the determient of this g tensor be +1 not -1 as this
%%%% is making things difficult I think.
Cr7Ni_g_tensor = [-Initial_Cr7Ni_gx;Initial_Cr7Ni_gy;Initial_Cr7Ni_gz];

% Add g-z to plot using mArrow3 function
h = mArrow3(0,Initial_Cr7Ni_gx(1,:),'color','red','stemWidth',0.05,...
    'tipWidth',0.1,'facealpha',0.2);
h = mArrow3(0,Initial_Cr7Ni_gy(1,:),'color','green','stemWidth',0.05,...
    'tipWidth',0.1,'facealpha',0.2);
h = mArrow3(0,Initial_Cr7Ni_gz(1,:),'color','blue','stemWidth',0.05,...
    'tipWidth',0.1,'facealpha',0.2);

title('Plot of vectors relative to Cr_{7}Ni ring')
legend('Cr_{7}Ni ring','Cu(II)','Nitroxide')
axis equal

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Position and orientate the coordinates at the end of the vectors
% Define zero matrices to store P_vector variables
% P = [phi, theta, r, -gamma, -beta, -alpha]; % active --> passive

% Cr7Ni g-frame
store_alpha_cu = zeros(length(x),1);
store_alpha_NO = zeros(length(x),1);
store_beta_cu = zeros(length(x),1);
store_beta_NO = zeros(length(x),1);
store_gamma_cu = zeros(length(x),1);
store_gamma_NO = zeros(length(x),1);

% Detection spin g-frame
store_phi_cu = zeros(length(x),1);
store_theta_cu = zeros(length(x),1);
store_phi_NO = zeros(length(x),1);
store_theta_NO = zeros(length(x),1);

store_rNO = zeros(length(x),1);
store_rCu = zeros(length(x),1);

store_alpha_rot_cu = zeros(length(x),1);
store_alpha_rot_NO = zeros(length(x),1);
store_beta_rot_cu = zeros(length(x),1);
store_beta_rot_NO = zeros(length(x),1);
store_gamma_rot_cu = zeros(length(x),1);
store_gamma_rot_NO = zeros(length(x),1);

% NO
NO_angle_to_vector = 17.5; % angle of NO bond (gx) to cone vector in degrees
input.NO.angle_to_vector  = NO_angle_to_vector;
% Cu
Cu_angle_to_vector = 90; % angle of Cu gz to cone vector in degrees
input.Cu.angle_to_vector  = Cu_angle_to_vector;

%% NO - generate p vector in inital Cr7Ni frame
counter = 1;
vary_yz_NO = 15; % you need to edit this for the size of vatiation in the yz plane -  it is +/- this calue so total flexibility is 2* value 
for x = 1:length(rotate_vecs_NO)
    for xx = 1:size(angle_vary_all_NO,2)
        %% NO spin as pump i.e. in g-frame of Cr7Ni
        % rand only gives values in 0 to +1 range to get a full range of
        % values you need to have the posibility to have -ve numbers too
        % in the vector. 
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
        [P_alpha_NO,P_beta_NO,P_gamma_NO]=eulang(NO_g_tensors(:,:,counter)); % active rotation matrix
        % Store alpha, beta and gamma each loop
        store_alpha_NO(counter) = [P_alpha_NO];
        store_beta_NO(counter) = [P_beta_NO];
        store_gamma_NO(counter) = [P_gamma_NO]; 
        % You need to translate the g tensor for NO to the postion of
        % the nitroxide then rotate then translate back if you need to
        % apply rotations in the frame of another center but in this
        % case you are already in the frame of Cr7Ni so you don't need
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
        p_pumpNO_detCr7Ni(counter,:) = [store_phi_NO(counter) store_theta_NO(counter) store_rNO(counter) -store_gamma_NO(counter) -store_beta_NO(counter) -store_alpha_NO(counter)]; % this way around because active rotation - easyspin language is passive needs to be this way around (see Easyspin euler angle page)
        % We will rotate both sets of coordinates into the other frames later at
        % the same time
        counter = counter+1;
        display(counter)
    end   
end

%% Test the p_vector created (these p-vectors are for a detection sequence on Cr7Ni)

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

%% Cu
counter = 1;
vary_xy_Cu = 20; % you need to edit this for the size of vatiation in the xy plane -  it is +/- this calue so total flexibility is 2* value 
for x = 1:length(rotate_vecs_Cu)
    for xx = 1:size(angle_vary_all_Cu,2)        
        random_vector = [rand-0.5 rand-0.5 rand-0.5];
        random_unitvector = random_vector./norm(random_vector);
        rotation_axis = cross(rotate_vecs_Cu(x,:),random_unitvector);
        unit_rotation_axis = rotation_axis./norm(rotation_axis);% axis to rotate about
        unit_rotate_vecs_Cu = rotate_vecs_Cu(x,:)./norm(rotate_vecs_Cu(x,:)); %vector to rotate
        theta = Cu_angle_to_vector+angle_vary_all_Cu(x,xx); % angle to rotate
        [rotatedUnitVector] = rotVecAroundArbAxis(unit_rotate_vecs_Cu,unit_rotation_axis,theta);      
        gz = rotatedUnitVector;
        % add in extra rotation for x-y plane
        theta_vary = (rand-0.5)*vary_xy_Cu;
        [rotatedUnitVector] = rotVecAroundArbAxis(unit_rotate_vecs_Cu,gz,theta_vary);
        gy = cross(gz,rotatedUnitVector);
        unit_gy = gy./norm(gy);
        gx = cross(gz,gy);
        unit_gx = gx./norm(gx);
        Cu_g_tensors(:,:,counter) = -[unit_gx;unit_gy;gz];
        [P_alpha_Cu,P_beta_Cu,P_gamma_Cu]=eulang(Cu_g_tensors(:,:,counter));
        % Store alpha, beta and gamma each loop
        store_alpha_Cu(counter) = [P_alpha_Cu];
        store_beta_Cu(counter) = [P_beta_Cu];
        store_gamma_Cu(counter) = [P_gamma_Cu];
        Cu_rotated = (Cu_g_tensors(:,:,counter)'*cu_det')';
        Cu_final_coord(:,:,counter) = Cu_rotated+(rotate_vecs_Cu(x,:).*ones(size(Cu_rotated,1),size(Cu_rotated,2)));
        r_Cu_norm = norm((Cu_final_coord(1,:,counter)));
        [phi_Cu, theta_Cu] = vec2ang(Cu_final_coord(1,:,counter));
        % Store distance and anges in each loop
        store_rCu(counter) = [r_Cu_norm];
        store_phi_Cu(counter) = phi_Cu;
        store_theta_Cu(counter) = theta_Cu;
        % Generate P_vector
        p_pumpCu_detCr7Ni(counter,:) = [store_phi_Cu(counter) store_theta_Cu(counter) store_rCu(counter) -store_gamma_Cu(counter) -store_beta_Cu(counter) -store_alpha_Cu(counter)];
        counter = counter+1;
        display(counter)
    end
end

%% Test the p_vector created (these p-vectors are for a detection sequence on Cr7Ni)

% Set up opt array
opt.Det_coordinates = Cr7Ni_pump_overlay';
opt.Det_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]);
opt.Det_spindensity=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei

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
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ko');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';
for plotting = 1:10  % do 5 test plots
    plot = ceil(rand*(counter-1));
    
    p.N1(1,:) = [store_phi_Cu(plot) store_theta_Cu(plot) store_rCu(plot) -store_gamma_Cu(plot) -store_beta_Cu(plot) -store_alpha_Cu(plot)];
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Plot coordinates
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ro');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(Cu_final_coord(:,1,plot),Cu_final_coord(:,2,plot),Cu_final_coord(:,3,plot),'rx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3([0 Cu_final_coord(1,1,plot)], [0 Cu_final_coord(1,2,plot)],[0 Cu_final_coord(1,3,plot)],'r-');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
   
    % Plot g axes for one center as a check
    pump_g1= erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(Cu_final_coord(1,:,plot),pump_g1(1,:)+Cu_final_coord(1,:,plot),'color','red','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Cu_final_coord(1,:,plot),pump_g1(2,:)+Cu_final_coord(1,:,plot),'color','green','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Cu_final_coord(1,:,plot),pump_g1(3,:)+Cu_final_coord(1,:,plot),'color','blue','stemWidth',0.02,...
            'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector      
end
axis equal
grid on
box on
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';
%% Generate P vectors for the rotations into the frames of the Cu and NO

%% Rotate into g-frame of NO
counter = 1;  
for m = 1:length(p_pumpNO_detCr7Ni)
    % Center system at center of coordinates for NO
    NO_center_NO(:,:,m) = NO_final_coord(:,:,m)-NO_final_coord(1,:,m);
    Cr7Ni_center_NO(:,:,m) = Cr7Ni_pump_overlay-NO_final_coord(1,:,m); 

    center_NO_center_NO(:,m) = NO_final_coord(1,:,m) - NO_final_coord(1,:,m); 
    center_Cr7Ni_center_NO(:,m) = [0 0 0] - NO_final_coord(1,:,m); 

    % Calcuate rotation matrix
    g_rot = erot(p_pumpNO_detCr7Ni(m,4:6)); % this is passive rotation
    
    % Rotate coordinates and ceters into g frame of ring 1 - using an active
    % rotation (g_rot transpose)
    NO_center_NO_gNO(:,:,m) = g_rot'*NO_center_NO(:,:,m)';
    Cr7Ni_center_NO_gNO(:,:,m) = g_rot'*Cr7Ni_center_NO(:,:,m)';

    center_NO_center_NO_gNO(:,m) = g_rot'*center_NO_center_NO(:,m);
    center_Cr7Ni_center_NO_gNO(:,m) = g_rot'*center_Cr7Ni_center_NO(:,m);
    
    % Translate g matrices to correct position, rotate g matrices into g frame of ring 1 and re-center
    gNO_gNO(:,:,m) = (g_rot'*(erot(p_pumpNO_detCr7Ni(m,4:6))+NO_final_coord(1,:,m)' - NO_final_coord(1,:,m)'))-center_NO_center_NO_gNO(:,m);
    gCr7Ni_gNO(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]' - NO_final_coord(1,:,m)'))-center_Cr7Ni_center_NO_gNO(:,m);
    
    % Generate p_vectors
    [phiCr7Ni_NO(m),thetaCr7Ni_NO(m)] = vec2ang(center_Cr7Ni_center_NO_gNO(:,m));
    g_angles_Cr7Ni_NO(:,m) = eulang(gCr7Ni_gNO(:,:,m));
    p_pumpCr7Ni_detNO(:,m) = [phiCr7Ni_NO(m),thetaCr7Ni_NO(m), norm(center_Cr7Ni_center_NO_gNO(:,m)), g_angles_Cr7Ni_NO(:,m)'];
    
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

opt.Pump_coordinates.N1 = Cr7Ni_pump_overlay';
opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation

opt.Print_coord = 0;

f=figure(4);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'bo');
m.Annotation.LegendInformation.IconDisplayStyle = 'off';

for plotting = 1:5  % do 5 test plots
    plot1 = ceil(rand*length(p_pumpNO_detCr7Ni));

    p.N1(1,:) = p_pumpCr7Ni_detNO(:,plot1)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Plot coordinates

    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ko');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(Cr7Ni_center_NO_gNO(1,:,plot1),Cr7Ni_center_NO_gNO(2,:,plot1),Cr7Ni_center_NO_gNO(3,:,plot1),'kx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    Pump_coordinates_centre=[sum(opt_out.Pump_coordinates.N1(1,:))/8;...
        sum(opt_out.Pump_coordinates.N1(2,:))/8;...
        sum(opt_out.Pump_coordinates.N1(3,:))/8]';
    % Plot g axes for one center as a check
    pump_g1= erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(Pump_coordinates_centre,pump_g1(1,:)+Pump_coordinates_centre,'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Pump_coordinates_centre,pump_g1(2,:)+Pump_coordinates_centre,'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Pump_coordinates_centre,pump_g1(3,:)+Pump_coordinates_centre,'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2);%+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
end

grid on
box on
axis equal
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Rotate into g-frame of Cu
counter = 1;
for m = 1:length(p_pumpCu_detCr7Ni)
    % Center system at center of coordinates for NO
    Cu_center_Cu(:,:,m) = Cu_final_coord(:,:,m)-Cu_final_coord(1,:,m);
    Cr7Ni_center_Cu(:,:,m) = Cr7Ni_pump_overlay-Cu_final_coord(1,:,m);
    
    center_Cu_center_Cu(:,m) = Cu_final_coord(1,:,m) - Cu_final_coord(1,:,m);
    center_Cr7Ni_center_Cu(:,m) = [0 0 0] - Cu_final_coord(1,:,m);
    
    % Calcuate rotation matrix
    g_rot = erot(p_pumpCu_detCr7Ni(m,4:6));
    
    % Rotate coordinates and ceters into g frame of ring 1
    Cu_center_Cu_gCu(:,:,m) = g_rot'*Cu_center_Cu(:,:,m)';
    Cr7Ni_center_Cu_gCu(:,:,m) = g_rot'*Cr7Ni_center_Cu(:,:,m)';
    
    center_Cu_center_Cu_gCu(:,m) = g_rot'*center_Cu_center_Cu(:,m);
    center_Cr7Ni_center_Cu_gCu(:,m) = g_rot'*center_Cr7Ni_center_Cu(:,m);
    
    % Translate g matrices to correct position, rotate g matrices into g frame of ring 1 and re center
    gCu_gCu(:,:,m) = (g_rot'*(erot(p_pumpCu_detCr7Ni(m,4:6))+Cu_final_coord(1,:,m)' - Cu_final_coord(1,:,m)'))-center_Cu_center_Cu_gCu(:,m);
    gCr7Ni_gCu(:,:,m) = (g_rot'*(erot([0 0 0])+[0 0 0]' - Cu_final_coord(1,:,m)'))-center_Cr7Ni_center_Cu_gCu(:,m);
    
    % Generate p_vectors
    [phiCr7Ni_Cu(m),thetaCr7Ni_Cu(m)] = vec2ang(center_Cr7Ni_center_Cu_gCu(:,m));
    g_angles_Cr7Ni_Cu(:,m) = eulang(gCr7Ni_gCu(:,:,m));
    p_pumpCr7Ni_detCu(:,m) = [phiCr7Ni_Cu(m),thetaCr7Ni_Cu(m), norm(center_Cr7Ni_center_Cu_gCu(:,m)), g_angles_Cr7Ni_Cu(:,m)'];
    
    counter = counter+1;
    display(counter)
end

%% Test the p_vectors created.

% Set up opt array
Cu_A = [-60.4642 -423.9036 -423.9036];
opt.Det_Ham=struct('S',1/2,'g',[2.0425 2.1604 2.1604],'gFrame',[0 0 0],...
    'HStrain',[55.7976 127.6371 127.6371],'Nucs','65Cu','A',[Cu_A],...
    'lwpp',[1.0594 1.6886]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Det_spindensity=[0.6 0.1 0.1 0.1 0.1]; % spin density for each point/nuclei
opt.Det_coordinates=cu_det'; %coordinates from BDP_Ala_frames in D frame ---> BDP_Dtensor_OriLaserIMDsim % coordinates of the first pump center - best to include these in the most anisotripc frame of the pump hamiltonian

% Pump system PARAMETERS
opt.Pump_number = 1; %define the number of pump centers

opt.Pump_coordinates.N1 = Cr7Ni_pump_overlay';
opt.Pump_Ham.N1=struct('S',1/2,'g',[1.7722 1.7699 1.7342],'gFrame',[0 0 0],...
    'gStrain',[0.0060 0.0032 0.0055],'lwpp',[18.5417]); %if gFrame, and g frame = 0 then corodinates must be inputted in the g frames of the molecules
opt.Pump_spindensity.N1=[1 1 1 1 1 1 1 1]./8;       % spin density for each point/nuclei
opt.Pump_coordinate_centre.N1=[0 0 0]'; % coordinates for the center of the moeity in the starting frame used as center of rotation and translation

opt.Print_coord = 0;

f=figure(5);
hold on
xlabel('X')
ylabel('Y')
zlabel('Z')
m=plot3(opt.Det_coordinates(1,:),opt.Det_coordinates(2,:),opt.Det_coordinates(3,:),'ro');
 m.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
for plotting = 1:5 % do 5 test plots
    plot1 = ceil(rand*length(p_pumpCu_detCr7Ni));
    
    p.N1(1,:) = p_pumpCr7Ni_detCu(:,plot1)';
    
    [opt_out]= OriDEER_multispin_coordinates(p,opt);
    
    % Plot coordinates
    
    m=plot3(opt_out.Pump_coordinates.N1(1,:),opt_out.Pump_coordinates.N1(2,:),opt_out.Pump_coordinates.N1(3,:),'ko');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';
    m=plot3(Cr7Ni_center_Cu_gCu(1,:,plot1),Cr7Ni_center_Cu_gCu(2,:,plot1),Cr7Ni_center_Cu_gCu(3,:,plot1),'kx');
    m.Annotation.LegendInformation.IconDisplayStyle = 'off';

    Pump_coordinates_centre=[sum(opt_out.Pump_coordinates.N1(1,:))/8;...
        sum(opt_out.Pump_coordinates.N1(2,:))/8;...
        sum(opt_out.Pump_coordinates.N1(3,:))/8]';
    % Plot g axes for one center as a check
    pump_g1= erot(opt_out.Pump_Ham.N1.gFrame);
    h = mArrow3(Pump_coordinates_centre,pump_g1(1,:)+Pump_coordinates_centre,'color','red','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Pump_coordinates_centre,pump_g1(2,:)+Pump_coordinates_centre,'color','green','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
    h = mArrow3(Pump_coordinates_centre,pump_g1(3,:)+Pump_coordinates_centre,'color','blue','stemWidth',0.02,...
        'tipWidth',0.05,'facealpha',0.2); %+opt_out.Pump_coordinates.(fn)(:,1) to translate the vector
end

grid on
box on
axis equal
legend('g_{x}','g_{y}','g_{z}')

f=make_plot_nice_p_vector(f);
f.Renderer='opengl';

%% Save the output - three sizes of variables - with information about the cone sizes etc. in input strucutre

%save('p_vectors_only','input','p_pumpCr7Ni_detNO','p_pumpCr7Ni_detCu')

