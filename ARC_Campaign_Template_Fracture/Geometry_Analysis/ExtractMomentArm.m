clear; 
clc;

% Store node positions of blocky geometry corresponding to theta_1=theta_2=0
% (no pre-rotation)
nodes=readmatrix("projectile_node_coordinates.csv");

%Theta1_Outliers = readmatrix("G:\My Drive\ActivityModules\NonFractureRebound\ACTIVITIES\MomentArm_Analysis\HalfSpheroid\Theta1_Outliers_HalfSphere.csv");
%Theta2_Outliers = readmatrix("G:\My Drive\ActivityModules\NonFractureRebound\ACTIVITIES\MomentArm_Analysis\HalfSpheroid\Theta2_Outliers_HalfSphere.csv");

%-------------------------------------------
% View geometry
%-------------------------------------------

% scatter3(nodes(:,1),nodes(:,2),nodes(:,3));
% axis equal

%-------------------------------------------
% Define the center of mass position
% Gotten using LSPP measure tool of type "inertia"
%-------------------------------------------

CoM_position = [0.6,100,60.46];

%-------------------------------------------
% Define the rotation angles vectors
%-------------------------------------------

theta_1_vector = linspace(0,336,15);
theta_2_vector = linspace(0,336,15);

%-------------------------------------------
% Identify the coordinates of the point on the 3D geometry corresponding to the lowest z-coordinate at which impact
% is assumed to first occur (assuming normal impact)
%-------------------------------------------

% Pre-allocate memory for the output matrix
% Rows equal to the total number of combinations of theta_1 and theta_2
output_matrix = zeros(length(theta_1_vector) * length(theta_2_vector), 5);

% Counter for rows in output matrix
row_counter = 1;

for i = 1:length(theta_1_vector)
    for j = 1:length(theta_2_vector)
        % Extract angles
        theta_1 = theta_1_vector(i);
        theta_2 = theta_2_vector(j);

        % Rotation matrices
        Rx = [1, 0, 0; 0, cos(theta_1*(pi/180)), -sin(theta_1*(pi/180)); 0, sin(theta_1*(pi/180)), cos(theta_1*(pi/180))];
        Ry = [cos(theta_2*(pi/180)), 0, sin(theta_2*(pi/180)); 0, 1, 0; -sin(theta_2*(pi/180)), 0, cos(theta_2*(pi/180))];

        % Shift nodes to rotation center, rotate, and shift back
        nodes_shifted = nodes - CoM_position;
        nodes_rotated = (Ry * (Rx * nodes_shifted.')).';

        % Add the rotation center back to get the actual positions
        nodes_rotated = nodes_rotated + CoM_position;

        % Find the node with the minimum z-coordinate
        tolerance = 0.01;

        [~, min_z_index] = min(nodes_rotated(:, 3));
        % Initialization before the loop
        minZ_nodes = nodes_rotated(min_z_index, :); % Start with the minimum z node
        sumZ_nodes = minZ_nodes; % Initialize sum with the coordinates of the node with the minimum z
        LoopCount = 1; % Counting the initial node
        
        % Loop through all nodes to find those within the tolerance
        for k = 1:size(nodes_rotated, 1)
            % Check if the node is within tolerance and not the initially found minimum
            if k ~= min_z_index && abs(nodes_rotated(k, 3) - nodes_rotated(min_z_index, 3)) < tolerance
                sumZ_nodes = sumZ_nodes + nodes_rotated(k, :); % Accumulate the coordinates
                LoopCount = LoopCount + 1; % Increment the count of nodes within tolerance
            end
        end
        
        % Calculate the average coordinates of the nodes within tolerance
        avg_minZ_node_coord = sumZ_nodes / LoopCount;
        
        % Continue with the rest of the code...


        % Store the coordinates of the node and the corresponding angles
        % --- output_matrix(row_counter, :) = [nodes_rotated(min_z_index, :), theta_1, theta_2];
        output_matrix(row_counter, :) = [avg_minZ_node_coord, theta_1, theta_2];


        % Increment the row counter
        row_counter = row_counter + 1;
    end
end


%-------------------------------------------
% Store the coordinate of the head of the vector that points from the
% estimated initial impact point of the particle geometry to the particle's
% center of mass.
%-------------------------------------------

output_matrix(:, 1:3) = CoM_position - output_matrix(:, 1:3);

%-------------------------------------------
% Calculate the moment arm of the CoM relative to the estimated impact point for each orientation and store 
% result in MomentArm_matrix
% MomentArm_matrix(:,1) - calculated moment arm
% MomentArm_matrix(:,2) - theta_1
% MomentArm_matrix(:,3) - theta_2
%-------------------------------------------

MomentArm_matrix = zeros(length(theta_1_vector) * length(theta_2_vector), 3);
MomentArm = sqrt(output_matrix(:,1).^2 + output_matrix(:,2).^2);
MomentArm_matrix(:,1) = MomentArm;
MomentArm_matrix(:,2) = output_matrix(:,4);
MomentArm_matrix(:,3) = output_matrix(:,5);

writematrix(MomentArm_matrix,"projectile_moment_arm_matrix.csv");


%-------------------------------------------
% Create surface of calculated moment arm spanning over theta_1, theta_2
% space.
%-------------------------------------------

zVec = MomentArm_matrix(:,1)./100;
xVec = MomentArm_matrix(:,2);
yVec = MomentArm_matrix(:,3);

[xGrid, yGrid] = meshgrid(linspace(min(xVec), max(xVec), 15), ...
                          linspace(min(yVec), max(yVec), 15));

% Interpolate Z values on the grid
zGrid = griddata(xVec, yVec, zVec, xGrid, yGrid, 'natural');

figure(2)
% Plotting the surface
%%% surf(xGrid, yGrid, zGrid);

contourf(xGrid, yGrid, zGrid,5);
colormap(jet);
colorbar;
% clim([0   0.34]);


xlabel('$\theta_1$','interpreter','latex','FontSize',15);
ylabel('$\theta_2$','interpreter','latex','FontSize',15);
zlabel('Moment arm $\mu\mathrm{m}$','interpreter','latex','FontSize',15);
title('Normalized Moment Arm Length: Half Spheroid Geometry','interpreter','latex','FontSize',15);

