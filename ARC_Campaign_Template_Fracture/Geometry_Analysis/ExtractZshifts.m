clear; 
clc;

%Store node positions of blocky geometry corresponding to theta_1=theta_2=theta3=0
%(no pre-rotation)
nodes = readmatrix("projectile_node_coordinates.csv");
theta_values = readmatrix("theta_values.csv");
DesiredDistanceAboveTarget=0.01; %[mm]

%-------------------------------------------
% View geometry
%-------------------------------------------

% scatter3(nodes(:,1),nodes(:,2),nodes(:,3));
% axis equal

%-------------------------------------------
% Define the "rotation center", the point about which rotation axes extend
% through. Generally defined as the center of mass gotten using LSPP measure tool of type "inertia".
%-------------------------------------------

rotation_center = [0 0 0.05574];

%-------------------------------------------
% Define the rotation angles vectors
%-------------------------------------------

theta_1_vector = theta_values(:,1);
theta_2_vector = theta_values(:,2);
theta_3_vector = theta_values(:,3);

%-------------------------------------------
% Identify the coordinates of the point on the 3D geometry corresponding to the lowest z-coordinate at which impact
% is assumed to first occur (assuming normal impact)
%-------------------------------------------

% Pre-allocate memory for the output matrix
% Rows equal to the total number of combinations of theta_1 and theta_2
min_zValue_matrix = zeros(length(theta_values), 4);

% Counter for rows in output matrix
row_counter = 1;

for i = 1:length(theta_values)
    % Extract angles
    theta_1 = theta_values(i,1);
    theta_2 = theta_values(i,2);
    theta_3 = theta_values(i,3);

    % Rotation matrices
    Rx = [1, 0, 0; 0, cos(theta_1*(pi/180)), -sin(theta_1*(pi/180)); 0, sin(theta_1*(pi/180)), cos(theta_1*(pi/180))];
    Ry = [cos(theta_2*(pi/180)), 0, sin(theta_2*(pi/180)); 0, 1, 0; -sin(theta_2*(pi/180)), 0, cos(theta_2*(pi/180))];
    Rz = [cosd(theta_3), -sind(theta_3), 0; sind(theta_3), cosd(theta_3), 0; 0, 0, 1];

    % Shift nodes to rotation center, rotate, and shift back
    nodes_shifted = nodes - rotation_center;
    nodes_rotated = (Rz * (Ry * (Rx * nodes_shifted.'))).';

    % Add the rotation center back to get the actual positions
    nodes_rotated = nodes_rotated + rotation_center;

    % Find the node with the minimum z-coordinate
    min_z_value = min(nodes_rotated(:, 3));

    % Store the coordinates of the node and the corresponding angles
    min_zValue_matrix(row_counter, :) = [min_z_value, theta_1, theta_2, theta_3];

    % Increment the row counter
    row_counter = row_counter + 1;
end


% Calculate the distance which the particle should be shifted in the
% z-direction to get as close to the target at the start of the simulation
% as possible.

zShift_Value_matrix = min_zValue_matrix;
zShift_Value_matrix(:,1) = DesiredDistanceAboveTarget -  min_zValue_matrix(:,1);

writematrix(zShift_Value_matrix,"Z_shift_matrix.csv");

