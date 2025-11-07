tic


% % % Move current directory to the location where .m file is being executed


%%%=========================================================
%%%=========================================================

%%% Read data from csv files and initialize parameters
node_positions = readmatrix(strcat(pwd,'/Input_files/node_positions.csv'));
node_velocities = readmatrix(strcat(pwd,'/Input_files/node_velocities.csv'));
element_data = readmatrix(strcat(pwd,'/Input_files/element_data.csv'));

%%% Deleted elements contains the list of elements that are eroded in LS-DYNA
%%% e.g., due to their associated timestep falling below minimums. Note that
%%% eroded element/node information is retained in node_positions and element_data
%%% and must be manually removed as is done here.

deleted_elements_data = readmatrix(strcat(pwd,'/Input_files/deleted_elements_data.csv'));



%%% Create folder to store results

folderName = 'Output_files';
% Check if the folder already exists
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

output_path=strcat(pwd,'/',folderName);

%%% Initialize parameters

dx = 0.0001;  % set your dx value here -- good starting guess is 5% or so of typical node separation distance in an element.
bin_size_1 = 0.01;  % partition size for finding immediate neighbors of each element -- good starting guess is the typical node separation distance in an element
bin_size_2=0.03; %partition size for binning "neighbor clusters" of each element -- good starting guess is 5 time the typical node separation distance in an element


%%%=========================================================
%%%=========================================================



%%% Note on the nature of the above spatial discretization levels

%%% CRITICAL:
%%% dx: the maximum separation distance between two nodes of separate elements to consider the elements to be stuck together

%%% CRITICAL:
%%% bin_size_1: the initial discretization of space which is used to find direct neighbors stuck to each element.
%%% the only way for the code to fail in identifying an immediate element neighbor is if the nodes connecting two elements together 
%%% happen to lie on either side of the space partitioning of width bin_size_1. This is already unlikely to happen given that the nodes
%%% are very close together and if the elements are connected by multiple nodes, this would have to happen at all node 
%%% connection points. In principle, however, the odds of this happening would increase with decreasing bin_size_1 because there 
%%% are more partitions created. This therefore presents a TRADE-OFF in selecting the size of bin_size_1:
%%% Increasing bin_size_1 ----> reduce the risk of missing nodes but increase redundancy in node checking (the number of nodes stored in each bin). Also, the matrix holding partitioning data is smaller (fewer, larger bins)
%%% Decreasing bin_size_1 ----> Reduce redunancy in node checking (less nodes stored per smaller-sized bin) but increase the risk of missing neighbor elements. Also, requires larger-sized matrix / cell array to store data (greater number of smaller bins) (this really can cause memory failure).


%%% sub-critical:
%%% bin_size_2: the initial discretization of space which is used to merge together elements into clusters.
%%% A poor choice of bin_size_2 makes the algorithm less efficient, but it won’t fundamentally break cluster identification.
%%% If bin_size_2 is too small, the script still finds the correct fragments, but it takes longer. More bins --> More initial clusters --> More merging required later.
%%% If bin_size_2 is too large, it processes larger bins unnecessarily, slowing down clustering. More elements grouped together initially, leading to larger clusters from the start ,but also more redundant comparisons within bins, increasing computational time in the bin-sorting step.






%%%=========================================================
%%%=========================================================


%%% Remove data pertaining to elements that erode during the course of the
%%% simulation. To find deleted elements, use the Python code I wrote
%%% called ExtractFailedElements.py

%%% 1. Efficiently remove rows from element_data corresponding to Deleted Elements
logical_index_to_remove = ismember(element_data(:, 1), deleted_elements_data);
element_data(logical_index_to_remove, :) = [];

%%% 2. Remove unused nodes from node_positions
logical_nodes_in_use = ismember(node_positions(:, 1), element_data(:, 2:5));
node_positions(~logical_nodes_in_use, :) = [];
node_velocities(~logical_nodes_in_use, :) = [];


%%% Remove rows with any NaN values
element_data(any(isnan(element_data), 2), :) = [];
node_positions(any(isnan(node_positions), 2), :) = [];
node_velocities(any(isnan(node_positions), 2), :) = [];


%%% Create a map for node/element IDs to their indexed positions
nodeID_to_index = containers.Map(node_positions(:, 1), 1:size(node_positions, 1));
elementID_to_index = containers.Map(element_data(:,1), 1:size(element_data, 1));

%%%=========================================================
%%%=========================================================





%%%=========================================================
% % % PHASE 1
%%%=========================================================

%%% Ultimate goal in the following block of code: find neighbors for each element

%%% First, partition 3D space into a grid with mesh spacing of bin_size_1.
%%% bin_size_1 should be made as small as possible but not so small as to
%%% make the associated matrices so large as to cause memory issues.

%%% Find the minimum position in each dimension
min_pos = min(node_positions(:,2:4));

%%% Adjust positions such that the smallest value in each dimension is greater than 0

node_positions_zeroed = node_positions;
node_positions_zeroed(:,2:4) = node_positions(:,2:4) - min_pos + bin_size_1;

%%% Create bins
bins = cell(ceil(max(node_positions_zeroed(:,2:4))/bin_size_1));

%%% Fill bins with node IDs
for i = 1:size(node_positions_zeroed,1)
    bin_indices = ceil(node_positions_zeroed(i,2:4)/bin_size_1);
    bins{bin_indices(1), bin_indices(2), bin_indices(3)} = [bins{bin_indices(1), bin_indices(2), bin_indices(3)}, node_positions_zeroed(i, 1)];
end


% Create a map for node IDs to their element IDs
nodeID_to_element = containers.Map('KeyType','double','ValueType','double');
for i = 1:size(element_data, 1)
    nodeIDs = element_data(i, 2:end);
    for nodeID = nodeIDs
        nodeID_to_element(nodeID) = element_data(i, 1);
    end
end


%%%=========================================================
%%%=========================================================

% Find neighbors for each element
num_elements = size(element_data, 1);
neighbors = cell(num_elements, 1);
parfor i = 1:num_elements
    % disp(i)
    nodes_i = element_data(i, 2:end);
    for nodeID_i = nodes_i
        index_i = nodeID_to_index(nodeID_i);
        pos_i = node_positions_zeroed(index_i, 2:4);
        bin_indices_i = ceil(pos_i/bin_size_1);
        search_indices = [max(1, bin_indices_i-1); min(bin_indices_i+1, size(bins))];
        for x = search_indices(1,1):search_indices(2,1)
            for y = search_indices(1,2):search_indices(2,2)
                for z = search_indices(1,3):search_indices(2,3)
                    potential_neighbors = bins{x,y,z};
                    for nodeID_j = potential_neighbors
                        index_j = nodeID_to_index(nodeID_j);
                        pos_j = node_positions_zeroed(index_j, 2:4);
                        distance = norm(pos_i - pos_j);
                        if distance <= dx
                            % Retrieve elements that nodeID_j belongs to
                            neighbor_elements = nodeID_to_element(nodeID_j);
                            % Add these elements to the neighbor list of current element, eliminating duplicates
                            neighbors{i} = unique([neighbors{i}, neighbor_elements]);
                        end
                    end
                end
            end
        end
    end
end


% Convert cell array to a matrix
num_neighbors = cellfun(@numel, neighbors);
max_neighbors = max(num_neighbors);
neighbor_array = nan(num_elements, max_neighbors + 1);
for i = 1:num_elements
    neighbor_array(i, 1:num_neighbors(i)) = neighbors{i};
end

% Write the result to a csv file
location = strcat(output_path,"/neighbors.csv");
writematrix(neighbor_array, location);

neighbor_matrix = neighbor_array;
neighbor_matrix=[element_data(:,1), neighbor_matrix];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%% RESULT: neighbor_matrix.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Each row corresponds to one element in the projectile. Each column
% of a given row specifies that element IDs in the "deformed" configuration
% of the particle that should be considered neighbors of the "base" element
% corresponding to the row.





%%%=========================================================
% % % PHASE 2
%%%=========================================================



%%% Now, we will bin sort the clusters of neighbors.



%%%=========================================================
%%%=========================================================
%%% Calculate the centroid location of each element and store it in the
%%% matrix called centroid_matrix


% Extract node coordinates from node_positions
node_ids = node_positions(:, 1);
node_coords = node_positions(:, 2:4);

% Extract element node IDs from element_data
element_ids = element_data(:, 1);
element_node_ids = element_data(:, 2:5);

% Map the element node IDs to its row position in node_positions
[~, loc] = ismember(element_node_ids, node_ids);

% Create an Nx8x3 matrix where N is the number of elements
% 8 is the number of nodes in each element
% 3 represents the x, y, and z coordinates
element_node_coords = reshape(node_coords(loc, :), [size(element_node_ids), 3]);
% Calculate the centroids for each element by averaging the coordinates of its nodes
centroids = squeeze(mean(element_node_coords, 2));

% Combining the centroids with their corresponding element IDs
centroids_matrix = [element_ids, centroids];



%%%==============================================
%%%==============================================
%Bin sort each element based on the position of its centroid.


% Extract centroid coordinates
centroids = centroids_matrix(:, 2:4);

% Compute min and max for each dimension
min_coords = min(centroids);
max_coords = max(centroids);

% Determine number of bins for each dimension
num_bins = ceil((max_coords - min_coords) / bin_size_2);

% Initialize the 3D bins to bin each element that makes up the sphere
bin_elements = cell(num_bins(1), num_bins(2), num_bins(3));

% Assign each element to its bin
for i = 1:size(centroids, 1)
    bin_indices = max(1, ceil((centroids(i,:) - min_coords) ./ bin_size_2));
    if isempty(bin_elements{bin_indices(1), bin_indices(2), bin_indices(3)})
        bin_elements{bin_indices(1), bin_indices(2), bin_indices(3)} = centroids_matrix(i, 1);
    else
        bin_elements{bin_indices(1), bin_indices(2), bin_indices(3)} = [bin_elements{bin_indices(1), bin_indices(2), bin_indices(3)}; centroids_matrix(i, 1)];
    end
end



%%% Within each bin, combine clusters of neighboring elements that have at
%%% least one common element ID.



% Iterate through each bin in bin_elements
[bin_dims1, bin_dims2, bin_dims3] = size(bin_elements);

numBins = bin_dims1 * bin_dims2 * bin_dims3;

% Precompute bin data and indices
binData = cell(numBins, 1);
binIndices = zeros(numBins, 3);

for linearIdx = 1:numBins
    [i, j, k] = ind2sub([bin_dims1, bin_dims2, bin_dims3], linearIdx);
    binData{linearIdx} = bin_elements{i, j, k};
    binIndices(linearIdx, :) = [i, j, k];
end

binIndices_x = binIndices(:,1);
binIndices_y = binIndices(:,2);
binIndices_z = binIndices(:,3);



parfor linearIdx = 1:numBins

    element_ids = binData{linearIdx};
    % Convert the linear index back to 3D indices
    i = binIndices_x(linearIdx);
    j = binIndices_y(linearIdx);
    k = binIndices_z(linearIdx);



    if isempty(element_ids)
        combined_clustersInBins_temp{linearIdx} = [];
        continue;
    end
    % disp(strcat("Cell Index: ", string(i),'...', string(j),'...', string(k)))
    % Construct NeighborList_inBin as a cell array for the current bin
    
    NeighborList_inBin = cell(length(element_ids), 1);

    %NieghborList_inBin is initially a cell array where each "row"
    %corresponds to one of the elements in the current bin
    %(discretized using bin_dims3). Each column of this row gives
    %the element IDs that are neighbors to the base element. The
    %goal of the following code is to group together elements and
    %their neighbors that have at least one neighbor in common. The
    %goal of course being to work toward identifying the collection of elements
    %that make up unique clusters. 
   

    for m = 1:length(element_ids)
        idx = elementID_to_index(element_ids(m));
        NeighborList_inBin{m} = neighbor_matrix(idx,~isnan(neighbor_matrix(idx,:)));
        
    end

    % Initialize reduced cell array
    clusters_in_bin = {};
    
    while ~isempty(NeighborList_inBin)
        % Start with the first cell
        current_cell = NeighborList_inBin{1};
        NeighborList_inBin(1) = [];
        
        % Check for overlap with remaining cells
        overlap = cellfun(@(x) any(ismember(x, current_cell)), NeighborList_inBin);
        
        % While there's overlap, merge and continue checking
        while any(overlap)
            % Add overlapping cells to current cell
            current_cell = unique([current_cell, NeighborList_inBin{overlap}]);
            
            % Remove added cells from the list
            NeighborList_inBin = NeighborList_inBin(~overlap);
            
            % Check for overlap again
            overlap = cellfun(@(x) any(ismember(x, current_cell)), NeighborList_inBin);
        end
        
        % Add the merged cell to the output
        clusters_in_bin{end+1} = current_cell;

        % Each row of reduced cells now contains a cluster of
        % elements that are all related by the fact that they share
        % at least one element as neighbors. At this point, for computational efficiency, only
        % common neighbors between elements in the same spatial bin
        % (where binning was done using bin_dims3) have een checked
        % and are included for a given "clusters_in_bin" cell array.
    end
    
    % Convert "clusters_in_bin" cell array associated with the current "bin" to a padded matrix
    max_len = max(cellfun('length', clusters_in_bin));
    padded_matrix = nan(length(clusters_in_bin), max_len);
    for m = 1:length(clusters_in_bin)
        current_len = length(clusters_in_bin{m});
        padded_matrix(m, 1:current_len) = clusters_in_bin{m};
    end
    
    % Create the cell array "combined_clustersInBins" to store
    % clusters_in_bin for each bin
    combined_clustersInBins_temp{linearIdx} = padded_matrix;
end


% Reassemble the results into the 3D array
combined_clustersInBins = cell(bin_dims1, bin_dims2, bin_dims3);
for linearIdx = 1:numBins
    i = binIndices_x(linearIdx);
    j = binIndices_y(linearIdx);
    k = binIndices_z(linearIdx);
    combined_clustersInBins{i, j, k} = combined_clustersInBins_temp{linearIdx};
end




%%% Turn combined_clustersInBins into a single matrix called combined_ClustersInBins_Matrix

% Step 1: Find the maximum number of columns
max_cols = 0;
[bin_dims1, bin_dims2, bin_dims3] = size(combined_clustersInBins);
for i = 1:bin_dims1
    for j = 1:bin_dims2
        for k = 1:bin_dims3
            current_matrix = combined_clustersInBins{i, j, k};
            if ~isempty(current_matrix)
                max_cols = max(max_cols, size(current_matrix, 2));
            end
        end
    end
end

% Step 2 and 3: Pad and concatenate
combined_ClustersInBins_Matrix = [];
for i = 1:bin_dims1
    for j = 1:bin_dims2
        for k = 1:bin_dims3
            current_matrix = combined_clustersInBins{i, j, k};
            if ~isempty(current_matrix)
                % Padding
                pad_size = max_cols - size(current_matrix, 2);
                padded_matrix = [current_matrix, nan(size(current_matrix, 1), pad_size)];
                
                % Concatenate
                combined_ClustersInBins_Matrix = [combined_ClustersInBins_Matrix; padded_matrix];
            end
        end
    end
end



%%%=========================================================
% % % PHASE 3
%%%=========================================================



%%% Finally, merge all rows of combined_ClustersInBins_Matrix that have at least one
%%% common element.

% Convert matrix back to a cell array
combined_ClustersInBin_Array = arrayfun(@(i) combined_ClustersInBins_Matrix(i,~isnan(combined_ClustersInBins_Matrix(i,:))), 1:size(combined_ClustersInBins_Matrix, 1), 'UniformOutput', false);

% Initialize a cell array to hold the identified, final particle clusters.
global_clusters = {};

while ~isempty(combined_ClustersInBin_Array)
    % Start with the first cell
    current_cell = unique(combined_ClustersInBin_Array{1});
    combined_ClustersInBin_Array(1) = [];
    % disp(strcat('combined_ClustersInBin_Array','...', string(size(combined_ClustersInBin_Array,2))))
    % Check for overlap with remaining cells
    overlap = cellfun(@(x) any(ismember(x, current_cell)), combined_ClustersInBin_Array);
    
    % While there's overlap, merge and continue checking
    while any(overlap)
        % Add overlapping cells to current cell
        current_cell = unique([current_cell, combined_ClustersInBin_Array{overlap}]);
        
        % Remove added cells from the list
        combined_ClustersInBin_Array = combined_ClustersInBin_Array(~overlap);
        
        % Check for overlap again
        overlap = cellfun(@(x) any(ismember(x, current_cell)), combined_ClustersInBin_Array);
    end
    
    % Add the merged cell to the output
    global_clusters{end+1} = current_cell;
end




% Number of columns in the cell array
nCols = length(global_clusters);

% Find the maximum length of row vectors in the cell array
maxLen = max(cellfun(@length, global_clusters));

% Initialize an empty sparse matrix
S = sparse(nCols, maxLen);

% Populate the sparse matrix
for i = 1:nCols
    data = global_clusters{i};
    S(i, 1:length(data)) = data;
end



location = strcat(output_path,"/saved_sparse_cluster.mat");
save(location,"S");


%%%=========================================================
%%%=========================================================
%%% Get information on cluster sizes
%%%=========================================================
%%%=========================================================


%%%-------------------------------------------------------
% Get the element number density of each cluster.
%%%------------------------------------------------------
clusterElCount = zeros(1, size(S, 1));

% Calculate the length of each row
for i = 1:size(S, 1)
    clusterElCount(i) = nnz(S(i,:));
end

location = strcat(output_path,"/cluster_element_count.csv");
writematrix(clusterElCount, location);




%%%-------------------------------------------------------
% Get the volume of each cluster.
%%%------------------------------------------------------

% Compute volumes for each element


% Initialize array for volumes
volumes = zeros(size(element_data, 1), 1);

% Compute volumes of the hexahedrons
for i = 1:size(element_data, 1)
    nodes = element_data(i, 2:end);
    positions = zeros(4, 3);
    for j = 1:4
        positions(j, :) = node_positions(nodeID_to_index(nodes(j)), 2:4);
    end
    
    volumes(i) = 1/6 * abs(dot((positions(2,:)-positions(1,:)), cross((positions(3,:)-positions(1,:)), (positions(4,:)-positions(1,:)))));
    
end

% Create a map for element IDs to their volumes
elID_to_volume = containers.Map(element_data(:, 1), volumes);



% Get the volume of each cluster

% Initialize vector to store cluster volumes
num_clusters = size(S, 1);
cluster_volumes = zeros(num_clusters, 1);

% Calculate cluster volumes
for i = 1:num_clusters
    % Get element IDs for this cluster
    nonZeroElements = find(S(i,:));
    element_IDs = full(S(i,nonZeroElements));
    % disp(strcat("cluster volume calculation", string(i)));
    % Sum volumes of elements in this cluster
    for j = 1:length(element_IDs)
        if elID_to_volume.isKey(element_IDs(j))
            cluster_volumes(i) = cluster_volumes(i) + elID_to_volume(element_IDs(j));
        else
            disp("error on mapping to volume ID")
        end
    end
end


% Write the results to a csv file

cluster_sizes = (3/(4*pi)*cluster_volumes).^(1/3).*1000;  %length units: um

total_volume=sum(cluster_volumes);
sorted_cluster_volumes = sort(cluster_volumes);
cumulative_volume_ratio = cumsum(sorted_cluster_volumes)./total_volume;


location = strcat(output_path,"/cluster_volumes.csv");
writematrix(cluster_volumes, location);

location = strcat(output_path,"/cluster_sizes.csv");
writematrix(cluster_sizes, location);

location = strcat(output_path,"/cumulative_volume_ratio.csv");
writematrix(cumulative_volume_ratio, location);


%%%=========================================================
%%%=========================================================




%%%=========================================================
%%%=========================================================
%%% Extract cluster translational velocities
%%%=========================================================
%%%=========================================================



%%% Approximate each cluster's center of mass velocity as the average
%%% of the node velocities that make up the cluster.

% Extract just the velocity values into a separate matrix.
velocity_values = node_velocities(:, 2:4);

% Preallocate a matrix to store cluster velocities.
cluster_velocities = zeros(size(S, 1), 3);

% Loop over each cluster in S.
for i = 1:size(S, 1)
    nonZeroElements = find(S(i,:));
    cluster = full(S(i,nonZeroElements));

    % Convert element IDs to indices in the element_data matrix.
    indices = cell2mat(values(elementID_to_index, num2cell(cluster)));

    % Extract all nodes associated with this cluster.
    nodes_in_cluster = unique(element_data(indices, 2:end));

    % Convert node IDs to indices in the velocity_values matrix.
    velocity_indices = cell2mat(values(nodeID_to_index, num2cell(nodes_in_cluster)));

    % Sum velocities of all nodes in the cluster.
    sum_velocity = sum(velocity_values(velocity_indices, :), 1);

    % Compute the average velocity for the cluster.
    cluster_velocities(i, :) = sum_velocity / length(nodes_in_cluster);
end

location = strcat(output_path,"/cluster_velocities.csv");
writematrix(cluster_velocities, location);




%%%=========================================================
%%%=========================================================
%%% Extract cluster angular velocities
%%%=========================================================
%%%=========================================================


numClusters = size(S, 1);
cluster_AngularVel = zeros(numClusters, 3); % Preallocating for speed

for i = 1:numClusters  
    % Extracting elements of the current cluster
    nonZeroElements = find(S(i,:));
    elementIDs = full(S(i,nonZeroElements));
    
    % Extracting the node IDs of these elements
    allNodeIDs = unique(element_data(ismember(element_data(:, 1), elementIDs), 2:end));
    
    % Sampling a subset of nodes if there are too many nodes
    maxNodes = 100; % for example
    if length(allNodeIDs) > maxNodes
        sampledNodeIDs = datasample(setdiff(allNodeIDs, allNodeIDs(1)), maxNodes, 'Replace', false);
        allNodeIDs = [allNodeIDs(1); sampledNodeIDs];
    end
    
    referenceNodePos = node_positions(nodeID_to_index(allNodeIDs(1)), 2:end);
    referenceNodeVel = node_velocities(nodeID_to_index(allNodeIDs(1)), 2:end);

    A = [];  % system matrix
    b = [];  % RHS vector
    
    for j = 2:length(allNodeIDs)
        nodeID = allNodeIDs(j);
        nodePos = (node_positions(nodeID_to_index(nodeID), 2:end) - referenceNodePos)';
        relativeVel = (node_velocities(nodeID_to_index(nodeID), 2:end) - referenceNodeVel)';
        
        skewMat = [0, nodePos(3), -nodePos(2); 
                   -nodePos(3), 0, nodePos(1);
                   nodePos(2), -nodePos(1), 0];
        
        A = [A; skewMat];
        b = [b; relativeVel];
    end
    
    angularVel = A \ b;  % Solving the system
    cluster_AngularVel(i, :) = angularVel';
    % disp(strcat(string(i), " of ", string(numClusters)))
end

writematrix(cluster_AngularVel, strcat(pwd,'/Output_files/cluster_angular_velocities.csv'));



fprintf ( 1, 'Run Time:  %.2f\n', toc );



%%%=========================================================
%%%=========================================================
%%% Extract cluster kinetic energies
%%%=========================================================
%%%=========================================================




%%%=========================================================
%%% USER INPUT
%%%=========================================================

ParticleMassDensity = 2.44000E-6;
NodesPerElement = 4; % % % Code must be modified if this is not true.





%%%=========================================================
%%% Read data from csv files and initialize parameters
%%%=========================================================

element_data = readmatrix(strcat(pwd,'/Input_files/element_data.csv'));
load(strcat(pwd,'/Output_files/saved_sparse_cluster.mat'));
cluster_matrix = full(S);
node_positions = readmatrix(strcat(pwd,'/Input_files/node_positions.csv'));
node_velocities = readmatrix(strcat(pwd,'/Input_files/node_velocities.csv'));
cluster_velocities = readmatrix(strcat(pwd,'/Output_files/cluster_velocities.csv'));
cluster_angular_velocities = readmatrix(strcat(pwd,'/Output_files/cluster_angular_velocities.csv'));
cluster_volumes = readmatrix(strcat(pwd,'/Output_files/cluster_volumes.csv'));




%%% Create a map for node/element IDs to their indexed positions
nodeID_to_index = containers.Map(node_positions(:, 1), 1:size(node_positions, 1));
elementID_to_index = containers.Map(element_data(:,1), 1:size(element_data, 1));


% Number of clusters
num_clusters = size(cluster_matrix, 1);




% % % ---------------------------------------------------
% % % Method for calculating rotational kinetic energy w/o explicitly calculating the cluster axis of rotation.
% % % ---------------------------------------------------


cluster_TotalKE = zeros(num_clusters,1);
cluster_RotKE = zeros(num_clusters,1);
cluster_transKE = zeros(num_clusters,1);


for i = 1:num_clusters
    
    cluster_elements = nonzeros(cluster_matrix(i, :));
    cluster_elements = cluster_elements(~isnan(cluster_elements));  % Removing NaNs

    % Iterate through each element of cluster i
    for j = 1:length(cluster_elements)

        element_id = cluster_elements(j);
        % Extract node IDs for the element
        node_ids = element_data(elementID_to_index(element_id), 2:5);
        
        % Calculate element volume
        positions = zeros(4, 3);
        for k = 1:NodesPerElement
            positions(k, :) = node_positions(nodeID_to_index(node_ids(k)), 2:4);
        end
    
        element_volume = 1/6 * abs(dot((positions(2,:)-positions(1,:)), cross((positions(3,:)-positions(1,:)), (positions(4,:)-positions(1,:)))));
        element_mass = element_volume * ParticleMassDensity;
        node_mass = element_mass / NodesPerElement;


        % Iterate through each node of element j
        for k = 1:NodesPerElement
            node_id = node_ids(k);
            
            v = norm(node_velocities(nodeID_to_index(node_id), 2:end));  % Subtracting translational velocity
            cluster_TotalKE(i) = cluster_TotalKE(i) + (1/2)*node_mass*v^2;

        end

    end

    cluster_mass = cluster_volumes(i)*ParticleMassDensity;
    cluster_speed = norm(cluster_velocities(i,:));
    cluster_transKE(i) = 1/2*cluster_mass*cluster_speed^2;

    cluster_RotKE(i) = cluster_TotalKE(i) - cluster_transKE(i); %cluster_TotalKE was calculated in the method 1 section. 

end


location = strcat(output_path,"/cluster_trans_KEs.csv");
writematrix(cluster_transKE, location);

location = strcat(output_path,"/cluster_rot_KEs.csv");
writematrix(cluster_RotKE, location);


% % % =========================================================
% % % Evaluate whether a fragment of significant size has a chance to strike the target. If low chance, we set a flag to terminate the simulation early.
% % % =========================================================


% SET PARAMETERS
Min_Cluster_Els = 100;
MinDownwardV = -10; % [m/s]
MinRotV = 15; % [m/s]


% initialize TERMINATE status variable
TERMINATE = 'y';


TotalElements = length(element_data);

% Extract the node positions into a separate matrix.
position_values = node_positions(:, 2:4);

% Loop over each cluster in S.
for i = 1:size(S, 1)
    
    nonZeroElements = find(S(i,:));
    cluster = full(S(i,nonZeroElements));

    % Convert element IDs to indices in the element_data matrix.
    indices = cell2mat(values(elementID_to_index, num2cell(cluster)));

    % Extract all nodes associated with this cluster.
    nodes_in_cluster = unique(element_data(indices, 2:end));

    % Convert node IDs to indices in the position_values matrix.
    position_indices = cell2mat(values(nodeID_to_index, num2cell(nodes_in_cluster)));

    % Matrix storing the positions of each point in the cluster
    points_in_cluster = position_values(position_indices, :);

    clusterElementCount = nnz(S(i,:));

    if clusterElementCount > Min_Cluster_Els

        % ------------------------------------------------------
        % If a sufficiently large particle is moving toward the target, the simulation should not yet be terminated.
        % ------------------------------------------------------
        
        % Compute the center of mass
        center_of_mass = mean(points_in_cluster, 1);
        
        if cluster_velocities(i,3) < MinDownwardV && center_of_mass(3) > 0
            TERMINATE = 'n';
            break;  % Exit the loop

        else
            % ------------------------------------------------------
            % If the maximum distance of a point from the particle CoM is smaller than the CoM z coordinate, further contact could happen as the fragment rotates.
            % so the simulation should not yet be terminated
            % ------------------------------------------------------

            
            % Compute Euclidean distances of all points from the center of mass 
            distances = vecnorm(points_in_cluster - center_of_mass, 2, 2);
            % Find the maximum distance
            cluster_reach = max(distances);

            if cluster_reach > center_of_mass(3) && cluster_reach*norm(cluster_AngularVel(i, :)) > MinRotV && center_of_mass(3) > 0
                TERMINATE = 'n';
                break;
            end  
        end
    end
end




% % % Get the path to the MAIN_*_*_* scratch directory holding the simulation results.

% Get current working directory
currentDir = pwd;

% Traverse up the directory tree until we find "Analyze"
while true
    % List all items in the current directory
    contents = dir(currentDir);
    
    % Check if "Analyze" is a folder in the current directory
    isAnalyzePresent = any(strcmp({contents.name}, 'Analyze') & [contents.isdir]);
    
    if isAnalyzePresent
        parentOfAnalyze = currentDir; % The parent folder of "Analyze"
        break; % Stop searching once found
    end
    
    % Move up one directory level
    parentDir = fileparts(currentDir);
    
    % If we reach the root and haven't found "Analyze", exit
    if strcmp(parentDir, currentDir)
        error('No "Analyze" folder found in any parent directory.');
    end
    
    currentDir = parentDir;
end


% Define the output file path in the parent directory of "Analyze"
outputFile = fullfile(parentOfAnalyze, 'TERMINATION_STATUS.txt');

% Open the file for writing
fileID = fopen(outputFile, 'w');

% Write the value of TERMINATE to the file
fprintf(fileID, '%s\n', TERMINATE);

% Close the file
fclose(fileID);







