clear;
clc;

% NB - this file must be placed adjacent to Collected_Output and
% Collected_GLSTAT.
% Ensure Vi, Ri, and index for figure 1 generation are correctly specified.

% Get current folder
output_folder = 'Collected_Output';
glstat_folder = 'Collected_GLSTAT';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% cluster_sizes extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_sizes.csv'));

cluster_sizes = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_sizes{i} = data;
end

%
% cluster_element_count extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_element_count.csv'));

cluster_element_counts = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_element_counts{i} = data;
end

%
% cluster_volumes extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_volumes.csv'));

cluster_volumes = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_volumes{i} = data;
end

%
% cumulative_volume_ratio extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cumulative_volume_ratio.csv'));

cumulative_volume_ratios = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cumulative_volume_ratios{i} = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% cluster_velocities extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_velocities.csv'));

cluster_velocities = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_velocities{i} = data;
end

%
% cluster_angular_velocities extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_angular_velocities.csv'));

cluster_angular_velocities = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_angular_velocities{i} = data;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% cluster_trans_KEs extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_trans_KEs.csv'));

cluster_trans_KEs = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_trans_KEs{i} = data;
end

%
% cluster_rot_KEs extraction
%

% List subfolders
subfolders = dir(fullfile(output_folder, 'MAIN*', 'cluster_rot_KEs.csv'));

cluster_rot_KEs = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    % Extract subfolder name
    [~, subfolderName, ~] = fileparts(subfolders(i).folder);
    folder_name = matlab.lang.makeValidName(subfolderName);

    % Read the CSV file from the current subfolder
    filePath = fullfile(subfolders(i).folder, subfolders(i).name);
    data = readmatrix(filePath);
    cluster_rot_KEs{i} = data;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Global IE, KE extraction
%

% List subfolders
glstat_csv = 'Collected_GLSTAT/glstat_summary.csv';

glstat = readmatrix(glstat_csv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Box plot, kinematic variables distribution v. cluster size (Requiv/Ri)
%

% Initial velocity
Vi = 100;

% Reference total volume
total_volume = sum(cluster_volumes{1}(:,1));
Ri = ((3 * total_volume)/(4 * pi)) ^ 1/3;

% Calculate normalized radii
equivalent_radius = cell(length(subfolders), 1);
normalized_radius = cell(length(subfolders), 1);

for i = 1:length(subfolders)
    equivalent_radius{i} = ((3 * cluster_volumes{i})/(4 * pi)) .^ 1/3;
    normalized_radius{i} = equivalent_radius{i} ./ Ri;
end

% Calculate kinematic variables
for i = 1:length(subfolders)
    CoR_normal{i} = abs(cluster_velocities{i}(:,3) / Vi);
    CoR_tangential{i} = sqrt(cluster_velocities{i}(:,1) .^2 + cluster_velocities{i}(:,2) .^2) / Vi;
    vel_rotational{i} = sqrt(cluster_angular_velocities{i}(:,1) .^ 2 + cluster_angular_velocities{i}(:,2) .^ 2 + cluster_angular_velocities{i}(:,3) .^ 2) * 1000 / (2 * pi) / 1000000; % rad/ms to 10^6 rot/s
    CoR_rotational{i} = sqrt(cluster_angular_velocities{i}(:,1) .^ 2 + cluster_angular_velocities{i}(:,2) .^ 2 + cluster_angular_velocities{i}(:,3) .^ 2) * 1000 .* equivalent_radius{i} / Vi;
end

% Index of interest (see subfolders)
index = 1;

% Bin particles according to size
bin = discretize(normalized_radius{index}, 6);

averageSizes = zeros(1, 6); % Preallocate for 6 bins
for i = 1:6
    averageSizes(i) = mean(normalized_radius{index}(bin == i)); % Calculate mean for each bin
end

figure(1); tiledlayout(1,3);

nexttile;
boxplot(CoR_normal{index}, bin);
ylabel('normal CoR'); xlabel('Requiv/Ri');
ylim([0 1.2]);
xticklabels(arrayfun(@num2str, averageSizes, 'UniformOutput', false));

nexttile;
boxplot(CoR_tangential{index}, bin);
ylabel('tangential CoR'); xlabel('Requiv/Ri');
ylim([0 2]);
xticklabels(arrayfun(@num2str, averageSizes, 'UniformOutput', false));

nexttile;
boxplot(vel_rotational{index}, bin);
ylabel('rotational velocity'); xlabel('Requiv/Ri');
ylim([0 1]);
xticklabels(arrayfun(@num2str, averageSizes, 'UniformOutput', false));


%
% t10, D dist. v. Vi
%

t10 = [];

for i = 1:length(subfolders)
    t10_cluster_volume = 0;
    for j = 1:length(cluster_volumes{i})
        if equivalent_radius{i}(j) < 0.1 * Ri
            t10_cluster_volume = t10_cluster_volume + cluster_volumes{i}(j);
        end
    end
    t10(i) = t10_cluster_volume/total_volume;
end

D = [];

for i = 1:length(subfolders)
    D_cluster_volume = 0;
    for j = 1:length(cluster_volumes{i})
        D_cluster_volume = D_cluster_volume + (cluster_volumes{i}(j) / total_volume) ^ 2;
    end
    D(i) = 1 - D_cluster_volume;
end

figure(2); tiledlayout(1,2);

nexttile;
boxplot(t10);
ylabel('t10'); xlabel('Vi');
ylim([0 0.5]);

nexttile;
boxplot(D);
ylabel('D'); xlabel('Vi');
ylim([0 1]);

%
% IE, KE dist. v. Vi
%

% Calculated normalized final IE & KE
normalized_IE = glstat(:,5) ./ glstat(:,4);
normalized_KE = glstat(:,6) ./ glstat(:,4);

figure(3); tiledlayout(1,2);

nexttile;
boxplot(normalized_IE);
ylabel('IE'); xlabel('Vi');
ylim([0 1]);

nexttile;
boxplot(normalized_KE);
ylabel('KE'); xlabel('Vi');
ylim([0 1]);

%
% Kinematic variables dist. v. Vi, mass-averaged
%

% Calculate mass adjusted kinematic variables & collapse into single array

normalized_cluster_volumes = cell(length(subfolders), 1);

CoR_normal_mass_adjusted_full = [];
CoR_tangential_mass_adjusted_full = [];
CoR_rotational_mass_adjusted_full = [];

for i = 1:length(subfolders)
    normalized_cluster_volumes{i} = cluster_volumes{i} / sum(cluster_volumes{i});

    CoR_normal_mass_adjusted_full(i) = sum(CoR_normal{i} .* normalized_cluster_volumes{i});
    CoR_tangential_mass_adjusted_full(i) = sum(CoR_tangential{i} .* normalized_cluster_volumes{i});
    CoR_rotational_mass_adjusted_full(i) = sum(CoR_rotational{i} .* normalized_cluster_volumes{i});
end


figure(4); tiledlayout(1,3);

nexttile;
boxplot(CoR_normal_mass_adjusted_full);
ylabel('normal CoR'); xlabel('Vi = 80');
ylim([0 1]);

nexttile;
boxplot(CoR_tangential_mass_adjusted_full);
ylabel('tangential CoR'); xlabel('Vi = 80');
ylim([0 1]);

nexttile;
boxplot(CoR_rotational_mass_adjusted_full);
ylabel('rotational velocity'); xlabel('Vi = 80');
ylim([0 0.5]);
