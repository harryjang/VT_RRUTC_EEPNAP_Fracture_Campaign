clear all; clc;


% % ==============================================
% % Set constants for polysuperquadric geometry B
% % ==============================================

actual_particle_volume=0.000826242;
R = (3/(4*pi)*actual_particle_volume)^(1/3)*1000; %Initial impactor radius in microns
rho = 2.44000E-6; %particle mass density [kg/mm^3]

% ==========================================================
% Calculate t_10 values for polysuperquadric geometry B
% ==========================================================

[t10_vector_Vni_60, Orientations_Vni_60_t10, NumBodyFragments_Vni_60, LargestFragment_Vni_60] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_60\N_100\collected_output");
[t10_vector_Vni_80, Orientations_Vni_80_t10, NumBodyFragments_Vni_80, LargestFragment_Vni_80] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_80\N_100\collected_output");
[t10_vector_Vni_100, Orientations_Vni_100_t10, NumBodyFragments_Vni_100, LargestFragment_Vni_100] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_100\N_100\collected_output");
[t10_vector_Vni_120, Orientations_Vni_120_t10, NumBodyFragments_Vni_120, LargestFragment_Vni_120] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_120\N_100\collected_output");
[t10_vector_Vni_140, Orientations_Vni_140_t10, NumBodyFragments_Vni_140, LargestFragment_Vni_140] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_140\N_100\collected_output");
[t10_vector_Vni_160, Orientations_Vni_160_t10, NumBodyFragments_Vni_160, LargestFragment_Vni_160] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_160\N_100\collected_output");
[t10_vector_Vni_180, Orientations_Vni_180_t10, NumBodyFragments_Vni_180, LargestFragment_Vni_180] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_180\N_100\collected_output");
[t10_vector_Vni_200, Orientations_Vni_200_t10, NumBodyFragments_Vni_200, LargestFragment_Vni_200] = Extract_t10_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_200\N_100\collected_output");

NumBodyFragments_all = [NumBodyFragments_Vni_60; NumBodyFragments_Vni_80; NumBodyFragments_Vni_100; NumBodyFragments_Vni_120; NumBodyFragments_Vni_140; NumBodyFragments_Vni_160; NumBodyFragments_Vni_180; NumBodyFragments_Vni_200];
LargestFragments_all = [LargestFragment_Vni_60; LargestFragment_Vni_80; LargestFragment_Vni_100; LargestFragment_Vni_120; LargestFragment_Vni_140; LargestFragment_Vni_160; LargestFragment_Vni_180; LargestFragment_Vni_200];
t10_all = [t10_vector_Vni_60; t10_vector_Vni_80; t10_vector_Vni_100; t10_vector_Vni_120; t10_vector_Vni_140; t10_vector_Vni_160; t10_vector_Vni_180; t10_vector_Vni_200];

% % % --- Calculate quantiles of t_10 distribution over full range of impact orientations.
Vi_60mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_60);
Vi_80mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_80);
Vi_100mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_100);
Vi_120mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_120);
Vi_140mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_140);
Vi_160mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_160);
Vi_180mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_180);
Vi_200mps_quartiles_t10_PolySup = GetQuartiles(t10_vector_Vni_200);
Q0_t10_PolySup =  [Vi_60mps_quartiles_t10_PolySup(1); Vi_80mps_quartiles_t10_PolySup(1); Vi_100mps_quartiles_t10_PolySup(1); Vi_120mps_quartiles_t10_PolySup(1); Vi_140mps_quartiles_t10_PolySup(1); Vi_160mps_quartiles_t10_PolySup(1); Vi_180mps_quartiles_t10_PolySup(1); Vi_200mps_quartiles_t10_PolySup(1)];
Q1_t10_PolySup =  [Vi_60mps_quartiles_t10_PolySup(2); Vi_80mps_quartiles_t10_PolySup(2); Vi_100mps_quartiles_t10_PolySup(2); Vi_120mps_quartiles_t10_PolySup(2); Vi_140mps_quartiles_t10_PolySup(2); Vi_160mps_quartiles_t10_PolySup(2); Vi_180mps_quartiles_t10_PolySup(2); Vi_200mps_quartiles_t10_PolySup(2)];
Q2_t10_PolySup =  [Vi_60mps_quartiles_t10_PolySup(3); Vi_80mps_quartiles_t10_PolySup(3); Vi_100mps_quartiles_t10_PolySup(3); Vi_120mps_quartiles_t10_PolySup(3); Vi_140mps_quartiles_t10_PolySup(3); Vi_160mps_quartiles_t10_PolySup(3); Vi_180mps_quartiles_t10_PolySup(3); Vi_200mps_quartiles_t10_PolySup(3)];
Q3_t10_PolySup =  [Vi_60mps_quartiles_t10_PolySup(4); Vi_80mps_quartiles_t10_PolySup(4); Vi_100mps_quartiles_t10_PolySup(4); Vi_120mps_quartiles_t10_PolySup(4); Vi_140mps_quartiles_t10_PolySup(4); Vi_160mps_quartiles_t10_PolySup(4); Vi_180mps_quartiles_t10_PolySup(4); Vi_200mps_quartiles_t10_PolySup(4)];
Q4_t10_PolySup =  [Vi_60mps_quartiles_t10_PolySup(5); Vi_80mps_quartiles_t10_PolySup(5); Vi_100mps_quartiles_t10_PolySup(5); Vi_120mps_quartiles_t10_PolySup(5); Vi_140mps_quartiles_t10_PolySup(5); Vi_160mps_quartiles_t10_PolySup(5); Vi_180mps_quartiles_t10_PolySup(5); Vi_200mps_quartiles_t10_PolySup(5)];













% % ==========================================================
% % Calculate Gini indices for polysuperquadric geometry B
% % ==========================================================

[Gini_index_vector_Vni_60, Orientations_Vni_60_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_60\N_100\collected_output");
[Gini_index_vector_Vni_80, Orientations_Vni_80_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_80\N_100\collected_output");
[Gini_index_vector_Vni_100, Orientations_Vni_100_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_100\N_100\collected_output");
[Gini_index_vector_Vni_120, Orientations_Vni_120_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_120\N_100\collected_output");
[Gini_index_vector_Vni_140, Orientations_Vni_140_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_140\N_100\collected_output");
[Gini_index_vector_Vni_160, Orientations_Vni_160_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_160\N_100\collected_output");
[Gini_index_vector_Vni_180, Orientations_Vni_180_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_180\N_100\collected_output");
[Gini_index_vector_Vni_200, Orientations_Vni_200_Gini] = Extract_GiniIndex_set("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Poly_SuperQuadric\GEOMETRY_B\CompareReboundStatsToSphere\G_10_Vi_200\N_100\collected_output");

% % % --- Calculate quantiles of the mass dispersion index over full range of impact orientations.
Vi_60mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_60);
Vi_80mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_80);
Vi_100mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_100);
Vi_120mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_120);
Vi_140mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_140);
Vi_160mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_160);
Vi_180mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_180);
Vi_200mps_quartiles_Gini_PolySup = GetQuartiles(Gini_index_vector_Vni_200);
Q0_Gini_PolySup =  [Vi_60mps_quartiles_Gini_PolySup(1); Vi_80mps_quartiles_Gini_PolySup(1); Vi_100mps_quartiles_Gini_PolySup(1); Vi_120mps_quartiles_Gini_PolySup(1); Vi_140mps_quartiles_Gini_PolySup(1); Vi_160mps_quartiles_Gini_PolySup(1); Vi_180mps_quartiles_Gini_PolySup(1); Vi_200mps_quartiles_Gini_PolySup(1)];
Q1_Gini_PolySup =  [Vi_60mps_quartiles_Gini_PolySup(2); Vi_80mps_quartiles_Gini_PolySup(2); Vi_100mps_quartiles_Gini_PolySup(2); Vi_120mps_quartiles_Gini_PolySup(2); Vi_140mps_quartiles_Gini_PolySup(2); Vi_160mps_quartiles_Gini_PolySup(2); Vi_180mps_quartiles_Gini_PolySup(2); Vi_200mps_quartiles_Gini_PolySup(2)];
Q2_Gini_PolySup =  [Vi_60mps_quartiles_Gini_PolySup(3); Vi_80mps_quartiles_Gini_PolySup(3); Vi_100mps_quartiles_Gini_PolySup(3); Vi_120mps_quartiles_Gini_PolySup(3); Vi_140mps_quartiles_Gini_PolySup(3); Vi_160mps_quartiles_Gini_PolySup(3); Vi_180mps_quartiles_Gini_PolySup(3); Vi_200mps_quartiles_Gini_PolySup(3)];
Q3_Gini_PolySup =  [Vi_60mps_quartiles_Gini_PolySup(4); Vi_80mps_quartiles_Gini_PolySup(4); Vi_100mps_quartiles_Gini_PolySup(4); Vi_120mps_quartiles_Gini_PolySup(4); Vi_140mps_quartiles_Gini_PolySup(4); Vi_160mps_quartiles_Gini_PolySup(4); Vi_180mps_quartiles_Gini_PolySup(4); Vi_200mps_quartiles_Gini_PolySup(4)];
Q4_Gini_PolySup =  [Vi_60mps_quartiles_Gini_PolySup(5); Vi_80mps_quartiles_Gini_PolySup(5); Vi_100mps_quartiles_Gini_PolySup(5); Vi_120mps_quartiles_Gini_PolySup(5); Vi_140mps_quartiles_Gini_PolySup(5); Vi_160mps_quartiles_Gini_PolySup(5); Vi_180mps_quartiles_Gini_PolySup(5); Vi_200mps_quartiles_Gini_PolySup(5)];

Gini_index_all = [Gini_index_vector_Vni_60; Gini_index_vector_Vni_80; Gini_index_vector_Vni_100; Gini_index_vector_Vni_120; Gini_index_vector_Vni_140; Gini_index_vector_Vni_160; Gini_index_vector_Vni_180; Gini_index_vector_Vni_200];



% % ==============================================
% % Set constants for SPHERE geometry
% % ==============================================

R = 58; %Initial impactor radius in microns


% % ==========================
% % Load data for brittle-elastic sphere impact
% % =========================

clustersize_60mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_60_G_10\Output_files\cluster_sizes.csv");
clustersize_60mps_Sphere = sort(clustersize_60mps_Sphere);
cumulative_mass_60mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_60_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_80mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_80_G_10\Output_files\cluster_sizes.csv");
clustersize_80mps_Sphere = sort(clustersize_80mps_Sphere);
cumulative_mass_80mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_80_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_100mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_100_G_10\Output_files\cluster_sizes.csv");
clustersize_100mps_Sphere = sort(clustersize_100mps_Sphere);
cumulative_mass_100mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_100_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_120mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_120_G_10\Output_files\cluster_sizes.csv");
clustersize_120mps_Sphere = sort(clustersize_120mps_Sphere);
cumulative_mass_120mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_120_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_140mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_140_G_10\Output_files\cluster_sizes.csv");
clustersize_140mps_Sphere = sort(clustersize_140mps_Sphere);
cumulative_mass_140mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_140_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_160mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_160_G_10\Output_files\cluster_sizes.csv");
clustersize_160mps_Sphere = sort(clustersize_160mps_Sphere);
cumulative_mass_160mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_160_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_180mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_180_G_10\Output_files\cluster_sizes.csv");
clustersize_180mps_Sphere = sort(clustersize_180mps_Sphere);
cumulative_mass_180mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_180_G_10\Output_files\cumulative_volume_ratio.csv");

clustersize_200mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_200_G_10\Output_files\cluster_sizes.csv");
clustersize_200mps_Sphere = sort(clustersize_200mps_Sphere);
cumulative_mass_200mps_Sphere = readmatrix("G:\My Drive\ActivityModules\Peridynamics\LSDYNA\RESUME_Summer2024\Consider_Sphere_Geometry\CompareWithPolySupQuad_B\R_58um_dx_6um\G_10_Vary_Vni\Vi_200_G_10\Output_files\cumulative_volume_ratio.csv");






% % % --------------------------------------------------------------------
% % % -------------- Calculate t10 for SPHERE geometry -------------------
% % % --------------------------------------------------------------------

Vni_60mps_t10_Sphere=Get_tn(clustersize_60mps_Sphere,cumulative_mass_60mps_Sphere,R,0.1);
Vni_80mps_t10_Sphere=Get_tn(clustersize_80mps_Sphere,cumulative_mass_80mps_Sphere,R,0.1);
Vni_100mps_t10_Sphere=Get_tn(clustersize_100mps_Sphere,cumulative_mass_100mps_Sphere,R,0.1);
Vni_120mps_t10_Sphere=Get_tn(clustersize_120mps_Sphere,cumulative_mass_120mps_Sphere,R,0.1);
Vni_140mps_t10_Sphere=Get_tn(clustersize_140mps_Sphere,cumulative_mass_140mps_Sphere,R,0.1);
Vni_160mps_t10_Sphere=Get_tn(clustersize_160mps_Sphere,cumulative_mass_160mps_Sphere,R,0.1);
Vni_180mps_t10_Sphere=Get_tn(clustersize_180mps_Sphere,cumulative_mass_180mps_Sphere,R,0.1);
Vni_200mps_t10_Sphere=Get_tn(clustersize_200mps_Sphere,cumulative_mass_200mps_Sphere,R,0.1);

t10_values_Sphere=[Vni_60mps_t10_Sphere; Vni_80mps_t10_Sphere; Vni_100mps_t10_Sphere; Vni_120mps_t10_Sphere; Vni_140mps_t10_Sphere; Vni_160mps_t10_Sphere; Vni_180mps_t10_Sphere; Vni_200mps_t10_Sphere];




% % % ------------------------------------------------------------------
% % % -------------- Calculate Gini indices for SPHERE geometry --------
% % % ------------------------------------------------------------------

Vni_60mps_Gini_Sphere = Get_Gini_index(clustersize_60mps_Sphere);
Vni_80mps_Gini_Sphere = Get_Gini_index(clustersize_80mps_Sphere);
Vni_100mps_Gini_Sphere = Get_Gini_index(clustersize_100mps_Sphere);
Vni_120mps_Gini_Sphere = Get_Gini_index(clustersize_120mps_Sphere);
Vni_140mps_Gini_Sphere = Get_Gini_index(clustersize_140mps_Sphere);
Vni_160mps_Gini_Sphere = Get_Gini_index(clustersize_160mps_Sphere);
Vni_180mps_Gini_Sphere = Get_Gini_index(clustersize_180mps_Sphere);
Vni_200mps_Gini_Sphere = Get_Gini_index(clustersize_200mps_Sphere);

Gini_values_Sphere=[Vni_60mps_Gini_Sphere; Vni_80mps_Gini_Sphere; Vni_100mps_Gini_Sphere; Vni_120mps_Gini_Sphere; Vni_140mps_Gini_Sphere; Vni_160mps_Gini_Sphere; Vni_180mps_Gini_Sphere; Vni_200mps_Gini_Sphere];




% % % -----------------------------------------------------------------------------------------------
% % % -------------- Calculate other fracture extent measures for SPHERE geometry -------------------
% % % -----------------------------------------------------------------------------------------------

[LargestFragment_Vni_60_Sphere, NumBodyFragments_Vni_60_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_60mps_Sphere.^3);
[LargestFragment_Vni_80_Sphere, NumBodyFragments_Vni_80_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_80mps_Sphere.^3);
[LargestFragment_Vni_100_Sphere, NumBodyFragments_Vni_100_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_100mps_Sphere.^3);
[LargestFragment_Vni_120_Sphere, NumBodyFragments_Vni_120_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_120mps_Sphere.^3);
[LargestFragment_Vni_140_Sphere, NumBodyFragments_Vni_140_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_140mps_Sphere.^3);
[LargestFragment_Vni_160_Sphere, NumBodyFragments_Vni_160_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_160mps_Sphere.^3);
[LargestFragment_Vni_180_Sphere, NumBodyFragments_Vni_180_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_180mps_Sphere.^3);
[LargestFragment_Vni_200_Sphere, NumBodyFragments_Vni_200_Sphere] = Get_OtherFractureMeasures(4/3*pi*clustersize_200mps_Sphere.^3);

NumBodyFragments_values_Sphere = [NumBodyFragments_Vni_60_Sphere; NumBodyFragments_Vni_80_Sphere; NumBodyFragments_Vni_100_Sphere; NumBodyFragments_Vni_120_Sphere; NumBodyFragments_Vni_140_Sphere; NumBodyFragments_Vni_160_Sphere; NumBodyFragments_Vni_180_Sphere; NumBodyFragments_Vni_200_Sphere];
LargestFragment_values_Sphere = [LargestFragment_Vni_60_Sphere; LargestFragment_Vni_80_Sphere; LargestFragment_Vni_100_Sphere; LargestFragment_Vni_120_Sphere; LargestFragment_Vni_140_Sphere; LargestFragment_Vni_160_Sphere; LargestFragment_Vni_180_Sphere; LargestFragment_Vni_200_Sphere];









% % % MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
% % % ======================================================
% % % Create manuscript plots
% % % =====================================================
% % % MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM




% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% % % % % -------------- VISUALIZE: natural Polysuperquadric B fracture measures relative to t_10, Gini index.
% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% 
% fullmap = turbo(256);
% LargestFragments_all_dim = (3/(4*pi)*(actual_particle_volume*LargestFragments_all)).^(1/3)*1000;
% scatter(Gini_index_all, LargestFragments_all_dim, 50, fullmap(100,:), 'filled', 'o');
% hold on;
% scatter(t10_all, LargestFragments_all_dim, 50, fullmap(16,:), 'filled', 'o');
% 
% set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'Box', 'off');
% lgd = legend('$\mathcal{D}$','$\mathrm{t}_{10}$', 'Interpreter', 'latex', 'Box','off','FontSize',20);
% ylabel('Largest fragment','interpreter','latex','FontSize',16);
% set(gcf, 'Color', 'w')  % Set background to white
% % exportgraphics(gcf, 'D_vs_t10.png', 'ContentType', 'auto', 'Resolution', 300)


 






% % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% % % % -------------- VISUALIZE: Polysuperquadric B t_10 distribution relative to sphere over variable impact speeds: quantile interval curves  -------
% % % % ------------------------------------------------------------------------------------------------------------------------------------------------------

fullmap = turbo(256);

V = [60; 80; 100; 120; 140; 160; 180; 200];
% % % polysuperquadric B result
% Shade the region between the 25% and 75% quartiles
x = [V; flip(V)]; % x-coordinates: forward then backward
y = [Q1_t10_PolySup; flip(Q3_t10_PolySup)]; % y-coordinates: lower bound then upper bound
hold on
color = fullmap(100,:);
fill(x, y, color, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % light blue fill with no edge
% Then overlay the curves
% plot(V,Q0_t10_PolySup,'b-','LineWidth',0.5,'MarkerFaceColor','b');
plot(V,Q1_t10_PolySup,'-','Color',color,'LineWidth',2,'MarkerFaceColor',color);
plot(V,Q2_t10_PolySup,'o-','Color',color,'LineWidth',3,'MarkerFaceColor',color);
plot(V,Q3_t10_PolySup,'-','Color',color,'LineWidth',2,'MarkerFaceColor',color);
% plot(V,Q4_t10_PolySup,'b-','LineWidth',0.5,'MarkerFaceColor','b');

% % % plot brittle-elastic sphere result for comparison
color = fullmap(16,:);
V = [60; 80; 100; 120; 140; 160; 180; 200];
plot(V, t10_values_Sphere,'o','MarkerFaceColor',color,'Color',color,'linewidth',3);

% lgd = legend(l,'$25\%  / 75\%$','$50\%$','interpreter','latex','FontSize',14);
% title(lgd,'mass-weighted quartiles','interpreter','latex','FontSize',15);

% % % % Plot markers corresponding to orientations visualized
% V_200mps = [200; 200; 200];
% t_10_visualized_200mps = [0.2; 0.3; 0.5];
% scatter(V_200mps, t_10_visualized_200mps,200,'k','x','linewidth',2);
% 
% V_100mps = [100; 100; 100];
% t_10_visualized_100mps = [0.05; 0.09; 0.13];
% scatter(V_100mps, t_10_visualized_100mps,200,'k','x','linewidth',2);

ylim([0 0.5])

set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'Box', 'off');
xlabel('$V_{\mathrm{ni}}$ $\left(\mathrm{m/s}\right)$','interpreter','latex','FontSize',20);
ylabel('$\mathrm{t}_{10}$','interpreter','latex','FontSize',26);
set(gcf, 'Color', 'w')  % Set background to white
% exportgraphics(gcf, 't10_sphere_vs_PolySup.png', 'ContentType', 'auto', 'Resolution', 300)










% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% % % % % -------------- VISUALIZE: Polysuperquadric B Gini index distribution relative to sphere over variable impact speeds: quantile interval curves  -------
% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% 
% fullmap = turbo(256);
% 
% V = [60; 80; 100; 120; 140; 160; 180; 200];
% % % % polysuperquadric B result
% % Shade the region between the 25% and 75% quartiles
% x = [V; flip(V)]; % x-coordinates: forward then backward
% y = [Q1_Gini_PolySup; flip(Q3_Gini_PolySup)]; % y-coordinates: lower bound then upper bound
% hold on
% color = fullmap(100,:);
% fill(x, y, color, 'FaceAlpha', 0.1, 'EdgeColor', 'none'); % light blue fill with no edge
% % Then overlay the curves
% % plot(V,Q0_Gini_PolySup,'b-','LineWidth',0.5,'MarkerFaceColor','b');
% plot(V,Q1_Gini_PolySup,'-','Color',color,'LineWidth',2,'MarkerFaceColor', color);
% plot(V,Q2_Gini_PolySup,'o-','Color',color,'LineWidth',3,'MarkerFaceColor', color);
% plot(V,Q3_Gini_PolySup,'-','Color',color,'LineWidth',2,'MarkerFaceColor', color);
% % plot(V,Q4_Gini_PolySup,'b-','LineWidth',0.5,'MarkerFaceColor','b');
% 
% 
% % % % plot brittle-elastic sphere result for comparison
% 
% color = fullmap(16,:);
% V = [60; 80; 100; 120; 140; 160; 180; 200];
% plot(V, Gini_values_Sphere,'o','MarkerFaceColor',color,'Color',color,'linewidth',3);
% 
% % lgd = legend(l,'$25\%  / 75\%$','$50\%$','interpreter','latex','FontSize',14);
% % title(lgd,'mass-weighted quartiles','interpreter','latex','FontSize',15);
% 
% % % % % Plot markers corresponding to orientations visualized
% % V_200mps = [200; 200; 200];
% % D_visualized_200mps = [0.87; 0.94; 0.99];
% % scatter(V_200mps, D_visualized_200mps,200,'k','x','linewidth',2);
% % 
% % V_100mps = [100; 100; 100];
% % D_visualized_100mps = [0.23; 0.58; 0.76];
% % scatter(V_100mps, D_visualized_100mps,200,'k','x','linewidth',2);
% 
% ylim([0 1])
% 
% set(gca,'FontSize',14, 'TickLabelInterpreter', 'latex', 'TickDir', 'out', 'Box', 'off');
% xlabel('$V_{\mathrm{ni}}$ $\left(\mathrm{m/s}\right)$','interpreter','latex','FontSize',20);
% ylabel('$\mathcal{D}$','interpreter','latex','FontSize',26);
% set(gcf, 'Color', 'w')  % Set background to white
% % exportgraphics(gcf, 'D_sphere_vs_PolySup.png', 'ContentType', 'auto', 'Resolution', 300)








% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% % % % % -------------- VISUALIZE: natural Polysuperquadric B fracture measures relative to t_10, Gini index, and entropy.
% % % % % ------------------------------------------------------------------------------------------------------------------------------------------------------
% 
% scatter(Gini_index_all, LargestFragments_all, 50, 'k', 'filled', 'o');
% hold on;
% scatter(t10_all, LargestFragments_all, 50, 'r', 'filled', 'o');
% legend('$\mathcal{D}$','$\mathrm{t}_{10}$','interpreter','latex','FontSize',16);
% ylabel('Normalized, largest fragment','interpreter','latex','FontSize',16);



% scatter(Gini_index_all, NumBodyFragments_all, 50, 'k', 'filled', 'o');
% hold on;
% scatter(t10_all, NumBodyFragments_all, 50, 'r', 'filled', 'o');
% legend('$\mathcal{D}$','$\mathrm{t}_{10}$','interpreter','latex','FontSize',16);
% ylabel('Fragment count','interpreter','latex','FontSize',16);


% scatter(Gini_index_all, LargestFragments_all, 50, 'k', 'filled', 'o');
% hold on;
% scatter(t10_all, LargestFragments_all, 50, 'r', 'filled', 'o');
% scatter(Entropy_all./max(Entropy_all), LargestFragments_all, 50, 'b', 'filled', 'o');
% legend('$\mathcal{D}$','$\mathrm{t}_{10}$','interpreter','latex','FontSize',16);
% ylabel('Normalized, largest fragment','interpreter','latex','FontSize',16);


% scatter(Gini_index_all, NumBodyFragments_all, 50, 'k', 'filled', 'o');
% hold on;
% scatter(t10_all, NumBodyFragments_all, 50, 'r', 'filled', 'o');
% % scatter(Entropy_all./max(Entropy_all), NumBodyFragments_all, 50, 'b', 'filled', 'o');
% legend('$\mathcal{D}$','$\mathrm{t}_{10}$','interpreter','latex','FontSize',16);
% ylabel('Fragment count','interpreter','latex','FontSize',16);











function [t10_vector, Orientations, NumBodyFragments_vector, LargestFragment_vector] = Extract_t10_set(baseDir)


    % % % input arguments
    % % % ------------------------
    % % % baseDir: path location to the folder collected_output which stores all output for each impact condition.

    % Set constants and parameters
    actual_particle_volume=0.000826242;
    R = (3/(4*pi)*actual_particle_volume)^(1/3)*1000; %Initial impactor radius in microns
    rho = 2.44000E-6; %particle mass density [kg/mm^3]
    folders = dir(fullfile(baseDir, 'MAIN_*_*_*'));

    % initialize vectors
    t10_vector = nan(length(folders),1);
    ClusterVolumeCheck = zeros(length(folders),1);
    SizeDist_cell = cell(1,length(folders));
    Orientations = nan(length(folders),3);
    NumBodyFragments_vector = nan(length(folders),1);
    LargestFragment_vector = nan(length(folders),1);


    for i = 1:length(folders)

        folderPath = fullfile(folders(i).folder, folders(i).name);
        disp(i);

        folderName = folders(i).name;  % e.g., 'MAIN_5p00_15p00_0p00'
        % Extract parts after 'MAIN_'
        parts = split(folderName, '_');  % {'MAIN', '5p00', '15p00', '0p00'}
        % Replace 'p' with '.' and convert to numeric
        theta_1 = str2double(strrep(parts{2}, 'p', '.'));
        theta_2 = str2double(strrep(parts{3}, 'p', '.'));
        theta_3 = str2double(strrep(parts{4}, 'p', '.'));

        clusterVolumesFile = fullfile(folderPath, 'clusterVolumes.csv');
        clusterVolumes = readmatrix(clusterVolumesFile);
        totalClusterVolume=sum(clusterVolumes);

        % Extract metric for the extent of "body-level" particle breakage
        [LargestFragment, NumBodyFragments] = Get_OtherFractureMeasures(clusterVolumes);
        LargestFragment_vector(i) = LargestFragment;
        NumBodyFragments_vector(i) = NumBodyFragments;

        clusterSizesFile = fullfile(folderPath, 'cluster_sizes.csv');
        ClusterSizes = readmatrix(clusterSizesFile);
        ClusterSizes = sort(ClusterSizes);

        CumulativeMassFile = fullfile(folderPath, 'cumulative_volume_ratio.csv');
        CumulativeMass = readmatrix(CumulativeMassFile);

        t10 = Get_tn(ClusterSizes, CumulativeMass, R, 0.1);

        % Store results
        Orientations(i,:) = [theta_1 theta_2 theta_3];
        ClusterVolumeCheck(i)=totalClusterVolume/actual_particle_volume;
        SizeDist_cell{i} = [ClusterSizes./R,CumulativeMass];
        t10_vector(i) = t10;

    end

end











function [Gini_vector, Orientations] = Extract_GiniIndex_set(baseDir)

    
    % % % input arguments
    % % % ------------------------
    % % % baseDir: path location to the folder collected_output which stores all output for each impact condition.

    % Set constants and parameters
    actual_particle_volume=0.000826242;
    R = (3/(4*pi)*actual_particle_volume)^(1/3)*1000; %Initial impactor radius in microns
    rho = 2.44000E-6; %particle mass density [kg/mm^3]
    folders = dir(fullfile(baseDir, 'MAIN_*_*_*'));
    
    % initialize vectors
    Gini_vector = nan(length(folders),1);
    ClusterVolumeCheck = zeros(length(folders),1);
    SizeDist_cell = cell(1,length(folders));
    Orientations = nan(length(folders),3);
    
    
    for i = 1:length(folders)
    
        folderPath = fullfile(folders(i).folder, folders(i).name);
    
        folderName = folders(i).name;  % e.g., 'MAIN_5p00_15p00_0p00'
        % Extract parts after 'MAIN_'
        parts = split(folderName, '_');  % {'MAIN', '5p00', '15p00', '0p00'}
        % Replace 'p' with '.' and convert to numeric
        theta_1 = str2double(strrep(parts{2}, 'p', '.'));
        theta_2 = str2double(strrep(parts{3}, 'p', '.'));
        theta_3 = str2double(strrep(parts{4}, 'p', '.'));
    
        clusterVolumesFile = fullfile(folderPath, 'clusterVolumes.csv');
        clusterVolumes = readmatrix(clusterVolumesFile);
        totalClusterVolume=sum(clusterVolumes);
    
        clusterSizesFile = fullfile(folderPath, 'cluster_sizes.csv');
        ClusterSizes = readmatrix(clusterSizesFile);
        ClusterSizes = sort(ClusterSizes);
    
        CumulativeMassFile = fullfile(folderPath, 'cumulative_volume_ratio.csv');
        CumulativeMass = readmatrix(CumulativeMassFile);
       
        Gini_index = Get_Gini_index(ClusterSizes);
    
        % Store results
        Orientations(i,:) = [theta_1 theta_2 theta_3];
        ClusterVolumeCheck(i)=totalClusterVolume/actual_particle_volume;
        SizeDist_cell{i} = [ClusterSizes./R,CumulativeMass];
        Gini_vector(i) = Gini_index;
        
    end

end












function [LargestFragment, NumBodyFragments] = Get_OtherFractureMeasures(clusterVolumes)
    totalClusterVolume=sum(clusterVolumes);
    % Extract metric for the extent of "body-level" particle breakage
    LargestFragmentVolume = max(clusterVolumes);
    LargestFragment = LargestFragmentVolume / totalClusterVolume;
    NumBodyFragments = nnz(clusterVolumes > (0.1).*LargestFragmentVolume);
end




function tn = Get_tn(ClusterSizes, CumulativeMass, R, n)
    % Extract t_n as the cumulative mass corresponding to the closest normalized cluster size to n that is less than n.
    Normed_clusterSizes = ClusterSizes./R;
    SmallestFragment_indices = find(Normed_clusterSizes <= n); % Find indices where x is less than C
    if isempty(SmallestFragment_indices)
        tn = 0; % Handle case where no values satisfy the condition
    else
        [~, SmallestList_idx] = min(abs(Normed_clusterSizes(SmallestFragment_indices) - n)); % Find index of closest value within valid indices
        closest_index = SmallestFragment_indices(SmallestList_idx); % Get original index in x
        tn = CumulativeMass(closest_index); % Access corresponding y value
    end
end



function Gini_index = Get_Gini_index(ClusterSizes_sorted)
    % Gini (mass distribution) index
    ClusterVolumes_sorted = (4/3)*pi*ClusterSizes_sorted.^3;
    total_volume = sum(ClusterVolumes_sorted);
    Gini_index = 1 - sum((ClusterVolumes_sorted./total_volume).^(2));
end





function quartiles = GetQuartiles(data)

    if length(data) < 4
        
        q2 = mean(data);
        q1 = q2;
        q3 = q2;
        lower_whisker = q2;
        upper_whisker = q2;
        quartiles = [lower_whisker; q1; q2; q3; upper_whisker];

    else

        sorted_data = sort(data);
        quantile_levels = [0.15, 0.25, 0.50, 0.75, 0.85];
        % Calculate the quartiles using interpolation
        quartiles = interp1(linspace(0,1,length(sorted_data)), sorted_data, quantile_levels).';

    end

end






function Gini_index_vector_Vni_100_reordered = reorder_gini_vector(...
    Gini_index_vector_Vni_100, Orientations_Vni_100_Gini, Orientations_Vni_80_Gini)

    % Ensure inputs are of compatible sizes
    assert(size(Orientations_Vni_100_Gini,2) == 3 && size(Orientations_Vni_80_Gini,2) == 3, ...
        'Orientation matrices must have 3 columns.');
    assert(size(Orientations_Vni_100_Gini,1) == length(Gini_index_vector_Vni_100), ...
        'Number of Gini index entries must match number of orientations in Vni_100.');

    N = size(Orientations_Vni_80_Gini, 1);
    reordered_Gini_100 = nan(N, 1);
    unmatched_Gini_100 = [];

    % Convert to string keys for comparison
    str_80 = string(Orientations_Vni_80_Gini);
    str_100 = string(Orientations_Vni_100_Gini);
    keys_80 = join(str_80, ",", 2);
    keys_100 = join(str_100, ",", 2);

    % Create map from key to Gini index
    keyMap = containers.Map(keys_100, Gini_index_vector_Vni_100);

    % Fill matched entries
    for i = 1:N
        key = keys_80(i);
        if isKey(keyMap, key)
            reordered_Gini_100(i) = keyMap(key);
        end
    end

    % Add unmatched entries at the end
    unmatched_keys = setdiff(keys_100, keys_80);
    for k = 1:length(unmatched_keys)
        unmatched_Gini_100(k,1) = keyMap(unmatched_keys(k));
    end

    Gini_index_vector_Vni_100_reordered = [reordered_Gini_100; unmatched_Gini_100];
end








function reordered_matrix = Reorder_matrix_by_orientation(input_matrix, target_orientation_order)
% Reorder_matrix_by_orientation
% Reorders the rows of IE_matrix so that columns 2–4 align with target_orientations.
%
% Inputs:
%   input_matrix           - NxM matrix, where columns 2 to 4 represent orientations
%   target_orientation_order - Nx3 matrix, each row is a [theta1, theta2, theta3] triplet
%
% Output:
%   reordered_matrix - NxM matrix, reordered to match target_orientations

    IE_orientations = input_matrix(:, 2:4);
    num_targets = size(target_orientation_order, 1);
    reordered_matrix = nan(num_targets, size(input_matrix,2));

    for i = 1:num_targets
        idx = find(ismember(IE_orientations, target_orientation_order(i,:), 'rows'), 1);
        if ~isempty(idx)
            reordered_matrix(i, :) = input_matrix(idx, :);
        else
            warning("No match found for row %d in target_orientations.", i);
        end
    end
end



