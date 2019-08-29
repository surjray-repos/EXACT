%% example of how to use Canopy

% first install CVX, if it is not install yet
% we recommend downloading, and installing, the latest CVX version from
% www.cvxr.com/cvx
run([pwd, '/../EXACT/distribution/cvx/cvx/cvx_setup']);

% adding path for the Canopy folders
addpath(genpath(pwd));

% to find functions to produce plots and compute errors
addpath([pwd, '/../all_tests/']);

% input file
pwd_Canopy = pwd;
input_file = [pwd_Canopy, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input'];
path_to_folder = [pwd_Canopy, '/demo_code/'];

% load groud truth
ground_truth_file = [pwd_Canopy, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true'];
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};

% parameters
burnin_val = 10;
thin_val = 5;
K_min_val = 3;
K_max_val = 6;
numchains_val = 15;
maxsimrun_val = 100000;
minsimrun_val = 10000;
writeskip_val = 200;
cluster_number_start = 2;
cluster_number_end = 9;

% run canopy
canopy_output = canopy_wrapper(input_file, path_to_folder, burnin_val, thin_val, K_min_val, K_max_val, numchains_val, maxsimrun_val, minsimrun_val, writeskip_val, cluster_number_start, cluster_number_end);

% get output
Mutant_Frequencies_M = canopy_output{2};

% this function reads read the output from canopy and extracts a tree and a
% clustering
[U, clust] = extract_U_mat_and_clust_from_canopy_output(canopy_output);

% generate figure with an optimal tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);