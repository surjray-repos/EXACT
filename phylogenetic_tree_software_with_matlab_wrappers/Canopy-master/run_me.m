%% example of how to use Canopy

% input file
input_file = '/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_4.input';
path_to_folder = [pwd, '/demo_code/'];

% parameters
burnin_val = 10;
thin_val = 5;
K_min_val = 3;
K_max_val = 3;
numchains_val = 15;
maxsimrun_val = 100000;
minsimrun_val = 10000;
writeskip_val = 200;
cluster_number_start = 2;
cluster_number_end = 9;

% run canopy
canopy_output = canopy_wrapper(input_file, path_to_folder, burnin_val, thin_val, K_min_val, K_max_val, numchains_val, maxsimrun_val, minsimrun_val, writeskip_val, cluster_number_start, cluster_number_end);

% get output
Ancestry_Matrix_U = canopy_output{1};
Mutant_Frequencies_M = canopy_output{2};
