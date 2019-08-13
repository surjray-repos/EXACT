%% example of how to use CITUP

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% parameters
wrapper_working_directory = [pwd,'/distribution/'];
CITUP_executable_path = [pwd,'/distribution/bin/'];
min_cluster_no = 5; % here we make the number of clusters be equal to the size of the input tree
max_cluster_no = 9;
citup_error_rate = 0.03;

% run Citup
[~, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate );

% get output
Tree_Matrix_T = citup_output{1};
Mutant_Frequencies_M = citup_output{3}; 