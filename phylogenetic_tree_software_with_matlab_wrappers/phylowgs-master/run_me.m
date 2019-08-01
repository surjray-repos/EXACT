%% example of how to use PhyloWGS

% input file
input_file = '/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% parameters
scale = scaling;

%Constants for PhyloWGS, should not be changed for sequencing data
mu_r = 0.999;
mu_v = 0.5;
wrapper_working_directory = [pwd,'/distribution/'];
phylowgs_executable_path = [pwd,'/distribution/'];

% run PhyloWGS
[phylowgs_output] = PhyloWGS_wrapper(F_from_SampleData, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);

% get output
Tree_Matrix_T = phylowgs_output{1};