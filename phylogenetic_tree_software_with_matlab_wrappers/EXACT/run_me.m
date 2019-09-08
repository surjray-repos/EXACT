%% example of how to use EXACT

% first install CVX, if it is not install yet
% we recommend downloading, and installing, the latest CVX version from
% www.cvxr.com/cvx
run([pwd, '/distribution/cvx/cvx/cvx_setup']);

% to find functions to produce plots and compute errors
addpath('../../all_tests/');

% adding path for the EXACT folders
addpath(genpath(pwd));

% input file
pwd_EXACT = pwd;
input_file = [pwd_EXACT, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input'];
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = [pwd_EXACT, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true'];
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};

% EXACT parameters
path_to_folder = [pwd_EXACT, '/distribution/'];
exec_name = 'EXACT_executable_x64_CUDA.out';
cost_function = 'cost1'; %possible options: cost1, cost2, cost3, cost4
cpu_gpu = 'cpu_multithread'; %possible options: cpu, cpu_multithread, gpu
GPU_id = 0;
min_tree_size = 6;
max_tree_size = 9;
num_CPU_cores = 50;
error_rate = 0.03;
max_num_partitions = 1;
device_tree_subset_value = 1;
CUDA_threads_per_block = 32;
CUDA_blocks = 128;
top_k_value = 20;

% run EXACT
[ourcode_output, ~, ourcode_all_Ms] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);

% extra Tree matrix and clustering
Tree_Matrix_T = ourcode_output{3};
Mutant_Frequencies_M = ourcode_output{4}; 

U = inv(eye(size(ourcode_output{3})) - ourcode_output{3});
clust = [ [1: length(ourcode_output{6})]' , ourcode_output{6}];

% generate figure with best tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);

% generate a list of top 10 trees ranked by their probabilities according to our assumed model
generate_list_of_top_k_trees_with_probabilities(10, top_k_value, max_tree_size, min_tree_size, error_rate, ourcode_all_Ms, ourcode_output);