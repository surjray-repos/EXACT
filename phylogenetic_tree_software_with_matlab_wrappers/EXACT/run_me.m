%% example of how to use our brute force phylogeny code EXACT!

% first install CVX, if it is not install yet
% we recommend downloading, and installing, the latest CVX version from
% www.cvxr.com/cvx
run([pwd, '/distribution/cvx/cvx/cvx_setup']);

% to find functions to produce plots and compute errors
addpath('../../all_tests/');

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true';
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};


% parameters
path_to_folder = [pwd, '/distribution/'];
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
[ourcode_output] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);

% get output
% ourcode_output is a matlab cell object having 7 components, namely, 
% where
%	ourcode_output{1} = likelihood score of the best tree as computed by the BIC criterion
%	ourcode_output{2} = Bayesian information criteria score
%	ourcode_output{3} = adjacency matrix for the best tree. This is a directed tree. If we can this matrix T, the U = inv(I - T), where U appears in the PPM model as F = UM.
%	ourcode_output{4} = recovered (clean) frequencies of mutations
%	ourcode_output{5} = clustered frequencies of mutants
%	ourcode_output{6} = cluster membership information, the number in column 1 of each row designating which cluster/node each mutation (row number) belongs to.
%	ourcode_output{7} = run time (in seconds) for the executable to infer the best tree among all trees of the same size while keeping track of all k_best trees
Tree_Matrix_T = ourcode_output{3};
Mutant_Frequencies_M = ourcode_output{4}; 

U = inv(eye(size(ourcode_output{3})) - ourcode_output{3});
clust = [ [1: length(ourcode_output{6})]' , ourcode_output{6}];

generate_joint_plot_of_tree_donut_muller_and_errors(U,clust, Mutant_Frequencies_M,  Ugt, clustgt);

%% generate a list of top 10 trees

