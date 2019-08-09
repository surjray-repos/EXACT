%% example of how to use our brute force phylogeny code EXACT!

% first install CVX, if it is not install yet
% we recommend downloading, and installing, the latest CVX version from
% www.cvxr.com/cvx
run([pwd, '/distribution/cvx/cvx/cvx_setup']);

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

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
Tree_Matrix_T = ourcode_output{3};
Mutant_Frequencies_M = ourcode_output{4}; 