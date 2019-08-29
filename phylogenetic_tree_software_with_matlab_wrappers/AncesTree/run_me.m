%% example of how to use Ancestree

% first install CVX, if it is not install yet
% we recommend downloading, and installing, the latest CVX version from
% www.cvxr.com/cvx
run([pwd, '/../EXACT/distribution/cvx/cvx/cvx_setup']);

% adding path for the AncesTree folders
addpath(genpath(pwd));

% to find functions to produce plots and compute errors
addpath([pwd, '/../all_tests/']);

% input file
pwd_AncesTree = pwd;
input_file = [pwd_AncesTree, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_1.input'];
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = [pwd_AncesTree, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_1.true'];
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};

% parameters
n_mutations = size(F_from_SampleData,1);
T_samples = size(F_from_SampleData,2);
disp(n_mutations);
disp(T_samples);

scale = scaling;
alpha = 0.3; % if alpha is big, lots of things will be clustered together
beta = 0.8; % to choose a larger beta we need more samples , i.e. larger T_samples
gamma = 0.01; % small gamma means larger confidence on the data
wrapper_working_directory = [pwd_AncesTree, '/distribution/'];
ancestree_executable_path = [pwd_AncesTree, '/distribution/build/ancestree'];

% run AncesTree
ancestree_output = ancestree_wrapper(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);

% get output
% since AncesTree might have multiple solutions, we shall check how many
% solutions are inferred and reference them by sol_id
sol_id = 1; %Considering the first output tree
Tree_Matrix_T = ancestree_output{3}{sol_id};
Mutant_Frequencies_M = ancestree_output{2}{sol_id}';

U = ancestree_output{3}{sol_id};
U = inv(eye(length(U)) - U);
clust = [];
for i = 1:length(ancestree_output{4}{sol_id})
	for j = ancestree_output{4}{sol_id}{i}'
		clust = [clust; [j,i]];
	end
end

% generate figure with an optimal tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);

% generate figure with all the output trees inferred by AncesTree, and their fitting costs
generate_list_of_AncesTree_trees_with_fitting_costs(ancestree_output, F_from_SampleData);