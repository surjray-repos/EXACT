%% example of how to use PhyloWGS

% to find functions to produce plots and compute errors
addpath('../../all_tests/');

% input file
pwd_PhyloWGS = pwd;
input_file = [pwd_PhyloWGS, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input'];
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = [pwd_PhyloWGS, '/../all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true'];
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};

% parameters
scale = scaling;

%Constants for PhyloWGS, should not be changed for sequencing data
mu_r = 0.999;
mu_v = 0.5;
wrapper_working_directory = [pwd_PhyloWGS,'/distribution/'];
phylowgs_executable_path = [pwd_PhyloWGS,'/distribution/'];

% run PhyloWGS
[phylowgs_output] = PhyloWGS_wrapper(F_from_SampleData, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);

% get output
% since PhyloWGS might have multiple solutions, we shall check how many
% solutions are inferred and reference them by sol_id
sol_id = 1; %Considering the first output tree
Tree_Matrix_T = phylowgs_output{1}{sol_id};

U = phylowgs_output{1}{sol_id};
root_phylowgs = 1; % the root of the treee is always node 1.
n_virt_nodes_phylowgs = size(U,1);
AdjLrecon = {};
for j =1:n_virt_nodes_phylowgs
	AdjLrecon{j} = find(U(:,j));
end
Treerecon = BFS(AdjLrecon,root_phylowgs); %root is the node 0 in unlabled tree
% get adj mat of tree
AdjTrecon = zeros(n_virt_nodes_phylowgs);
for j =1:n_virt_nodes_phylowgs
	for r = Treerecon{j}
		AdjTrecon(j,r) = 1;
	end
end
U = inv(eye(n_virt_nodes_phylowgs) - AdjTrecon);
num_mutants = size(phylowgs_output{2}{sol_id},1);

clust = zeros(num_mutants,2);
for i = 1:size(phylowgs_output{2}{sol_id},1)
	clust(1 + phylowgs_output{2}{sol_id}{i,1},:) = [1 + phylowgs_output{2}{sol_id}{i,1},1 + phylowgs_output{2}{sol_id}{i,2}];
end

Mutant_Frequencies_M = inv(U)*phylowgs_output{3}{sol_id};

% generate figure with an optimal tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);

% generate figure with all the output trees inferred by PhyloWGS
% we note that, often, many of these trees are actually the same
generate_list_of_PhyloWGS_trees(phylowgs_output);