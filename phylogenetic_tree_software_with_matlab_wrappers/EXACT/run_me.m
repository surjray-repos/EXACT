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
[ourcode_output, ~, ourcode_all_Ms] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);

Tree_Matrix_T = ourcode_output{3};
Mutant_Frequencies_M = ourcode_output{4}; 

U = inv(eye(size(ourcode_output{3})) - ourcode_output{3});
clust = [ [1: length(ourcode_output{6})]' , ourcode_output{6}];

% generate figure with best tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);

% generate a list of top 10 trees ranked by their probabilities according to our assumed model

%%
figure;

tree_probs = [];
sorted_ids = nan(top_k_value, 1 + max_tree_size - min_tree_size);
sorted_vals = nan(top_k_value, 1 + max_tree_size - min_tree_size);

for j = 1:1 + max_tree_size - min_tree_size
	for i = 1:top_k_value
		sorted_vals(i,j) = ourcode_all_Ms{j}{i}{1};
	end
	[sorted_vals(:,j), sorted_ids(:,j)] = sort(sorted_vals(:,j));
end

for t = 1:10 % this will show us the top 10 trees
	best_size = size(ourcode_output{3}, 1); % we can choose which size to work with here
	best_tree_ix = sorted_ids(t, best_size - min_tree_size);
	best_ourcode_output = ourcode_all_Ms{best_size - min_tree_size}{best_tree_ix};

	% the probability decays fast, so that we can normalized by the sum over the first few trees
	tree_probs(t) = (exp(-best_ourcode_output{1}./(2*error_rate*error_rate)))/sum(exp(-sorted_vals(:,best_size - min_tree_size)./(2*error_rate*error_rate)));

	M_target = best_ourcode_output{4};

	U2 = inv(eye(size(best_ourcode_output{3})) - best_ourcode_output{3});
	clust2 = [ [1: length(best_ourcode_output{6})]' , best_ourcode_output{6}];

	nodelbs = cell(1,length(U2));
	for i = 1:length(U2)
		nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
	end

	ax = subplot(2,5,t);
	h = plot(digraph(  eye(length(U2)) - inv(U2)      ));
	labelnode(h, [1:length(U2)],nodelbs);
	if (t==1)
		title(ax, ['List of probabilities: ' num2str(tree_probs(t))]);
	else
		title(ax, num2str(tree_probs(t)));
	end
	set(ax,'visible','off');
	set(findall(ax, 'type', 'text'), 'visible', 'on');

end
	
% maximize window
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);