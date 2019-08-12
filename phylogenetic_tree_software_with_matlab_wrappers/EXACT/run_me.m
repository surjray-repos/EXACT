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

% build tree labels from the cluster information
nodelbs = cell(1,length(U));
for i = 1:length(U)
	nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
end

% compute colors for each mutant
node_col = cell(1,length(U));
cum_nodelbs = cell(1,length(U));
for i = 1:length(U)
	for j = find(U(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	
	rng(mod(  sum(((length(clust)).^[1:length(array_rep)]).*array_rep )    ,(2^32)-1)); % use a simple hash to generate a seed to then generate a random color for each mutant
	clrix = randi(64);
	node_col{i} = colormap_array(clrix,:); 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end
%% Tree diagram with information about which nodes contain which mutations
figure;
subplot(1,3,1);
h = plot(digraph( eye(length(U)) - inv(U) ));
labelnode(h,[1:length(U)],nodelbs);
highlight(h,[1:length(U)],'MarkerSize',20);
for i = 1:length(U)
	highlight(h,i,'NodeColor',node_col{i});
end
set(gca,'visible','off');
%% draw donut plot with mutant mixing ratios
subplot(1,3,2);
[~, lgd] = donut(Mutant_Frequencies_M(:,1:5)', cum_nodelbs, node_col);
set(lgd,'visible','off');
set(gca,'visible','off');

%% assuming that the different samples are obtained in time, we can draw a muller plot
subplot(1,3,3);
[~, lgd] = generate_simple_muller_plots(U, Mutant_Frequencies_M, cum_nodelbs, node_col);
set(lgd,'visible','off');
set(gca,'visible','off');

%% compute four different error types between synthetic data ground truth and EXACT inferred tree
%See section 4.1 in paper
ground_truth_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true';
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
U1 = true_tree_data{3}';
clust1 = true_tree_data{5};

U2 = U;
clust2 = clust;

%Calculate the 4 different error types 
[error_rates_EXACT] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);
figure;
bar(error_rates_EXACT);
set(gca, 'xticklabel', {'error1', 'error2', 'error3', 'error4'})