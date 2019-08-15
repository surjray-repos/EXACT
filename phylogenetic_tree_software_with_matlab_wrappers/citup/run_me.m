%% example of how to use CITUP

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_6_Mut_100_Clone_10_PCR_Removed/sim_4.true';
[true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
Ugt = true_tree_data{3}';
clustgt = true_tree_data{5};

% parameters
wrapper_working_directory = [pwd,'/distribution/'];
CITUP_executable_path = [pwd,'/distribution/bin/'];
min_cluster_no = 5; % here we make the number of clusters be equal to the size of the input tree
max_cluster_no = 8;
citup_error_rate = 0.03;

% run Citup
[~, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate );

% get output
% since CITUP might have multiple solutions, we shall check how many
% solutions are inferred and reference them by sol_id
tree_size_id = 1; %Considering the first output tree
Tree_Matrix_T = citup_output{1}{tree_size_id};
Mutant_Frequencies_M = citup_output{3}{tree_size_id}';

U = citup_output{1}{tree_size_id};
root_citup = 1; % the root of the treee is always node 1.
n_virt_nodes_citup = size(U,1);
AdjLrecon = {};
for j =1:n_virt_nodes_citup
	AdjLrecon{j} = find(U(:,j));
end
Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
% get adj mat of tree
AdjTrecon = zeros(n_virt_nodes_citup);
for j =1:n_virt_nodes_citup
	for r = Treerecon{j}
		AdjTrecon(j,r) = 1;
	end
end
U = inv(eye(n_virt_nodes_citup) - AdjTrecon);
clust = citup_output{2};
clust(:,2) = 1 + clust(:,2);

% generate figure with an optimal tree, mutant frequencies, and errors to ground truth
generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M,  Ugt, clustgt);

% generate figure with all the output trees inferred by CITUP, and their fitting costs
figure;

num_sizes_tested = size(all_inter_sol, 2);

for tree_size_id = 1:num_sizes_tested % this will show us all output trees
	
	U = all_inter_sol{tree_size_id}{1}{1};
	root_citup = 1; % the root of the treee is always node 1.
	n_virt_nodes_citup = size(U,1);
	AdjLrecon = {};
	for j =1:n_virt_nodes_citup
		AdjLrecon{j} = find(U(:,j));
	end
	Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
	% get adj mat of tree
	AdjTrecon = zeros(n_virt_nodes_citup);
	for j =1:n_virt_nodes_citup
		for r = Treerecon{j}
			AdjTrecon(j,r) = 1;
		end
	end
	U = inv(eye(n_virt_nodes_citup) - AdjTrecon);
	clust = all_inter_sol{tree_size_id}{2};
	clust(:,2) = 1 + clust(:,2);
	
	cost = norm(all_inter_sol{tree_size_id}{5}' - all_inter_sol{tree_size_id}{4}{1} ,'fro')^2; 
	
	nodelbs = cell(1,length(U));
	for i = 1:length(U)
		nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
	end
	
	ax = subplot(2, num_sizes_tested/2, tree_size_id);
	h = plot(digraph(  eye(length(U)) - inv(U)  ));
	labelnode(h, [1:length(U)],nodelbs);
	if (tree_size_id == 1)
		title(ax, ['CITUP (fitting cost, tree size) for the best tree for each size : (' , num2str(cost), ',', num2str(min_cluster_no + tree_size_id - 1),')']);
	else
		title(ax,['(' , num2str(cost), ',', num2str(min_cluster_no + tree_size_id - 1),')']);
	end
	set(ax,'visible','off');
	set(findall(ax, 'type', 'text'), 'visible', 'on');
end

% maximize window
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);