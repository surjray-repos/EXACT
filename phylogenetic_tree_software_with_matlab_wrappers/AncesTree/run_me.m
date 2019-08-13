%% example of how to use Ancestree

% input file
input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_1.input';
[F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

% load groud truth
ground_truth_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_1.true';
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
wrapper_working_directory = [pwd, '/distribution/'];
ancestree_executable_path = [pwd, '/distribution/build/ancestree'];

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

% generate figure with all the output trees inferred by AncesTree
%%
figure;
num_solutions = size(ancestree_output{1}, 1);

for sol_id = 1:num_solutions % this will show us all output trees
	Mutant_Frequencies_M = ancestree_output{2}{sol_id}';
	
	U = ancestree_output{3}{sol_id};
	U = inv(eye(length(U)) - U);
	clust = [];
	for i = 1:length(ancestree_output{4}{sol_id})
		for j = ancestree_output{4}{sol_id}{i}'
			clust = [clust; [j,i]];
		end
	end
	
	cost = compute_ancestree_fitting_cost(F_from_SampleData, clust,ancestree_output{1}{sol_id});
	
	nodelbs = cell(1,length(U));
	for i = 1:length(U)
		nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
	end
	
	ax = subplot(2, num_solutions/2, sol_id);
	h = plot(digraph(  eye(length(U)) - inv(U)  ));
	labelnode(h, [1:length(U)],nodelbs);
	if (sol_id == 1)
		title(ax, ['AncesTree fitting cost of different trees: ' , num2str(cost)]);
	else
		title(ax, [num2str(cost)]);
	end
	set(ax,'visible','off');
	set(findall(ax, 'type', 'text'), 'visible', 'on');
end
% maximize window
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
