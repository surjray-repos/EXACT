%% running all tools on Ralph's simple data set


%% load the data set

%input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/ralphs_simple_data_set_all_mutations_except_mut_5';

input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/ralphs_second_simple_data_set_all_mutations';

%input_file = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/ralphs_third_simple_data_set_all_mutations';


[F, scale] = transform_elkebir_input_data_into_F_matrix(input_file);

figure;
%% run our code
base_folder = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/EXACT';

run([base_folder, '/distribution/cvx/cvx/cvx_setup']);


% parameters
path_to_folder = [base_folder, '/distribution/'];
exec_name = 'b.out';
cost_function = 'cost1'; %possible options: cost1, cost2, cost3, cost4
cpu_gpu = 'cpu_multithread'; %possible options: cpu, cpu_multithread, gpu
GPU_id = 0;
min_tree_size = round(size(F,1)*2/3);
max_tree_size = size(F,1);
num_CPU_cores = 50;
error_rate = 0.03;
max_num_partitions = 1;
device_tree_subset_value = 1;
CUDA_threads_per_block = 32;
CUDA_blocks = 128;
top_k_value = min(100,(min_tree_size+1)^(min_tree_size-1));

% run EXACT
[ourcode_output, ourcode_best_bic, ourcode_all_Ms] = EXACT_wrapper_diff_tree_size(F, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);

sorted_ids = nan(top_k_value,1+max_tree_size-min_tree_size);
sorted_vals = nan(top_k_value,1+max_tree_size-min_tree_size);

for j = 1:1+max_tree_size-min_tree_size
	for i = 1:top_k_value
		sorted_vals(i,j) = ourcode_all_Ms{j}{i}{1};
	end
	[sorted_vals(:,j), sorted_ids(:,j)] = sort(sorted_vals(:,j));
end

best_ourcode_output = ourcode_output;

M_target = best_ourcode_output{4};


U2 = inv(eye(size(best_ourcode_output{3})) - best_ourcode_output{3});
clust2 = [ [1: length(best_ourcode_output{6})]' , best_ourcode_output{6}];

nodelbs = cell(1,length(U2));
for i = 1:length(U2)
	nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
end

subplot(2,5,1);
h = plot(digraph(  eye(length(U2)) - inv(U2)      ));

labelnode(h,[1:length(U2)],nodelbs);

cum_nodelbs = cell(1,length(U2));
node_col = cell(1,length(U2));
for i = 1:length(U2)
	for j = find(U2(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
	clrix = randi(64);
	node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end

subplot(2,5,6);
%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col);
%donut(M_target(:,1:5)',cum_nodelbs,node_col);
M_target = M_target(:,[1,4,2,5,3,6]);
generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
set(gca,'visible','off')

%% run PhyloWGS
base_folder = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/phylowgs';

%Constants for PhyloWGS, should not be changed for sequencing data
mu_r = 0.999;
mu_v = 0.5;
wrapper_working_directory = [base_folder,'/distribution/'];
phylowgs_executable_path = [base_folder,'/distribution/'];

% run PhyloWGS
[phyloWGS_output] = PhyloWGS_wrapper(F, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);

U2 = phyloWGS_output{1}{1};
root_phyloWGS = 1; % the root of the treee is always node 1.
n_virt_nodes_phyloWGS = size(U2,1);
AdjLrecon = {};
for j =1:n_virt_nodes_phyloWGS
	AdjLrecon{j} = find(U2(:,j));
end
Treerecon = BFS(AdjLrecon,root_phyloWGS); %root is the node 0 in unlabled tree
% get adj mat of tree
AdjTrecon = zeros(n_virt_nodes_phyloWGS);
for j =1:n_virt_nodes_phyloWGS
	for r = Treerecon{j}
		AdjTrecon(j,r) = 1;
	end
end
U2 = inv(eye(n_virt_nodes_phyloWGS) - AdjTrecon);
num_mutants = size(phyloWGS_output{2}{1},1);
clust2 = zeros(num_mutants,2);
for i = 1:size(phyloWGS_output{2}{1},1)
	clust2(1 + phyloWGS_output{2}{1}{i,1},:) = [1 + phyloWGS_output{2}{1}{i,1},1 + phyloWGS_output{2}{1}{i,2}];
end

M_target = inv(U2)*phyloWGS_output{3}{1};


nodelbs = cell(1,length(U2));
for i = 1:length(U2)
	nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
end

subplot(2,5,2);
h = plot(digraph(  eye(length(U2)) - inv(U2)      ));

labelnode(h,[1:length(U2)],nodelbs);

cum_nodelbs = cell(1,length(U2));
node_col = cell(1,length(U2));
for i = 1:length(U2)
	for j = find(U2(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
	clrix = randi(64);
	node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end

subplot(2,5,7);
%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col);
%donut(M_target(:,1:5)',cum_nodelbs,node_col);
M_target = M_target(:,[1,4,2,5,3,6]);
generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
set(gca,'visible','off')
%% run AncesTree

base_folder = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/AncesTree';

alpha = 0.3; % if alpha is big, lots of things will be clustered together
beta = 0.8; % to choose a larger beta we need more samples , i.e. larger T_samples
gamma = 0.01; % small gamma means larger confidence on the data
wrapper_working_directory = [base_folder, '/distribution/'];
ancestree_executable_path = [base_folder, '/distribution/build/ancestree'];

% run AncesTree
ancestree_output = ancestree_wrapper(F, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);

M_target = ancestree_output{2}{1}';


U2 = ancestree_output{3}{1};
U2 = inv(eye(length(U2)) - U2);
clust2 = [];
for i = 1:length(ancestree_output{4}{1})
	for j = ancestree_output{4}{1}{i}'
		clust2 = [clust2; [j,i]];
	end
end

nodelbs = cell(1,length(U2));
for i = 1:length(U2)
	nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
end

subplot(2,5,3);
h = plot(digraph(  eye(length(U2)) - inv(U2)      ));

labelnode(h,[1:length(U2)],nodelbs);

cum_nodelbs = cell(1,length(U2));
node_col = cell(1,length(U2));
for i = 1:length(U2)
	for j = find(U2(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
	clrix = randi(64);
	node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end

subplot(2,5,8);
%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col);
%donut(M_target(:,1:5)',cum_nodelbs,node_col);
M_target = M_target(:,[1,4,2,5,3,6]);

generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
set(gca,'visible','off')
%% run Canopy

base_folder =  '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/Canopy';

path_to_folder = [base_folder, '/demo_code/'];

% parameters
burnin_val = 10;
thin_val = 5;
K_min_val = 3;
K_max_val = 3;
numchains_val = 15;
maxsimrun_val = 100000;
minsimrun_val = 10000;
writeskip_val = 200;
cluster_number_start = round(size(F,1)*2/3);
cluster_number_end = size(F,1)-1;

% why doesn't canopy output multiple tree sizes? 
canopy_output = canopy_wrapper(input_file, path_to_folder, burnin_val, thin_val, K_min_val, K_max_val, numchains_val, maxsimrun_val, minsimrun_val, writeskip_val, cluster_number_start, cluster_number_end);

[U2 , clust2] = extract_U_mat_and_clust_from_canopy_output(canopy_output);

M_target = canopy_output{2};


nodelbs = cell(1,length(U2));
for i = 1:length(U2)
	nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
end

subplot(2,5,4);
h = plot(digraph(  eye(length(U2)) - inv(U2)      ));

labelnode(h,[1:length(U2)],nodelbs);

cum_nodelbs = cell(1,length(U2));
node_col = cell(1,length(U2));
for i = 1:length(U2)
	for j = find(U2(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
	clrix = randi(64);
	node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end

subplot(2,5,9);
%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col);
%donut(M_target(:,1:5)',cum_nodelbs,node_col);
M_target = M_target(:,[1,4,2,5,3,6]);

generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
set(gca,'visible','off')
%% run CITUP
base_folder = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/citup';

% parameters
wrapper_working_directory = [base_folder,'/distribution/'];
CITUP_executable_path = [base_folder,'/distribution/bin/'];
min_cluster_no = round(size(F,1)*2/3); % here we make the number of clusters be equal to the size of the input tree
max_cluster_no = size(F,1);
citup_error_rate = 0.03;

% run Citup
[sol_min_bic, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate );

U2 = citup_output{1}{1};
root_citup = 1; % the root of the treee is always node 1.
n_virt_nodes_citup = size(U2,1);
AdjLrecon = {};
for j =1:n_virt_nodes_citup
	AdjLrecon{j} = find(U2(:,j));
end
Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
% get adj mat of tree
AdjTrecon = zeros(n_virt_nodes_citup);
for j =1:n_virt_nodes_citup
	for r = Treerecon{j}
		AdjTrecon(j,r) = 1;
	end
end
U2 = inv(eye(n_virt_nodes_citup) - AdjTrecon);
clust2 = citup_output{2};
clust2(:,2) = 1 + clust2(:,2);

M_target = citup_output{3}{1}';


nodelbs = cell(1,length(U2));
for i = 1:length(U2)
	nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
end

subplot(2,5,5);
h = plot(digraph(  eye(length(U2)) - inv(U2)      ));

labelnode(h,[1:length(U2)],nodelbs);

cum_nodelbs = cell(1,length(U2));
node_col = cell(1,length(U2));
for i = 1:length(U2)
	for j = find(U2(:,i))'
		cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
	end
	array_rep = [0, sort(str2num(cum_nodelbs{i}))];
	colormap_array = colormap(colorcube);
	rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
	clrix = randi(64);
	node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
	cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
end

subplot(2,5,10);
%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col);
%donut(M_target(:,1:5)',cum_nodelbs,node_col);
M_target = M_target(:,[1,4,2,5,3,6]);

generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
set(gca,'visible','off')

%% we can also run EXACT to generate not one, but multiple trees with the corresponding probabilities

figure;
tree_probs = [];
for t = 1:10 % this will show us the top 10 trees
	best_size = size(ourcode_output{3},1); % we can choose which size to work with here
	best_tree_ix = sorted_ids( t   ,best_size-min_tree_size);
	best_ourcode_output = ourcode_all_Ms{best_size-min_tree_size}{best_tree_ix};

	tree_probs(t) = (exp(-best_ourcode_output{1}./(2*0.03*0.03)))/sum(exp(-sorted_vals(:,best_size-min_tree_size)./(2*0.03*0.03)));

	M_target = best_ourcode_output{4};

	U2 = inv(eye(size(best_ourcode_output{3})) - best_ourcode_output{3});
	clust2 = [ [1: length(best_ourcode_output{6})]' , best_ourcode_output{6}];

	nodelbs = cell(1,length(U2));
	for i = 1:length(U2)
		nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
	end

	ax=subplot(2,10,t);
	h = plot(digraph(  eye(length(U2)) - inv(U2)      ));
	labelnode(h,[1:length(U2)],nodelbs);
	title(ax,num2str(tree_probs(t)));
	
	cum_nodelbs = cell(1,length(U2));
	node_col = cell(1,length(U2));
	for i = 1:length(U2)
		for j = find(U2(:,i))'
			cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
		end
		array_rep = [0, sort(str2num(cum_nodelbs{i}))];
		colormap_array = colormap(colorcube);
		rng(round((sum(((7.534).^[1:length(array_rep)]).*array_rep))));
		clrix = randi(64);
		node_col{i} = colormap_array(clrix   ,:); %[1/(1+var(array_rep)/(1+mean(array_rep))), (1+median(array_rep))/(1+norm(array_rep)),min(1,(1+length(array_rep))/(1+max(array_rep)))]; 
		cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
	end

	subplot(2,10,t+10);
	%donut((M_target(:,1:3)+M_target(:,4:6))',cum_nodelbs,node_col)
	%donut(M_target(:,1:5)',cum_nodelbs,node_col)
	M_target = M_target(:,[1,4,2,5,3,6]);

	generate_simple_muller_plots(U2,M_target,cum_nodelbs,node_col);
	set(gca,'visible','off')
end