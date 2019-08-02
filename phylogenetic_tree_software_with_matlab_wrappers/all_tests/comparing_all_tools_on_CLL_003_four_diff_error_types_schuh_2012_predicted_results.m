%% comparing all tools on Schuh et al predicted output for CLL 003 data

% these are the results that Schuh et al give for the CLL data sets
CLL3_T = [
0 1 1 0 0
0 0 0 0 0
0 0 0 1 0
0 0 0 0 0
1 0 0 0 0];

CLL3_T_clust = [
1 1 
2 2
3 2
4 3
5 2
6 2
7 4
8 1
9 4 
10 1
11 1
12 2
13 3
14 4
15 4
16 2
17 3
18 1
19 1
20 2
];

U1 = CLL3_T;
clust1 = CLL3_T_clust;

all_error_of_ancestry_relations = [ ];

file_count_id = 23;

%getting EXACT data for real file no. 23
ourcode_output = load('/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_on_real_data_tree_sizes_6_9_top_20_16-Apr-2019-20-52-36.mat','all_our_code_outputs');
ourcode_output_check = ourcode_output.all_our_code_outputs{file_count_id};

U2 = inv(eye(size(ourcode_output_check{3})) - ourcode_output_check{3});
clust2 = [ [1: length(ourcode_output_check{6})]' , ourcode_output_check{6}];

%Computing error rates of EXACT run on file 23 vs CLL_003 Schuh et.al tree
[exact_error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);

%getting AncesTree data for real file no. 23
ancestree_output = load('/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/ancestree_all_files_minus_a_few_ELKEBIR_real_data_17-Apr-2019_01_21_08.mat','all_ancestree_output');
ancestree_output = ancestree_output.all_ancestree_output{file_count_id};

U2 = ancestree_output{3}{1};
U2 = inv(eye(length(U2)) - U2);
clust2 = [];
for i = 1:length(ancestree_output{4}{1})
	for j = ancestree_output{4}{1}{i}'
		clust2 = [clust2; [j,i]];
	end
end

%Computing error rates of AncesTree run on file 23 vs CLL_003 Schuh et.al tree
[ancestree_error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);

%getting PhyloWGS data for real file no. 23
phyloWGS_output = load('/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/phyloWGS_all_files_ELKEBIR_real_data_16-April-2019 18_27_30.mat','all_phylosub_outputs');
phyloWGS_output = phyloWGS_output.all_phylosub_outputs{file_count_id};

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

%Computing error rates of PhyloWGS run on file 23 vs CLL_003 Schuh et.al tree
[phylowgs_error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);

%getting CITUP data for real file no. 23
citup_output = load('/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/citup_all_files_ELKEBIR_real_data_17-Apr-2019_01_07_16.mat','all_citup_outputs');
citup_output = citup_output.all_citup_outputs{file_count_id};

U2 = citup_output{1}{1}{1};
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

clust2 = citup_output{1}{2};
clust2(:,2) = 1 + clust2(:,2);

%Computing error rates of CITUP run on file 23 vs CLL_003 Schuh et.al tree
[citup_error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);

%getting CANOPY data for real file no. 23
canopy_output = load('/home/surjray/Dropbox/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/all_canopy_outputs_on_real_Elkebir_data_24-Apr-2019.mat','all_canopy_outputs');
canopy_output = canopy_output.all_canopy_outputs{file_count_id};

[U2, clust2] = extract_U_mat_and_clust_from_canopy_output(canopy_output);

%Computing error rates of Canopy run on file 23 vs CLL_003 Schuh et.al tree
[canopy_error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);