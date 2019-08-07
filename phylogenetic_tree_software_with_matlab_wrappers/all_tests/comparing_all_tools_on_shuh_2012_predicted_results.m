load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_on_real_data_tree_sizes_6_9_top_20_16-Apr-2019-20-52-36.mat','all_our_code_outputs');
load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/ancestree_all_files_minus_a_few_ELKEBIR_real_data_17-Apr-2019_01_21_08.mat','all_ancestree_output');
load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/phyloWGS_all_files_ELKEBIR_real_data_16-April-2019 18_27_30.mat','all_phylosub_outputs');
load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/citup_all_files_ELKEBIR_real_data_17-Apr-2019_01_07_16.mat','all_citup_outputs');
load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/all_canopy_outputs_on_real_Elkebir_data_24-Apr-2019.mat','all_canopy_outputs');


%% comparing all tools on Schuh et al predicted outputs

% these are the results that Schuh et al give for the CLL data sets
CLL3_M = [3 0 0 0 0
82 91 13 0 0 
9 6 12 5 3
9 1 44 89 96
6 2 30 6 1];
M_CLL3_M_normalized = CLL3_M*diag(1./sum(CLL3_M));

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

CLL77_M =[20 15 16 17 8
39 27 30 26 13
33 52 50 39 10
20 3 4 14 31
8 3 0 4 38];
M_CLL77_M_normalized = CLL77_M*diag(1./sum(CLL77_M));


CLL77_T = [0 1 0 1 0
0 0 1 0 0
0 0 0 0 0
0 0 0 0 0
1 0 0 0 0];

CLL6_M = [32 9 30 6 9
31 31 32 22 17
11 8 11 15 4
24 43 24 44 57
0 0 0 0 0
2 9 3 13 13];
M_CLL6_M_normalized = CLL6_M*diag(1./sum(CLL6_M));


CLL6_T = [0 1 0 0 0 0
0 0 1 0 0 0
0 0 0 1 0 0
0 0 0 0 0 1
1 0 0 0 0 0
0 0 0 0 0 0];

%% comparison with CLL3 

CLLfile = 2;
switch CLLfile
	case 1
		M_CLL_normalized = M_CLL3_M_normalized;
		file_ix = 23;
	case 2
		M_CLL_normalized = M_CLL6_M_normalized;
		file_ix = 25;
	case 3
		M_CLL_normalized = M_CLL77_M_normalized;
		file_ix = 27;
end

alg = 5;
switch alg
	case 1
		% Comparing with EXACT
		M_target = all_our_code_outputs{file_ix}{4};
	case 2
		sol_ix = 1;
		% Comparing with PhyloWGS
		U2 = all_phylosub_outputs{file_ix}{1}{sol_ix};
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
		
		M_target = inv(U2)*all_phylosub_outputs{file_ix}{3}{sol_ix};
	case 3
		% Comparing with AncesTree
		output_1  = all_ancestree_output{file_ix};
		sol_ix = 1;
        U2 = output_1{3}{sol_ix};
        U2 = inv(eye(length(U2)) - U2);
		
		M_target = output_1{2}{sol_ix}';
	case 4
		% Comparing with Canopy
		output_1  = all_canopy_outputs{file_ix};
		%[U2, clust2]= extract_U_mat_and_clust_from_canopy_output(output_2);
		M_target = output_1{2};
	case 5
		%Comparing with CITUP
		output_1  = all_citup_outputs{file_ix};
		sol_ix = 1;
		U2 = output_1{sol_ix}{1}{1};
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
		
		M_target = output_1{sol_ix}{3}{1}';
end


num_muts = size(M_target,1);
num_muts_other = size(M_CLL_normalized,1);

num_sampl = size(M_target,2);
all_perms = perms([1:num_muts]);

if (num_muts > num_muts_other)
	M_CLL_normalized_extended = M_CLL_normalized;
	M_CLL_normalized_extended(end+1:end+num_muts-num_muts_other,:) = 0;
else
	M_CLL_normalized_extended = M_CLL_normalized;
	M_target(end+1:end+num_muts_other-num_muts,:) = 0;
	num_muts = num_muts_other;
	all_perms = perms([1:num_muts]);
end

small_dist = inf;
for i = 1:size(all_perms,1)
	
	M_target_shuffled = M_target(all_perms(i,:),:);
	
	dist = sum(sum(abs(M_CLL_normalized_extended - M_target_shuffled)))/(num_muts*num_sampl);
	if (dist < small_dist)
		small_dist = dist;
		best_ix = i;
	end
end
disp(small_dist);

