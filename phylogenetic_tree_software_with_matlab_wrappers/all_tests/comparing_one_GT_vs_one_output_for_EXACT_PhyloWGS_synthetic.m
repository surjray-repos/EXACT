% Copyright (c) 2019 Surjyendu Ray
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%% Comparing ground truth and output of EXACT across all 90 files
% Compares all U matrices and clusters

start_folder = pwd;

all_error_of_ancestry_relations_EXACT = [ ];
all_our_code_outputs_new = {};

figure;
for file_count_id = 1:90
    %% read input data from file

    path_to_sim_ances_tree_data = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/';
    list_of_dirs = dir(path_to_sim_ances_tree_data);
    count = 0;
    found_flag = 0;
    for i = 3:length(list_of_dirs)
        list_of_fildes = dir([path_to_sim_ances_tree_data,list_of_dirs(i).name,'/*.input']);
        list_of_fildes_truth = dir([path_to_sim_ances_tree_data,list_of_dirs(i).name,'/*.true']);

        for j = 1:length(list_of_fildes)
            count = count + 1;
            if (count == file_count_id)
                full_input_file_name = [path_to_sim_ances_tree_data,list_of_dirs(i).name,'/',list_of_fildes(j).name];
                ground_truth_file = [path_to_sim_ances_tree_data,list_of_dirs(i).name,'/',list_of_fildes_truth(j).name];
                found_flag = 1;
                break;
            end
            if (found_flag == 1)
                break;
            end
        end
    end

    [F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(full_input_file_name);

    n_mutations = size(F_from_SampleData,1);
    T_samples = size(F_from_SampleData,2);

    %% read ground truth

    [true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
    U1 = true_tree_data{3}';
    clust1 = true_tree_data{5};
    
    %% get correct form of matrices before comparing
    ourcode_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_all_files_synthetic_data_top_20_size_7_10.mat','all_our_code_outputs_top_k');
	%ourcode_output_check = ourcode_output.all_our_code_outputs_new{file_count_id};
    ourcode_output_check = ourcode_output.all_our_code_outputs_top_k{file_count_id};
	
    U2 = inv(eye(size(ourcode_output_check{3})) - ourcode_output_check{3} );
    clust2 = zeros(100,2);
    clust2(:,1) = 1:100;
    clust2(:,2) = ourcode_output_check{6};

    %%
    [error_rates_EXACT] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);   
    all_error_of_ancestry_relations_EXACT = [all_error_of_ancestry_relations_EXACT , error_rates_EXACT(2)];
    plot(all_error_of_ancestry_relations_EXACT,'*');
    drawnow; 
end

%% Comparing ground truth and output of EXACT across all 90 files
% Compares all U matrices and clusters and storing the errors for each tree
% size separately, with the top k trees also

start_folder = pwd;

all_error_of_ancestry_relations_EXACT_1 = [ ];
all_error_of_ancestry_relations_EXACT_2 = [ ];
all_error_of_ancestry_relations_EXACT_3 = [ ];
all_error_of_ancestry_relations_EXACT_4 = [ ];
all_our_code_outputs_new = {};

figure;
for file_count_id = 1:90
    %% read input data from file

    path_to_sim_ances_tree_data = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/';
    list_of_dirs = dir(path_to_sim_ances_tree_data);
    count = 0;
    found_flag = 0;
    for i = 3:length(list_of_dirs)
        list_of_fildes = dir([path_to_sim_ances_tree_data,list_of_dirs(i).name,'/*.input']);
        list_of_fildes_truth = dir([path_to_sim_ances_tree_data,list_of_dirs(i).name,'/*.true']);

        for j = 1:length(list_of_fildes)
            count = count + 1;
            if (count == file_count_id)
                full_input_file_name = [path_to_sim_ances_tree_data,list_of_dirs(i).name,'/',list_of_fildes(j).name];
                ground_truth_file = [path_to_sim_ances_tree_data,list_of_dirs(i).name,'/',list_of_fildes_truth(j).name];
                found_flag = 1;
                break;
            end
            if (found_flag == 1)
                break;
            end
        end
    end

    [F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(full_input_file_name);

    n_mutations = size(F_from_SampleData,1);
    T_samples = size(F_from_SampleData,2);

    %% read ground truth

    [true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
    U1 = true_tree_data{3}';
    clust1 = true_tree_data{5};
    
    %% get correct form of matrices before comparing
    ourcode_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_all_files_synthetic_data_top_20_size_7_10.mat','all_our_code_outputs_top_k');
	%ourcode_output_check = ourcode_output.all_our_code_outputs_new{file_count_id};
    ourcode_output_check = ourcode_output.all_our_code_outputs_top_k{file_count_id};
	
    U2 = inv(eye(size(ourcode_output_check{3})) - ourcode_output_check{3} );
    clust2 = zeros(100,2);
    clust2(:,1) = 1:100;
    clust2(:,2) = ourcode_output_check{6};

    %%
    [error_rates_EXACT] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);   
	all_error_of_ancestry_relations_EXACT_1 = [all_error_of_ancestry_relations_EXACT_1 , error_rates_EXACT(1)];
    all_error_of_ancestry_relations_EXACT_2 = [all_error_of_ancestry_relations_EXACT_2 , error_rates_EXACT(2)];
	all_error_of_ancestry_relations_EXACT_3 = [all_error_of_ancestry_relations_EXACT_3 , error_rates_EXACT(3)];
	all_error_of_ancestry_relations_EXACT_4 = [all_error_of_ancestry_relations_EXACT_4 , error_rates_EXACT(4)];
    %plot(all_error_of_ancestry_relations_EXACT,'*');
    %drawnow; 
end

%% Comparing ground truth and output of phyloWGS across all 90 files
% Compares all U matrices and clusters
all_error_of_ancestry_relations_phyloWGS_1 = [ ];
all_error_of_ancestry_relations_phyloWGS_2 = [ ];
all_error_of_ancestry_relations_phyloWGS_3 = [ ];
all_error_of_ancestry_relations_phyloWGS_4 = [ ];

figure;
for file_count_id = 1:90
    
    if (file_count_id ~= 35)
    
        % this identifies the file with number equal to file_count_id
        path_to_sim_ances_tree_data = '/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/Sample_test_data/AncesTree_data/simulated/';
        list_of_dirs = dir(path_to_sim_ances_tree_data);
        count = 0;
        found_flag = 0;
        for i = 3:length(list_of_dirs)
            list_of_fildes = dir([path_to_sim_ances_tree_data,list_of_dirs(i).name,'/*.true']);
            for j = 1:length(list_of_fildes)
                count = count + 1;
                if (count == file_count_id)
                    ground_truth_file = [path_to_sim_ances_tree_data,list_of_dirs(i).name,'/',list_of_fildes(j).name];
                    found_flag = 1;
                    break;
                end
                if (found_flag == 1)
                    break;
                end
            end
        end

        [true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_file);
        phyloWGS_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/phyloWGS_all_files_except_35_ELKEBIR_simulated_data.mat','all_phylosub_outputs');
        phyloWGS_output = phyloWGS_output.all_phylosub_outputs{file_count_id};

        U1 = true_tree_data{3}';
        clust1 = true_tree_data{5};

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
        clust2 = zeros(100,2);
        for i = 1:size(phyloWGS_output{2}{1},1)
            clust2(1+phyloWGS_output{2}{1}{i,1},:) = [1+phyloWGS_output{2}{1}{i,1},1+phyloWGS_output{2}{1}{i,2}];
        end

        [error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);
		
		all_error_of_ancestry_relations_phyloWGS_1 = [all_error_of_ancestry_relations_phyloWGS_1 , error_rates(1)];
        all_error_of_ancestry_relations_phyloWGS_2 = [all_error_of_ancestry_relations_phyloWGS_2 , error_rates(2)];
		all_error_of_ancestry_relations_phyloWGS_3 = [all_error_of_ancestry_relations_phyloWGS_3 , error_rates(3)];
		all_error_of_ancestry_relations_phyloWGS_4 = [all_error_of_ancestry_relations_phyloWGS_4 , error_rates(4)];
        %plot(all_error_of_ancestry_relations_phyloWGS, '*');
		%drawnow;
	end
end