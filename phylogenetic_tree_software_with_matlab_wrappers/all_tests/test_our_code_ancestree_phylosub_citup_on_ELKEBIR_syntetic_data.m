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

%% add paths for different tools

addpath(genpath('../EXACT/'));
addpath(genpath('../AncesTree/'));
addpath(genpath('../citup/'));
addpath(genpath('../phylowgs/'));
addpath(genpath('../Canopy/'));

%% process all 90 simulated files, provided in the AncesTree paper
% More information for the data is found at: 
% https://github.com/raphael-group/AncesTree/tree/master/data/simulated

pwd_start = pwd;

start_directory = [pwd_start, '/Sample_test_data/AncesTree_data/simulated'];
all_folders = dir(start_directory);

all_ancestree_errors = {};
all_citup_errors = {};
all_phylosub_errors = {};
all_our_code_errors = {};

all_ancestree_outputs = {};
all_citup_outputs = {};
all_phylosub_outputs = {};
all_our_code_outputs = {};

count_files_tested = 0;

% This option helps choose which tool we wish to run
% 1 = AncesTree
% 2 = CITUP
% 3 = PhyloWGS
% 4 = EXACT
% 5 = Canopy
tools_to_test = [4];

workspace_file_name = datestr(datetime);

for folder_ix = 3: size(all_folders,1)
    folder_name = all_folders(folder_ix).name;
    all_input_files_in_folder = dir([start_directory,'/',folder_name,'/*.input']);
    all_truth_files_in_folder = dir([start_directory,'/',folder_name,'/*.true']);

    for file_ix = 1: size(all_input_files_in_folder,1)
        
        plot(count_files_tested);
		drawnow;
        count_files_tested = count_files_tested + 1;
        
		% This if statement helps with the situation where some tools may
		% not run on certain files.
        if (count_files_tested >= 1)
        input_file_name = all_input_files_in_folder(file_ix).name;
        truth_file_name = all_truth_files_in_folder(file_ix).name;

        full_input_file_name = [start_directory,'/',folder_name,'/',input_file_name];
        full_truth_file_name = [start_directory,'/',folder_name,'/',truth_file_name];
        disp(full_input_file_name);
        disp(full_truth_file_name);

        
        %% read input data from file
        input_file = full_input_file_name;
        [F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);

        n_mutations = size(F_from_SampleData,1);
        T_samples = size(F_from_SampleData,2);

        %% read ground truth from file
        ground_truth_input_file = full_truth_file_name;
        [true_tree_data] =  read_ground_truth_from_elkebir_data(ground_truth_input_file);

        U_virtual_ground_truth= true_tree_data{3}';

        n_real_clus = size(U_virtual_ground_truth,1);

        U_real_ground_truth = zeros(n_mutations);

        for i = 1:n_real_clus
            list_of_real_nodes_corrs_to_currt_clust_node = true_tree_data{5}(true_tree_data{5}(:,2) == i,1);
            list_of_clust_ances_corrs_to_currt_clust_node = find(U_virtual_ground_truth(:,i));
            list_of_real_nodes_corrs_to_ancestors_of_currt_clust_node = [];
            for j = list_of_clust_ances_corrs_to_currt_clust_node'
                tmp  =  true_tree_data{5}(true_tree_data{5}(:,2) == j,1);
                list_of_real_nodes_corrs_to_ancestors_of_currt_clust_node = [list_of_real_nodes_corrs_to_ancestors_of_currt_clust_node, tmp'];
            end

            U_real_ground_truth(list_of_real_nodes_corrs_to_ancestors_of_currt_clust_node,list_of_real_nodes_corrs_to_currt_clust_node) = 1;
        end

        ground_truth_Tree_matrix_on_clustered_nodes = eye(n_real_clus)  - inv(U_virtual_ground_truth);
        ground_truth_clusters = true_tree_data{5};
        
        %% run ancestree on data

        scale = scaling;
        alpha = 0.3; % if alpha is big, lots of things will be clustered together
        beta = 0.8; % to choose a larger beta we need more samples , i.e. larger T_samples
        gamma = 0.01; % small gamma means larger confidence on the data
        wrapper_working_directory = [pwd_start, '/../AncesTree/distribution/'];
        ancestree_executable_path = [pwd_start, '/../AncesTree/distribution/build/ancestree'];

        % call the ancestree tool using ancestree_wrapper
		% ancestree_output is a matlab cell object having 6 components.
		% The first four components are in turn a cell object indexed by sol_id, which lists different good solutions.
		% The number of solutions that AncesTree outputs is given by how many different sol_id indices there are in the output
		%	ancestree_output{1}{sol_id} = recovered (clean) frequencies of clustered mutations.
		%	Each row is associated to a different sample, and each column to a different cluster of mutations
		%	ancestree_output{2}{sol_id} = clustered frequencies of mutants
		%	ancestree_output{3}{sol_id} = adjacency matrix for the optimal tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
		%	ancestree_output{4}{sol_id} = cluster membership information for the clustering
		%	associated to M{2}{sol_id} in the form of a cell array. The ith cell is
		%	an array that lists the nodes that belong to the ith cluster. These clusters are a subset of the clusters in M{6}
		%
		% The last two components are not indexed by sol_id, and are the same for all of the solutions that AncesTree outputs
		%	ancestree_output{5} = F_reduced'/2
		%	ancestree_output{6} = pre-clustering of mutations. An array with 2 columns, the 1st column designating the cluster ID, and the 2nd column designating the mutation that belongs to that cluster

        if (  ismember(1, tools_to_test) )
            ancestree_output = ancestree_wrapper(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
            all_ancestree_outputs{count_files_tested} = ancestree_output;
        else
            all_ancestree_outputs{count_files_tested} = nan;
		end
		
        %% run CITUP on data

        wrapper_working_directory = [pwd_start,'/../citup/distribution/'];
        CITUP_executable_path = [pwd_start,'/../citup/distribution/bin/'];
        min_cluster_no = 5;
        max_cluster_no = 9;
        citup_error_rate = 0.03;

        % call the CITUP tool using CITUP_wrapper_diff_clust_size
		% citup_output is a matlab cell object having 5 components, containg the best tree found using a
		% quadratic fitting cost and a Bayesian Information Criterion to select the
		% output tree size
		%	citup_output{1}{1} = adjacency matrix for the best output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
		%	citup_output{2} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
		%	citup_output{3}{1} = clustered frequencies of mutants. Rows are associated with samples, and columns with mutants
		%	citup_output{4}{1} = recovered (clean) frequencies of clustered mutations. Rows are associated with samples, and columns with mutations
		%	Note that citup_output{4}{1}' = U * citup_output{3}{1}'
		%	citup_output{5} = frequencies of clustered-mutations after the kmeans preclustering
		% all_inter_sol contains all the best tree found by CITUP for each tree size in the range min_clust_no and max_clust_no, indexed by tree_size_id here, each cell storing
		%	all_inter_sol{tree_size_id}{1}{1} = adjacency matrix for the optimal output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
		%	all_inter_sol{tree_size_id}{2} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
		%	all_inter_sol{tree_size_id}{3}{1} =  clustered frequencies of mutants. Rows are associated with samples, and columns with mutants
		%	all_inter_sol{tree_size_id}{4}{1} = recovered (clean) frequencies of clustered mutations. Rows are associated with samples, and columns with mutations
		%	all_inter_sol{tree_size_id}{5} =  frequencies of clustered-mutations after the kmeans preclustering
		
        if (ismember(2, tools_to_test))
            [~, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no,max_cluster_no, CITUP_executable_path, citup_error_rate )
            all_citup_outputs{count_files_tested} = {citup_output,all_inter_sol};
        else
            all_citup_outputs{count_files_tested} = nan;
        end
        
        %% here we run phyloWGS on the data
        scale = scaling;

        % Constants for PhyloWGS, should not be changed for sequencing data
        % one = 1;
        mu_r = 0.999;
        mu_v = 0.5;
        wrapper_working_directory = [pwd_start,'/../phylowgs/distribution/'];
        phylowgs_executable_path = [pwd_start,'/../phylowgs/distribution/'];
        
		% call the phyloWGS tool using PhyloWGS_wrapper
		% phylowgs_output is a matlab cell object having 3 components.
		% Each component is a cell object indexed by sol_id, which lists different good solutions.
		% The number of solutions that PhyloWGS outputs is given by how many different sol_id indices there are in the output
		% phylowgs_output{1}{sol_id} = adjacency matrix for the sol_id th output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
		% phylowgs_output{2}{sol_id} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
		% phylowgs_output{3}{sol_id} = frequencies of the input, after mutations have been clustered, helps in obtaining clustered frequencies of mutants.
		
        if (ismember(3, tools_to_test))
            [phylowgs_output] = PhyloWGS_wrapper(F_from_SampleData, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);
            all_phylosub_outputs{count_files_tested} = phylowgs_output;
        else
            all_phylosub_outputs{count_files_tested} = nan;
        end
        
        %% run our code on the data
        
        path_to_folder = [pwd_start, '/../EXACT/distribution/'];
        exec_name = 'EXACT_executable_x64_CUDA.out';
        cost_function = 'cost1';
        cpu_gpu = 'cpu_multithread';
        GPU_id = 0;
        min_tree_size = 7;
        max_tree_size = 10;
        num_CPU_cores = 50;
        error_rate = 0.03;
        max_num_partitions = 1;
        device_tree_subset_value = 1;
        CUDA_threads_per_block = 32;
        CUDA_blocks = 128;
        top_k_value = 20;
        
        
        if (ismember(4, tools_to_test))
            [ourcode_output] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);
            all_our_code_outputs{count_files_tested} = ourcode_output;
        else
            all_our_code_outputs{count_files_tested} = nan;
		end
        
		%% here we run Canopy on the simulated data
		%File name for the Ancestree format file to transform and input to Canopy		
		%Here input_file will serve as the variable containing the full
		%file name for the Ancestree input file
		
		path_to_folder = [pwd_start, '/../Canopy/demo_code/'];
		
		burnin_val = 10;
		thin_val = 5;
		K_min_val = 3;
		K_max_val = 6;
		numchains_val = 15;
		maxsimrun_val = 100000;
		minsimrun_val = 10000;
		writeskip_val = 200;
		
		cluster_number_start = 2;
		cluster_number_end = 7; % we use 9 for all files except 90
		
		%Calling the canopy_input_file_maker_function to transform input and run
		%Canopy inside a Rscript call
		if (ismember(5, tools_to_test))
            [canopy_output] = canopy_wrapper(input_file, path_to_folder, burnin_val, thin_val, K_min_val, K_max_val, numchains_val, maxsimrun_val, minsimrun_val, writeskip_val, cluster_number_start, cluster_number_end);
			all_canopy_outputs{count_files_tested} = canopy_output;
		else
			all_canopy_outputs{count_files_tested} = nan;
        end
		
        %% clear stuff on the command line and history to make sure there are no funny out of memory errors from java
        save(workspace_file_name);
        clc;
        com.mathworks.mlservices.MLCommandHistoryServices.removeAll;
		end
        
    end 
end