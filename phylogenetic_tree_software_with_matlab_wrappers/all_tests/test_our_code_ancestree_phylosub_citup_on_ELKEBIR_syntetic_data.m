%% add paths for different tools

addpath(genpath('../EXACT/'));
addpath(genpath('../AncesTree-master/'));
addpath(genpath('../citup-master'));
addpath(genpath('../phylowgs-master'));

%% process all the files

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

tools_to_test = [2];

workspace_file_name = datestr(datetime);

for folder_ix = 3: size(all_folders,1)
    folder_name = all_folders(folder_ix).name;
    all_input_files_in_folder = dir([start_directory,'/',folder_name,'/*.input']);
    all_truth_files_in_folder = dir([start_directory,'/',folder_name,'/*.true']);

    for file_ix = 1: size(all_input_files_in_folder,1)
        
        plot(count_files_tested);
		drawnow;
        count_files_tested = count_files_tested + 1;
        
        if (count_files_tested == 85)
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
        wrapper_working_directory = [pwd_start, '/../AncesTree-master/distribution/'];
        ancestree_executable_path = [pwd_start, '/../AncesTree-master/distribution/build/ancestree'];

        % call ancestree
        %ancestree_output = ancestree_wrapper_ll(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
        %ancestree_output = ancestree_wrapper_ll_noTransform(input_file,  alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
        
        if (  ismember(1,tools_to_test) )
            %ancestree_output = ancestree_wrapper_reading_only_sol_file(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
            ancestree_output = ancestree_wrapper(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
            all_ancestree_outputs{count_files_tested} = ancestree_output;
        else
            all_ancestree_outputs{count_files_tested} = nan;
        end
%         num_diff_sols = length(ancestree_output{1});
%         % go over all solutions and compute the error
%         ances_multi_sol_per_file = {};
%         for sol_id = 1:num_diff_sols
%              % get U from clustered data
%             clustered_ancestree = ancestree_output{3}{sol_id};
%             root_of_tree = find(sum(clustered_ancestree)==0);
% 
%             clustered_ancestree = clustered_ancestree + clustered_ancestree';
%             n_ances_clus = size(clustered_ancestree,1);
%             AdjLrecon = {};
%             for j =1:n_ances_clus
%                 AdjLrecon{j} = find(clustered_ancestree(:,j));
%             end
% 
%             Treerecon = BFS(AdjLrecon,root_of_tree);
%             % get adj mat of tree
%             AdjTrecon = zeros(n_ances_clus);
%             for j =1:n_ances_clus
%                 for r = Treerecon{j}
%                     AdjTrecon(j,r) = 1;
%                 end
%             end
% 
%             U_recon = inv(eye(n_ances_clus) - AdjTrecon);
%             
%             clusters_for_sol = ancestree_output{4}{sol_id};
%             all_relevant_muts = cell2mat(clusters_for_sol);
%             list_of_errors = [];
%             for ix = 1:length(all_relevant_muts)
%                 for jx = ix+1:length(all_relevant_muts)
%                     
%                         i = all_relevant_muts(ix);
%                         j = all_relevant_muts(jx);
% 
%                         i_gt_virt = ground_truth_clusters(i,2);
%                         j_gt_virt = ground_truth_clusters(j,2);
%                         
%                         ground_truth_flag = zeros(3,1); % 1 0 0 = clustered 0 1 0 = i_ances_j 0 -1 0 = j_ances_i 0 0 1 = incommp
%                         if (i_gt_virt==j_gt_virt)
%                             ground_truth_flag(1) = 1;
%                         else
%                             if ( ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 1)
%                                 ground_truth_flag(2) = 1;
%                             end
%                             if ( ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 1)
%                                 ground_truth_flag(2) = -1;
%                             end
%                             if(ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 0 && ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 0)
%                                 ground_truth_flag(3) = 1;
%                             end
% 
%                         end
%                     
%                         i_anc_virt = nan;
%                         j_anc_virt = nan;
%                         for r = 1:n_ances_clus
%                             if (ismember(i,clusters_for_sol{r}))
%                                 i_anc_virt = r;
%                             end
%                            if (ismember(j,clusters_for_sol{r}))
%                                 j_anc_virt = r;
%                             end
%                         end
%                         
%                         ances_flag = zeros(3,1);
%                         if (i_anc_virt==j_anc_virt)
%                             ances_flag(1) =  1;
%                         else
%                             if ( U_recon(i_anc_virt,j_anc_virt) == 1)
%                                 ances_flag(2) = 1;
%                             end
%                             if ( U_recon(j_anc_virt,i_anc_virt) == 1)
%                                 ances_flag(2) = -1;
%                             end
%                             if(U_recon(i_anc_virt,j_anc_virt) == 0 && U_recon(j_anc_virt,i_anc_virt) == 0)
%                                 ances_flag(3) = 1;
%                             end
%                         end
%                         list_of_errors =  [list_of_errors, [i;j;ances_flag ~= ground_truth_flag]];
%                         
%                 end
%             end
%             
%             error_rates = mean(list_of_errors(3:5,:),2);
%             ances_multi_sol_per_file{sol_id} = error_rates;
%         end
%         
%         all_ancestree_errors{count_files_tested} = ances_multi_sol_per_file;

      
    
        
        %% run CITUP on data

        wrapper_working_directory = [pwd_start,'/../citup-master/distribution/'];
        CITUP_executable_path = [pwd_start,'/../citup-master/distribution/bin/'];
        min_cluster_no = 5; % here we make the number of clusters be equal to the size of the input tree
        max_cluster_no = 9;
        citup_error_rate = 0.03;

        %call our CITUP wrapper
        %[citup_output] = count_general_wrapper_CITUP_bowie2_OC_ll(F_from_SampleData, wrapper_working_directory, cluster_no, CITUP_executable_path, citup_error_rate);

        if (ismember(2,tools_to_test))
            [~, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no,max_cluster_no, CITUP_executable_path, citup_error_rate )
            all_citup_outputs{count_files_tested} = {citup_output,all_inter_sol};
        else
            all_citup_outputs{count_files_tested} = nan;
        end
        
%         num_citup_sol = length(citup_output{1});
%         errors_CITUP_vs_groundtruth_on_U_matrix = nan(num_citup_sol,3);
% 
%         for sol_id = 1:num_citup_sol
%         
%             AdjT_citup = citup_output{1}{sol_id}; % note that citup already outputs an adjacancy matrix
%             root_citup = 1; % the root of the treee is always node 1.
%             n_citup_clus = size(AdjT_citup,1);
%             AdjLrecon = {};
%             for j =1:n_citup_clus
%                 AdjLrecon{j} = find(AdjT_citup(:,j));
%             end
% 
%             Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
%             % get adj mat of tree
%             AdjTrecon = zeros(n_citup_clus);
%             for j =1:n_citup_clus
%                 for r = Treerecon{j}
%                     AdjTrecon(j,r) = 1;
%                 end
%             end
% 
%             U_recon = inv(eye(n_citup_clus) - AdjTrecon);
%         
%             % compute the error accross the categories ancestral, clustered
%             % and incomparable accros different pairs of mutants i and j
%             all_relevant_muts = [1:n_mutations];
%             list_of_errors = [];
%             for ix = 1:length(all_relevant_muts)
%                 for jx = ix+1:length(all_relevant_muts)
%                     
%                         i = all_relevant_muts(ix);
%                         j = all_relevant_muts(jx);
% 
%                         i_gt_virt = ground_truth_clusters(i,2);
%                         j_gt_virt = ground_truth_clusters(j,2);
%                         
%                         ground_truth_flag = zeros(3,1); % 1 0 0 = clustered 0 1 0 = i_ances_j 0 -1 0 = j_ances_i 0 0 1 = incommp
%                         if (i_gt_virt==j_gt_virt)
%                             ground_truth_flag(1) = 1;
%                         else
%                             if ( ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 1)
%                                 ground_truth_flag(2) = 1;
%                             end
%                             if ( ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 1)
%                                 ground_truth_flag(2) = -1;
%                             end
%                             if(ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 0 && ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 0)
%                                 ground_truth_flag(3) = 1;
%                             end
% 
%                         end
%                     
%                         i_citup_virt = citup_output{6}{sol_id}(i);
%                         j_citup_virt = citup_output{6}{sol_id}(j);
%                         
%                         citup_flag = zeros(3,1);
%                         if (i_citup_virt==j_citup_virt)
%                             citup_flag(1) =  1;
%                         else
%                             if ( U_recon(i_citup_virt,j_citup_virt) == 1)
%                                 citup_flag(2) = 1;
%                             end
%                             if ( U_recon(j_citup_virt,i_citup_virt) == 1)
%                                 citup_flag(2) = -1;
%                             end
%                             if(U_recon(i_citup_virt,j_citup_virt) == 0 && U_recon(j_citup_virt,i_citup_virt) == 0)
%                                 citup_flag(3) = 1;
%                             end
%                         end
%                         list_of_errors =  [list_of_errors, [i;j;citup_flag ~= ground_truth_flag]];
%                         
%                 end
%             end
%             
%             error_rates = mean(list_of_errors(3:5,:),2);
%             errors_CITUP_vs_groundtruth_on_U_matrix(sol_id,:) = error_rates;
%             
%         end
%         all_citup_errors{count_files_tested} = errors_CITUP_vs_groundtruth_on_U_matrix;

     
        %% here we run phyloWGS on the data
        scale = scaling;

        %Constants for PhyloWGS, should not be changed for sequencing data
        %one = 1;
        mu_r = 0.999;
        mu_v = 0.5;
        wrapper_working_directory = [pwd_start,'/../phylowgs-master/distribution/'];
        phylowgs_executable_path = [pwd_start,'/../phylowgs-master/distribution/'];
        
        if (ismember(3,tools_to_test))
            [phylowgs_output] = PhyloWGS_wrapper(F_from_SampleData, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);
            all_phylosub_outputs{count_files_tested} = phylowgs_output;
        else
            all_phylosub_outputs{count_files_tested} = nan;
        end
        
        %% run our code on the data
        
        path_to_folder = [pwd_start, '/../EXACT/distribution/'];
        exec_name = 'b.out';
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
        
        
        if (ismember(4,tools_to_test))
            [ourcode_output] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks);
            all_our_code_outputs{count_files_tested} = ourcode_output;
        else
            all_our_code_outputs{count_files_tested} = nan;
		end
        
		%% here we run Canopy on the simulated data
		%File name for the Ancestree format file to transform and input to Canopy		
		%Here input_file will serve as the variable containing the full
		%file name for the Ancestree input file
		
		path_to_folder = [pwd_start, '/../Canopy-master/demo_code/'];
		
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
		if (ismember(5,tools_to_test))
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
        
        %plot(count_files_tested);
        %drawnow;
        
    end 
end


%% commpute average errors accross all test

% here we only test the first solution of the different solutions output
% all_ancestree_errors_array = []; 
% for i =1:length(all_ancestree_errors)
%     all_ancestree_errors_array = [all_ancestree_errors_array ,  all_ancestree_errors{i}{1}];
% end

