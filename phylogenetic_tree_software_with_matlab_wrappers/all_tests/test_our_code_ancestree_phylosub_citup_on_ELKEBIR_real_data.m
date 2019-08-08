%% add paths for different tools

addpath(genpath('../EXACT/'));
addpath(genpath('../AncesTree/'));
addpath(genpath('../citup'));
addpath(genpath('../phylowgs'));

%% process all the files

pwd_start = pwd;

start_directory = [pwd_start,'/Sample_test_data/AncesTree_data/real'];
all_folders = dir(start_directory);

all_ancestree_outputs = {};
all_citup_outputs = {};
all_phylosub_outputs = {};
all_our_code_outputs = {};

count_files_tested = 0;

tools_to_test = [4];

workspace_file_name = datestr(datetime);

all_input_files_in_folder = all_folders;

for file_ix = 3: size(all_input_files_in_folder,1)
        
    plot(count_files_tested);
    drawnow;


    count_files_tested = count_files_tested + 1;
    
    if ( count_files_tested >= 27 && count_files_tested <= 28)
        input_file_name = all_input_files_in_folder(file_ix).name;
        
        full_input_file_name = [start_directory,'/',input_file_name];
        disp(full_input_file_name);
        
        
        %% read input data from file
        input_file = full_input_file_name;
        [F_from_SampleData, scaling] =  transform_elkebir_input_data_into_F_matrix(input_file);
        
        n_mutations = size(F_from_SampleData,1);
        T_samples = size(F_from_SampleData,2);
        
        
        %% run ancestree on data
        
        scale = scaling;
        alpha = 0.3; % if alpha is big, lots of things will be clustered together
        beta = 0.8; % to choose a larger beta we need more samples , i.e. larger T_samples
        gamma = 0.01; % small gamma means larger confidence on the data
        wrapper_working_directory = [pwd_start, '/../AncesTree/distribution/'];
        ancestree_executable_path = [pwd_start, '/../AncesTree/distribution/build/ancestree'];
        
        % call ancestree
        %ancestree_output = ancestree_wrapper_ll(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
        %ancestree_output = ancestree_wrapper_ll_noTransform(input_file,  alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
        
        if (  ismember(1,tools_to_test) )
            ancestree_output = ancestree_wrapper(F_from_SampleData, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
            %ancestree_output = ancestree_wrapper_reading_only_sol_file_noTransform(F_from_SampleData, input_file, scale, alpha, beta, gamma, wrapper_working_directory, ancestree_executable_path);
            
            all_ancestree_outputs{count_files_tested} = ancestree_output;
        else
            all_ancestree_outputs{count_files_tested} = nan;
        end
        
        %% run CITUP on data
        
        wrapper_working_directory = [pwd_start,'/../citup/distribution/'];
        CITUP_executable_path = [pwd_start,'/../citup/distribution/bin/'];
        min_cluster_no = 5; % here we make the number of clusters be equal to the size of the input tree
        max_cluster_no = 9;
        %cluster_no = 8; % here we make the number of clusters be equal to the size of the input tree
        citup_error_rate = 0.03;
        
        %call our CITUP wrapper
        %[citup_output] = count_general_wrapper_CITUP_bowie2_OC_ll(F_from_SampleData, wrapper_working_directory, cluster_no, CITUP_executable_path, citup_error_rate);
        
        if (ismember(2,tools_to_test))
            %[~, citup_output] = count_general_wrapper_CITUP_BIC_ll( F_from_SampleData, wrapper_working_directory, 5,8, CITUP_executable_path, citup_error_rate )
            [~, citup_output, all_inter_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate )
            %all_citup_outputs{count_files_tested} = citup_output;
            all_citup_outputs{count_files_tested} = {citup_output, all_inter_sol};
        else
            all_citup_outputs{count_files_tested} = nan;
        end
        
        
        
        %% here we run phyloWGS on the data
        scale = scaling;
        
        %Constants for PhyloWGS, should not be changed for sequencing data
        %one = 1;
        mu_r = 0.999;
        mu_v = 0.5;
        wrapper_working_directory = [pwd_start,'/../phylowgs/'];
        phylowgs_executable_path = [pwd_start,'/../phylowgs/'];
        
        if (ismember(3,tools_to_test))
            [phylowgs_output] = count_general_wrapper_PhyloWGS_allow_duplicates_ll(F_from_SampleData, scale, wrapper_working_directory, phylowgs_executable_path, mu_r, mu_v);
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
        min_tree_size = 6;
        max_tree_size = 10;
        num_CPU_cores = 50;
        error_rate = 0.03;
        max_num_partitions = 1;
        device_tree_subset_value = 1;
        CUDA_threads_per_block = 32;
        CUDA_blocks = 128;
        top_k_value = 100;
        
        
        if (ismember(4,tools_to_test))
            [best_M, best_bic, all_Ms] = EXACT_wrapper_diff_tree_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks );
            %[best_M, best_bic, all_Ms] = wrapper_for_brute_force_phylo_with_pre_clust_auto_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks );

			all_our_code_outputs{count_files_tested} = [best_M, best_bic, all_Ms];
        else
            all_our_code_outputs{count_files_tested} = nan;
        end
        
        %% here we run Canopy on the real data
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
		cluster_number_end = 9; % we use 9 for all files except 90
		
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
end
