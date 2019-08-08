%%Comparing ground truth and output of EXACT across all 90 files
%Compares all U matrices and clusters

start_folder = pwd;

all_error_of_ancestry_relations = [ ];
all_our_code_outputs_new = {};

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

    %% test our code

%     path_to_folder = [pwd,'/../our_GPU_code/'];
%     exec_name = 'b.out';
%     cost_function = 'cost1';
%     cpu_gpu = 'gpu';
%     GPU_id = 0;
%     min_tree_size = 10;
%     max_tree_size = 10;
%     num_CPU_cores = 96;
%     error_rate = 0.03;
%     max_num_partitions = 1;
%     device_tree_subset_value = 1;
%     CUDA_threads_per_block = 32;
%     CUDA_blocks = 128;
%     top_k_value = 20;
% 
%     tic();
%     [ourcode_output_check] = wrapper_for_brute_force_phylo_with_pre_clust_auto_size(F_from_SampleData, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost_function, top_k_value, GPU_id, num_CPU_cores, max_num_partitions, device_tree_subset_value, CUDA_threads_per_block, CUDA_blocks );
%     runtime = toc();
% 
%     all_our_code_outputs_new{file_count_id} = ourcode_output_check;
%     
%     cd(start_folder);
    
    %% get correct form of matrices before comparing
    ourcode_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_all_files_ELKEBIR_synthetic_data_only_size_10_trees.mat','all_our_code_outputs_new');
    ourcode_output_check = ourcode_output.all_our_code_outputs_new{file_count_id};
    
    U2 = inv(eye(size(ourcode_output_check{3})) - ourcode_output_check{3} );
    clust2 = zeros(100,2);
    clust2(:,1) = 1:100;
    clust2(:,2) = ourcode_output_check{6};

    %%
    [error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);   
    all_error_of_ancestry_relations = [all_error_of_ancestry_relations , error_rates(2)];
    plot(all_error_of_ancestry_relations,'*');
    drawnow; 
end
