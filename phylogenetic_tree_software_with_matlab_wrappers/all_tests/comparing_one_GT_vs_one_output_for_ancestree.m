%% Comparing ground truth and output of Canopy across all 90 files
%Compares all U matrices and clusters

all_error_of_ancestry_relations = [ ];

for file_count_id = 1:90

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
    ancestree_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/all_ancestree_outputs_on_synthetic_Elkebir_data_05-Jun-2018_14:08:10.mat','all_ancestree_outputs');
    ancestree_output = ancestree_output.all_ancestree_outputs{file_count_id};

    U1 = true_tree_data{3}';
    clust1 = true_tree_data{5};

    U2 = ancestree_output{3}{1};
    U2 = inv(eye(length(U2)) - U2);
    clust2 = [];
    for i = 1:length(ancestree_output{4}{1})
        for j = ancestree_output{4}{1}{i}'
        clust2 = [clust2; [j,i]];
        end
    end

    [error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);

    all_error_of_ancestry_relations = [all_error_of_ancestry_relations , error_rates(2)];
	plot(all_error_of_ancestry_relations, '*');
    drawnow;
end