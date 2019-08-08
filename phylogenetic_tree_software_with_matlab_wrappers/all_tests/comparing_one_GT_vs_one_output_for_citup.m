%%Comparing ground truth and output of Canopy across all 90 files
%Compares all U matrices and clusters

all_error_of_ancestry_relations = [ ];

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
        citup_output = load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/citup_all_files_ELKEBIR_simulated_data.mat','all_citup_outputs');
        citup_output = citup_output.all_citup_outputs{file_count_id};



        U1 = true_tree_data{3}';
        clust1 = true_tree_data{5};


        % all_citup_outputs{1}{1}{1}{1}
        % all_citup_outputs{1}{1}{2}

        U2 = citup_output{1}{1}{1};
        root_citup = 1; % the root of the treee is always node 1.
        n_virt_nodes_citup = size(U2,1);
        AdjLrecon = {};
        for j =1:n_virt_nodes_citupc  
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

        [error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);
        disp(error_rates)
        all_error_of_ancestry_relations = [all_error_of_ancestry_relations , error_rates(2)];
        plot(all_error_of_ancestry_relations, '*');
		drawnow;
	end
end