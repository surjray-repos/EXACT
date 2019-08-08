clc
clear
all_error_of_type_II = [ ];
all_four_errors = [ ];

alg_1 = 5;
alg_2 = 5;
for file_count_id = 1:36
    
    switch alg_1
        case 1
            output_1 = load('phyloWGS_all_files_ELKEBIR_real_data.mat','all_phylosub_outputs');
        case 2
            output_1 = load('our_code_on_real_data_tree_sizes_6_9_top_20.mat','all_our_code_outputs');
        case 3
            output_1 = load('citup_all_files_ELKEBIR_real_data.mat','all_citup_outputs');
        case 4
            output_1  = load('ancestree_all_files_minus_a_few_ELKEBIR_real_data.mat','all_ancestree_output');
        case 5
            output_1  = load('all_canopy_outputs_on_real_Elkebir_data.mat','all_canopy_outputs');
    end
    switch alg_2
        case 1
            output_2 = load('phyloWGS_all_files_ELKEBIR_real_data.mat','all_phylosub_outputs');
        case 2
            output_2 = load('our_code_on_real_data_tree_sizes_6_9_top_20.mat','all_our_code_outputs');
        case 3
            output_2 = load('citup_all_files_ELKEBIR_real_data.mat','all_citup_outputs');
        case 4
            output_2  = load('ancestree_all_files_minus_a_few_ELKEBIR_real_data.mat','all_ancestree_output');
        case 5
            output_2  = load('all_canopy_outputs_on_real_Elkebir_data.mat','all_canopy_outputs');
    end
    
    switch alg_1
        case 1
            output_1 = output_1.all_phylosub_outputs{file_count_id};
        case 2
            output_1 = output_1.all_our_code_outputs{file_count_id};
        case 3
            output_1  = output_1.all_citup_outputs{file_count_id};
        case 4
            output_1  = output_1.all_ancestree_output{file_count_id};
        case 5
            output_1  = output_1.all_canopy_outputs{file_count_id};
    end
    switch alg_2
        case 1
            output_2 = output_2.all_phylosub_outputs{file_count_id};
        case 2
            output_2 = output_2.all_our_code_outputs{file_count_id};
        case 3
            output_2 = output_2.all_citup_outputs{file_count_id};
        case 4
            output_2 = output_2.all_ancestree_output{file_count_id};
        case 5
            output_2  = output_2.all_canopy_outputs{file_count_id};
    end
    
    
    if (length(output_1) >= 2 && length(output_2) >= 2)
        
        switch alg_1
            case 1
                %PhyloWGS
                solution_id = 1;
                U1 = output_1{1}{solution_id};
                root_phyloWGS = 1; % the root of the treee is always node 1.
                n_virt_nodes_phyloWGS = size(U1,1);
                AdjLrecon = {};
                for j =1:n_virt_nodes_phyloWGS
                    AdjLrecon{j} = find(U1(:,j));
                end
                Treerecon = BFS(AdjLrecon,root_phyloWGS); %root is the node 0 in unlabled tree
                % get adj mat of tree
                AdjTrecon = zeros(n_virt_nodes_phyloWGS);
                for j =1:n_virt_nodes_phyloWGS
                    for r = Treerecon{j}
                        AdjTrecon(j,r) = 1;
                    end
                end
                U1 = inv(eye(n_virt_nodes_phyloWGS) - AdjTrecon);
                num_mutants = size(output_1{2}{solution_id},1);
                clust1 = zeros(num_mutants,2);
                for i = 1:size(output_1{2}{solution_id},1)
                    clust1(1+output_1{2}{solution_id}{i,1},:) = [1+output_1{2}{solution_id}{i,1},1+output_1{2}{solution_id}{i,2}];
                end
                
            case 2
                % EXACT
                U1 = inv(eye(size(output_1{3})) - output_1{3});
                clust1 = [ [1: length(output_1{6})]', output_1{6}];
            case 3
                %CITUP
                solution_id = 1;
                U1 = output_1{1}{1}{solution_id};
                root_citup = 1; % the root of the treee is always node 1.
                n_virt_nodes_citup = size(U1,1);
                AdjLrecon = {};
                for j =1:n_virt_nodes_citup
                    AdjLrecon{j} = find(U1(:,j));
                end
                Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
                % get adj mat of tree
                AdjTrecon = zeros(n_virt_nodes_citup);
                for j =1:n_virt_nodes_citup
                    for r = Treerecon{j}
                        AdjTrecon(j,r) = 1;
                    end
                end
                U1 = inv(eye(n_virt_nodes_citup) - AdjTrecon);
                clust1 = output_1{1}{2};
                clust1(:,2) = 1 + clust1(:,2);
                
            case 4
                %AncesTree
                solution_id = 1;
                U1 = output_1{3}{solution_id};
                U1 = inv(eye(length(U1)) - U1);
                clust1 = [];
                for i = 1:length(output_1{4}{solution_id})
                    for j = output_1{4}{solution_id}{i}'
                        clust1 = [clust1; [j,i]]; % the clusters array does not have the columns ordered. This is not a problem because the function   compare_trees_using_U_matrices_and_clustering   is compatible with this
                    end
                end
            case 5
                %Canopy
                [U1, clust1]= extract_U_mat_and_clust_from_canopy_output(output_1);
        end
        
        switch alg_2
            case 1
                %PhyloWGS
                solution_id = 1;
                U2 = output_2{1}{solution_id};
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
                num_mutants = size(output_2{2}{solution_id},1);
                clust2 = zeros(num_mutants,2);
                for i = 1:size(output_2{2}{solution_id},1)
                    clust2(1+output_2{2}{solution_id}{i,1},:) = [1+output_2{2}{solution_id}{i,1},1+output_2{2}{solution_id}{i,2}];
                end
            case 2
                % EXACT
                U2 = inv(eye(size(output_2{3})) - output_2{3});
                clust2 = [ [1: length(output_2{6})]' , output_2{6}];
            case 3
                %CITUP
                solution_id = 1;
                U2 = output_2{1}{1}{solution_id};
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
                clust2 = output_2{1}{2};
                clust2(:,2) = 1 + clust2(:,2);
            case 4
                %AncesTree
                solution_id = 1;
                U2 = output_2{3}{solution_id};
                U2 = inv(eye(length(U2)) - U2);
                clust2 = [];
                for i = 1:length(output_2{4}{solution_id})
                    for j = output_2{4}{solution_id}{i}'
                        clust2 = [clust2; [j,i]]; % the clusters array does not have the columns ordered. This is not a problem because the function   compare_trees_using_U_matrices_and_clustering   is compatible with this
                    end
                end
            case 5
                %Canopy
                [U2, clust2]= extract_U_mat_and_clust_from_canopy_output(output_2);
        end
        
        [error_rates] = compare_trees_using_U_matrices_and_clustering(U1, clust1, U2, clust2);
        disp(error_rates)
        all_error_of_type_II = [all_error_of_type_II , error_rates(2)];
        all_four_errors = [all_four_errors error_rates];
        
    end
    
    
end