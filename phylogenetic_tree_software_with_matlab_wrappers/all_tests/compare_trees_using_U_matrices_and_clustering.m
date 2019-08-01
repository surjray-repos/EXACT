function [error_rates] = compare_trees_using_U_matrices_and_clustering(ground_truth_Tree_matrix_on_clustered_nodes, ground_truth_clusters, U_recon, clusters_for_sol)
    % compares observed tree with ground truth
    % go over one solution and compute the error

    all_relevant_muts = union(ground_truth_clusters(:,1),    clusters_for_sol(:,1));
    
   
    list_of_errors = [];
    for ix = 1:length(all_relevant_muts)
        for jx = ix+1:length(all_relevant_muts)
            
            i = all_relevant_muts(ix);
            j = all_relevant_muts(jx);
            
            i_gt_virt = ground_truth_clusters(ground_truth_clusters(:,1)==i,2);
            j_gt_virt = ground_truth_clusters(ground_truth_clusters(:,1)==j,2);
            
            ground_truth_flag = zeros(4,1); % 1 0 0 = clustered 0 1 0 = i_ances_j 0 -1 0 = j_ances_i 0 0 1 = incommp
            
            if (isempty(i_gt_virt) || isempty(j_gt_virt))
                ground_truth_flag(4) = 1;
            else
                if (i_gt_virt == j_gt_virt)
                    ground_truth_flag(1) = 1;
                else
                    if ( ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 1)
                        ground_truth_flag(2) = 1;
                    end
                    if ( ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 1)
                        ground_truth_flag(2) = -1;
                    end
                    if(ground_truth_Tree_matrix_on_clustered_nodes(i_gt_virt,j_gt_virt) == 0 && ground_truth_Tree_matrix_on_clustered_nodes(j_gt_virt,i_gt_virt) == 0)
                        ground_truth_flag(3) = 1;
                    end
                end
            end
            
          
            i_anc_virt = clusters_for_sol(clusters_for_sol(:,1)==i,2);
            j_anc_virt = clusters_for_sol(clusters_for_sol(:,1)==j,2);
            
            ances_flag = zeros(4,1);
            
             if (isempty(i_anc_virt) || isempty(j_anc_virt))
                ances_flag(4) = 1;
             else
                if (i_anc_virt==j_anc_virt)
                    ances_flag(1) =  1;
                else
                    if ( U_recon(i_anc_virt,j_anc_virt) == 1)
                        ances_flag(2) = 1;
                    end
                    if ( U_recon(j_anc_virt,i_anc_virt) == 1)
                        ances_flag(2) = -1;
                    end
                    if(U_recon(i_anc_virt,j_anc_virt) == 0 && U_recon(j_anc_virt,i_anc_virt) == 0)
                        ances_flag(3) = 1;
                    end
                end
            end

            list_of_errors =  [list_of_errors, [i;j;ances_flag ~= ground_truth_flag]];
            
        end
    end
    
    error_rates = mean(list_of_errors(3:6,:),2);

end

