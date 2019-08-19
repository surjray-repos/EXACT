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

%% Function to calculate 4 error types from U matrices and clustering information for inferred tree vs ground truth tree

% INPUTS:
% ground_truth_Tree_matrix_on_clustered_nodes = ancestry matrix, for the ground truth tree, showing parent-child relationships between mutants, U appears in the PPM model as F = UM.
% ground_truth_clusters = clustering information, for the ground truth tree, to show which mutations belong to which cluster/node.
% U_recon = ancestry matrix, for the inferred tree, showing parent-child relationships between mutants, U appears in the PPM model as F = UM.
% clusters_for_sol = clustering information, for the inferred tree, to show which mutations belong to which cluster/node.
% OUTPUTS:
% Details of computing the four error types are given in section 4.1 of our paper, "Exact inference under the perfect phylogeny model"
% when being contrasted with the ground truth tree, 
% error type 1 checks whether the mutations are correctly clustered or not, 
% error type 2 checks whether ancestral relations are being maintained correctly or not, 
% error type 3 checks whether nodes that are incomparable in the ground truth tree display similar behavior in the inferred tree or not, and 
% error type 4 checks whether all the input mutations show up in the inferred tree,
% i.e. whether the inferred tree considers all the available mutation information or a subset of the mutations.
% error_rates = error rates is an array object with 4 rows, each row containing the value of a particular error type. 


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

