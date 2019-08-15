%% Function to find the smallest BIC(Bayesian Information criterion) for
% multiple cluster sizes
% calls CITUP_wrapper.m
% INPUTS:
% F_from_SampleData = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% wrapper_working_directory = path to the folder where CITUP will create temporary files and folders
% min_cluster_no, max_cluster_no = minimum and maximum tree size that CITUP will explore during inference
% CITUP_executable_path = full path of the CITUP R scripts and executables
% Description of error_rate from https://github.com/sfu-compbio/citup
% citup_error_rate = used to determine how noisy the allelic frequencies are. For real data, this value is set low (0.03-0.05) for deep sequencing and higher (0.05-0.08) for lower coverage datasets.
% OUTPUTS:
% sol_min_bic = minimum value of the BIC score, signifying most optimal tree
% best_sol is a matlab cell object having 5 components, containg the best tree found using a
% quadratic fitting cost and a Bayesian Information Criterion to select the
% output tree size
%	best_sol{1}{1} = adjacency matrix for the best output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
%	best_sol{2} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
%	best_sol{3}{1} = clustered frequencies of mutants. Rows are associated with samples, and columns with mutants
%	best_sol{4}{1} = recovered (clean) frequencies of clustered mutations. Rows are associated with samples, and columns with mutations
%	Note that best_sol{4}{1}' = U*best_sol{3}{1}'
%	best_sol{5} = frequencies of clustered-mutations after the kmeans preclustering
% all_int_sol contains all the best tree found by CITUP for each tree size in the range min_clust_no and max_clust_no, indexed by tree_size_id here, each cell storing 
%	all_int_sol{tree_size_id}{1}{1} = adjacency matrix for the optimal output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
%	all_int_sol{tree_size_id}{2} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
%	all_int_sol{tree_size_id}{3}{1} =  clustered frequencies of mutants. Rows are associated with samples, and columns with mutants
%	all_int_sol{tree_size_id}{4}{1} = recovered (clean) frequencies of clustered mutations. Rows are associated with samples, and columns with mutations
%	all_int_sol{tree_size_id}{5} =  frequencies of clustered-mutations after the kmeans preclustering	


function [sol_min_bic, best_sol, all_int_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate )
    error_rate = citup_error_rate;
    num_mutations = size(F_from_SampleData,1);
    num_samples = size(F_from_SampleData,2);
    
    sol_min_bic = inf;
    best_sol = {};
    all_int_sol = {};
    
    %testing the CITUP wrapper run:

    count_sol_process = 0;
    for cluster_no = min_cluster_no:max_cluster_no
        count_sol_process = count_sol_process +1;
        disp(['trying clusters of size ',num2str(cluster_no)]);
        %Running the CITUP wrapper with different cluster sizes
        [citup_output_BIC] = CITUP_wrapper(F_from_SampleData, wrapper_working_directory, cluster_no, CITUP_executable_path, citup_error_rate);
        all_int_sol{count_sol_process} = citup_output_BIC;
        
        %Getting the number of solutions
        [num_solutions, ~] = size(citup_output_BIC{7});
        for n_sol = 1:num_solutions
            %Getting the number of nodes in the corresponding solution tree
            [num_nodes, ~] = size(citup_output_BIC{1}{n_sol});
            if (num_nodes ~= cluster_no )
                error('num_nodes ~= cluster_no');
            end
            
            objective_value = citup_output_BIC{7}{n_sol};
            
            %BIC formula from pickBestTree.py
            likelihood = objective_value / (2.0 * error_rate * error_rate);
            bic = likelihood + num_samples * (num_nodes - 1.0) * log(num_mutations);
            disp([bic,n_sol, cluster_no]);
            
            if( bic < sol_min_bic )
                sol_min_bic = bic;
                best_sol = citup_output_BIC;
            end
        end
    end
    
   
end