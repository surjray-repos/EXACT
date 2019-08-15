%% Function to find the smallest BIC(Bayesian Information criterion) for
% multiple cluster sizes


function [sol_min_bic, best_sol,all_int_sol] = CITUP_wrapper_diff_clust_size( F_from_SampleData, wrapper_working_directory, min_cluster_no, max_cluster_no, CITUP_executable_path, citup_error_rate )
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