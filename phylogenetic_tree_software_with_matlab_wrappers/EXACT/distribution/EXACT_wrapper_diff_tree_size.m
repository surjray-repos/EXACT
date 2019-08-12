%% Matlab wrapper for EXACT
% Calls EXACT_wrapper.m

% INPUTS:
% F_reduced = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% error_rate = 0.5 * sqrt (estimated variance of the samples)
% min_tree_size, max_tree_size = minimum and maximum tree size to explore during inference
% path_to_folder = path to the folder where EXACT will create temporary files
% exec_name = name of the compiled EXACT executable
% cpu_gpu = computing architecture to use. Possible options: "cpu" (use single CPU core), "cpu_multithread" (use multiple CPU cores), "gpu" (use GPU)
% cost = likelihood function to score each tree: possible options: "cost1", "cost2", "cost3", "cost4".
% In "cost1" we minimize | F_hat - F |_frobenious subject to F = UM, M >= 0 sum(M) == 1
% In "cost2" we minimize | F_hat - F |_frobenious subject to F = UM, M >= 0 sum(M) <= 1
% In "cost3" we minimize | U^-1 F_hat - M |_frobenious subject to M >= 0 sum(M) == 1
% In "cost4" we minimize | U^-1 F_hat - M |_frobenious subject to M >= 0 sum(M) <= 1
% k_best = how many top k best trees we want as output
% gpu_id = ID of the GPU to use, if gpu architecture was chosen
% cpu_cores = number of CPU cores to use if multithreading was chosen
% To explore multiple CPUs, or GPUs, even on different computers, the wrapper allows us to specify which portion of the space of all possible tree the current function call is going to explore.
% This is controlled using the parameters num_devices and tree_subset
% num_devices = specifies in how many equal parts the space number of possible trees is being divided by
% tree_subset = which particular subset of the tree space gets will be explored by the current function call.
% So, for example, if num_devices = 6 and tree_subset = 4, then the current calling function will scan all of the trees in the 4th subset of the partion of the space into 6 equal parts.
% CUDA_thread_block = number of CUDA threads per block, if gpu architecture was chosen
% CUDA_blocks = number of CUDA thread blocks, if gpu architecture was chosen
% OUTPUTS:
% best_M is a matlab cell object having 7 components, namely, 
% where
%	best_M{1} = likelihood score of the best tree as computed by the BIC criterion
%	best_M{2} = Bayesian information criteria score
%	best_M{3} = adjacency matrix for the best tree. This is a directed tree. If we call this matrix T, the U = inv(I - T), where U appears in the PPM model as F = UM.
%	best_M{4} = recovered (clean) frequencies of mutations
%	best_M{5} = clustered frequencies of mutants
%	best_M{6} = cluster membership information. Array where the ith element indicates the cluster ID to which the ith mutation belongs to.
%	best_M{7} = run time (in seconds) for the executable to infer the best tree among all trees of the same size while keeping track of all k_best trees
% best_bic = best_M{2}
% all_Ms is a matlab cell object with k_best cells (indexed by sol_id), each cell storing 
%	all_Ms{sol_id}{1} = likelihood score of the sold_id tree
%	all_Ms{sol_id}{2} = Bayesian information criteria score
%	all_Ms{sol_id}{3} = adjacency matrix for the sold_id tree. This is a directed tree. If we can this matrix T, the U = inv(I - T), where U appears in the PPM model as F = UM.
%	all_Ms{sol_id}{4} = recovered (clean) frequencies of mutations
%	all_Ms{sol_id}{5} = clustered frequencies of mutations
%	all_Ms{sol_id}{6} = cluster membership information. Array where the ith element indicates the cluster ID to which the ith mutation belongs to.
%	all_Ms{sol_id}{7} = run time (in seconds) for the executable to infer the sold_id tree among all trees of the same size while keeping track of all k_best trees


function [best_M, best_bic, all_Ms] = EXACT_wrapper_diff_tree_size(F_reduced, error_rate, min_tree_size, max_tree_size, path_to_folder, exec_name, cpu_gpu, cost, k_best, gpu_id, cpu_cores, num_devices, tree_subset, CUDA_thread_block, CUDA_blocks )

    best_bic = inf;
    best_M = {};
	curr_sol = {};
    all_Ms = {};
	
    for tree_size = min_tree_size:max_tree_size

        [clusters_ix, clustered_Fs] = kmeans(F_reduced, tree_size);
        
        clustered_Fs = [ones(1,size(clustered_Fs,2),1) ; clustered_Fs];
        clusters_ix = clusters_ix + 1;
        
        [M_tmp, run_time] = EXACT_wrapper(clustered_Fs, path_to_folder, exec_name, cpu_gpu, cost, k_best, gpu_id, cpu_cores, num_devices, tree_subset, CUDA_thread_block, CUDA_blocks );
        
        for sol_id = 1:k_best
            costfunctionid = str2num(cost(5));
            [cost_value, Mut_freqs] = compute_tree_cost_in_cvx(clustered_Fs, inv(eye(tree_size+1) - M_tmp{sol_id}.tree), tree_size+1, size(clustered_Fs,2), costfunctionid);
            
            %Checking our tree's cost value vs cvx cost_value^2
            if (abs(cost_value * cost_value - M_tmp{sol_id}.val) > 0.001)
                error('There is an error between the error produced by Matlab and our C code');
            end

            F_clust = inv(eye(tree_size+1) - M_tmp{sol_id}.tree) * Mut_freqs;

            likelihood = cost_value*cost_value / (2.0 * error_rate * error_rate);
            bic = likelihood + size(clustered_Fs, 2) * (tree_size+1 - 1.0) * log(size(F_reduced, 1));

            curr_sol{sol_id} = {M_tmp{sol_id}.val, bic, M_tmp{sol_id}.tree, Mut_freqs, F_clust, clusters_ix, run_time};

            %all_Ms{tree_size - min_tree_size + 1} = curr_sol;
			all_Ms{tree_size - min_tree_size + 1}{sol_id} = curr_sol{sol_id};

            if (bic < best_bic)
                best_bic = bic;
                best_M = curr_sol{sol_id};
            end
        end

    end
end
