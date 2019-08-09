%% Wrapper to test our code with large data fro El.Kebir simulated+real F matrices
% with k-means clustering and checking against cvx value

% this adds a row of ones to the input
% the cluster index starts at 1, which is a virtual root.
% all the real data gets clustered into nodes 2 to tree_size in the virtual root

% INPUTS:
% F_reduced = matrix with frequency of mutation values, each row is a
%mutated position, each column is a sample or time-point.
% error_rate = error rate to calculate likelihood for BIC 
% min_tree_size, max_tree_size = minimum and maximum size, respectively, of
%inferred tree that we want. Determines number of clusters.
% path_to_folder = path to folder where EXACT will create temporary files
% exec_name = name of the compiled EXACT executable
% cpu_gpu = architecture: possible options: cpu, cpu_multithread, gpu
% cost = cost function to measure each tree: possible options: cost1, cost2, cost3, cost4
% k_best = how many top k best trees we want as output
% gpu_id = ID of the GPU to use, if gpu architecture was chosen
% cpu_cores = number of CPU cores if multithreading was chosen
% num_devices = number of partitions of the tree space for multithreading/parallelism
% tree_subset = Which particular subset of the tree space gets assigned to
%this particular device we are running on, in a parallel, workload, related to num_devices partitions above.
% CUDA_thread_block = number of CUDA threads per block we want if using a NVIDIA CUDA GPU
% CUDA_blocks = number of CUDA thread blocks we want if using a NVIDIA CUDA GPU
% OUTPUTS:
% best_M is {M_tmp{sol_id}.val, bic, M_tmp{sol_id}.tree, Mut_freqs, F_clust, clusters_ix, run_time} 
%contains the tree and cluster information for the tree with the best BIC score 
%	{M_tmp{sol_id}.val = cost of the best tree
%	bic = Bayesian information criteria score
%	M_tmp{sol_id}.tree = ancestry matrix for the best tree
%	Mut_freqs = calculated mutation frequency values from cvx
%	F_clust = clustered frequency values 
%	clusters_ix = cluster membership information, which cluster/node each mutation belongs to
%	run_time = run time (in seconds) for the executable to infer this tree
% best_bic = best bic value determining the best tree
% all_Ms = storing the {M_tmp{sol_id}.val, bic, M_tmp{sol_id}.tree, Mut_freqs, F_clust, clusters_ix, run_time} structure 
% for all the output solution trees for different sizes.

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