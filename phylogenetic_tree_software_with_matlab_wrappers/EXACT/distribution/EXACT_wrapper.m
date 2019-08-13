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

%% this is a matlab wrapper for transforming the F_reduced matrix to our GPU code format
% and calling the EXACT executable to infer phylogenetic trees
% INPUTS:
% F_reduced = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% path_to_folder = path to folder where EXACT will create temporary files
% exec_name = name of the compiled EXACT executable
% cpu_gpu = computing architecture to use. Possible options: "cpu" (use single CPU core), "cpu_multithread" (use multiple CPU cores), "gpu" (use GPU)
% cost = likelihood function to score each tree: possible options: "cost1", "cost2", "cost3", "cost4".
% k_best = how many top k best trees we want as output
% gpu_id = ID of the GPU to use, if gpu architecture was chosen
% cpu_cores = number of CPU cores to use if multithreading was chosen
% to explore multiple CPUs, or GPUs, even on different computers, the wrapper allows us to specify which portion of the space of all possible tree the current call is going to explore. This is controlled using num_devices and tree_subset
% num_devices = specifies in how many equal parts the space number of possible trees is being divided by
% tree_subset = which particular subset of the tree space gets will be explored by the current function call.
% CUDA_thread_block = number of CUDA threads per block, if gpu architecture was chosen
% CUDA_blocks = number of CUDA thread blocks, if gpu architecture was chosen
% OUTPUTS:
% M contains structures with the k_best tree/s and associated tree cost/s for the output tree/s 
%	M.tree = adjacency matrix for the best tree. This is a directed tree. If we can this matrix T, the U = inv(I - T), where U appears in the PPM model as F = UM.
%	M.val = cost of the solution tree
% runtime_bruteForce = run time (in seconds) for the executable to infer this tree

function [M, runtime_bruteForce] = EXACT_wrapper(F_reduced, path_to_folder, exec_name, cpu_gpu, cost, k_best, gpu_id, cpu_cores, num_devices, tree_subset, CUDA_thread_block, CUDA_blocks )
    
    path_to_program = [path_to_folder, exec_name, ];

    %dimensions of the input F matrix
    n = size(F_reduced,1);
    T = size(F_reduced,2);
    
    rng shuffle; 
    file_number = num2str(randi(1000000000));
    fileF_GPU_code = [path_to_folder, 'fileF', file_number];  
    fID_GPU_code = fopen(fileF_GPU_code,'w');
    
    for i = 1:n
        for t = 1:T
            fprintf( fID_GPU_code, '%f\n', F_reduced(i,t) );
        end
    end
    
    numberNodes = n;
    %For the GPU code we can only have trees upto 14 nodes, so for testing:
    if ( numberNodes > 13 )
        numberNodes_temp = 13;
    else
        numberNodes_temp = numberNodes;
    end
    
    %Reconstructing the file name containing the F input
    fileF = [path_to_folder, 'fileF', file_number];
    fprintf( 'The current path is: %s\n', pwd );
    cd(path_to_folder);
    fprintf( 'The current path, after changing directory to the GPU code is: %s\n', pwd );
    
    %Platform and cost function variables
    platform = cpu_gpu;
    cost_function = cost;   
    k = num2str(k_best); %best k trees
    fileoutput_k_best_trees = [path_to_folder,'output_',file_number, '_', platform, '_', cost_function, '_', num2str(numberNodes), '_', num2str(T), '_', num2str(k), '_', 'best_trees', '.txt'];
    user_GPU_choice = num2str(gpu_id);
    number_CPU_cores = num2str(cpu_cores);
    
    number_GPU_devices = num2str(num_devices);
    tree_subset_current_device = num2str(tree_subset);
    CUDA_threads_each_block = num2str(CUDA_thread_block);
    number_CUDA_blocks = num2str(CUDA_blocks);

    command_to_exec = [path_to_program, ' ', platform, ' ', cost_function, ' ', num2str(numberNodes_temp), ' ', num2str(T), ' ', fileF, ' ', fileoutput_k_best_trees, ' ', k, ' ', user_GPU_choice, ' ', number_CPU_cores, ' ', number_GPU_devices, ' ', tree_subset_current_device, ' ', CUDA_threads_each_block, ' ', number_CUDA_blocks];
    fprintf( 'GPU command: %s\n', command_to_exec )
    
    tic();
    system(command_to_exec);
    runtime_bruteForce = toc();
    %f = msgbox(num2str(runtime_bruteForce));

    %Reading the top k trees file and creating a cell of ancestry matrices
    top_k_trees_sol_cell = {};
    top_k_count_num_sol = 1;
    %Reading the _best_trees.txt.txt file and extracting the top_k adjacency list from it
    fID_output_top_k_trees = fopen(fileoutput_k_best_trees,'r');

    number_trees = k; %number of the best trees in the file, k best trees
	
    while (~feof(fID_output_top_k_trees) )
		
        %Reading the first line, cost of this tree:
        least_cost = fgets(fID_output_top_k_trees); %# read line by line

        top_k_trees_ancestry_mat = zeros(n,n);

        top_k_text_line = fgets(fID_output_top_k_trees); %# read line by line

        while ( isempty(top_k_text_line) == 0 )
            top_k_mat_parent_child = sscanf(top_k_text_line,'(%d,%d)');
            parent = top_k_mat_parent_child(2);
            child = top_k_mat_parent_child(1);
            top_k_trees_ancestry_mat(parent, child) = 1;
            top_k_text_line = fgets(fID_output_top_k_trees); %# read line by line

            %If we have reached the end of this block
            block_end = strcmp( top_k_text_line, newline );
            if ( block_end == 1 )
                break
            end
        end 
        top_k_trees_sol_cell{top_k_count_num_sol}.tree = top_k_trees_ancestry_mat;
        top_k_trees_sol_cell{top_k_count_num_sol}.val = str2num(least_cost);
        top_k_count_num_sol = top_k_count_num_sol + 1;

    end
    M = top_k_trees_sol_cell;

    fclose(fID_GPU_code);
    fclose(fID_output_top_k_trees);
    
    %Removing the input fileF and fileoutput_k_best_trees generated by this run
    command_to_exec = ['rm ', fileF];
    fprintf( 'fileF removal: %s\n', command_to_exec);
    system(command_to_exec);
    command_to_exec = ['rm ', fileoutput_k_best_trees];
    fprintf( 'fileoutput_k_best_trees removal: %s\n', command_to_exec);
    system(command_to_exec);
    
end