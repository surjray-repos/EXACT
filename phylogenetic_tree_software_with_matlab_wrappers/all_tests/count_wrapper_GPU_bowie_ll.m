%% this is a matlab wrapper for transforming the F_reduced matrix to our GPU code format
%Sample execution: M = count_wrapper_GPU_bowie(F, '/home/jbento/phylogenetic_tree_software_with_matlab_wrappers/our_GPU_code/', 'a.out', 'cpu_multithread', 'cost1', 4, 1, 32, 8, 3, 32, 128)

function [M, runtime_bruteForce] = count_wrapper_GPU_bowie_ll(F_reduced, path_to_folder, exec_name, cpu_gpu, cost, k_best, gpu_id, cpu_cores, num_devices, tree_subset, CUDA_thread_block, CUDA_blocks )
    
    path_to_program = [path_to_folder, exec_name, ];

    %dimensions of the input F matrix
    n = size(F_reduced,1);
    T = size(F_reduced,2);
    
    fprintf('The data of the input F matrix is:\n');
    disp(F_reduced);
    
    rng shuffle; 
    file_number = num2str(randi(1000000000));
    fprintf('file random number = %ld\n', file_number);
    fileF_GPU_code = [path_to_folder, 'fileF', file_number];  
    fID_GPU_code = fopen(fileF_GPU_code,'w');
    
    fprintf('The number of vertices in tree = %d\n', n);
    fprintf('The time limit for simulation = %d\n', T);
    
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
    disp( fileF );
    type (fileF);
    fprintf( 'The current path is: %s\n', pwd );
    cd(path_to_folder);
    fprintf( 'The current path, after changing directory to the GPU code is: %s\n', pwd );
    
    %Platform and cost function variables
    %platform = 'cpu_multithread';
    platform = cpu_gpu;
    %cost_function = 'cost1';
    cost_function = cost;
    %fileoutput = [path_to_folder, 'output', '_', platform, '_', cost_function, '_', num2str(numberNodes), '_', num2str(T), '.txt'];
    %k = 4; %best k trees
    k = num2str(k_best);
    fileoutput_k_best_trees = [path_to_folder,'output_',file_number, '_', platform, '_', cost_function, '_', num2str(numberNodes), '_', num2str(T), '_', num2str(k), '_', 'best_trees', '.txt'];
    %user_GPU_choice = 0;
    user_GPU_choice = num2str(gpu_id);
    %number_CPU_cores = 32;
    number_CPU_cores = num2str(cpu_cores);
    
    number_GPU_devices = num2str(num_devices);
    tree_subset_current_device = num2str(tree_subset);
    CUDA_threads_each_block = num2str(CUDA_thread_block);
    number_CUDA_blocks = num2str(CUDA_blocks);

    command_to_exec = [path_to_program, ' ', platform, ' ', cost_function, ' ', num2str(numberNodes_temp), ' ', num2str(T), ' ', fileF, ' ', fileoutput_k_best_trees, ' ', k, ' ', user_GPU_choice, ' ', number_CPU_cores, ' ', number_GPU_devices, ' ', tree_subset_current_device, ' ', CUDA_threads_each_block, ' ', number_CUDA_blocks];
    fprintf( 'GPU command: %s\n', command_to_exec )
    fprintf('The number of mutations = %d\n', numberNodes);
    fprintf('The time limit for simulation = %d\n', T);
    
    tic();
    system(command_to_exec);
    runtime_bruteForce = toc();
    
    %f = msgbox(num2str(runtime_bruteForce));
    
%     sol_cell = {};
%     count_num_sol = 1;
%     %Reading the output.txt file and extracting the adjacency list from it
%     fprintf( 'Opening the adjacency list file output by the GPU code\n' );
%     fID_output_trees = fopen(fileoutput,'r');
%     %Reading the first line, cost of this tree:
%     least_cost = fscanf(fID_output_trees, '%s:\n');
%     fprintf( 'The cost of this least cost tree = %s\n', least_cost );
%     text_line = fgets(fID_output_trees); %# read line by line
% 
%     ancestry_mat = [];
% 
%     while (~feof(fID_output_trees) )
%         fprintf( 'Tree solution number = %d\n', count_num_sol );
% 
%         text_line = fgets(fID_output_trees); %# read line by line
%         fprintf( 'text_line = %s\n', text_line );
% 
%         %Checking for emptiness of read line, MATLAB 2016b still has no ""
%         %if ( text_line ~= '' )
%         %   mat_parent_child = sscanf(text_line,'(%d,%d)');
%         %end
%         if ( isempty(text_line) == 0 )
%             mat_parent_child = sscanf(text_line,'(%d,%d)');
%         end
% 
%         disp ( mat_parent_child );
%         parent = mat_parent_child(1);
%         child = mat_parent_child(2);
%         ancestry_mat(parent, child) = 1;
%         ancestry_mat(child, parent) = 1;
% 
%     end   
%     fprintf( 'The ancestry matrix for tree %d is:\n', count_num_sol );
%     disp(ancestry_mat);
% 
%     sol_cell{count_num_sol} = ancestry_mat;
%     count_num_sol = count_num_sol + 1;

    %Reading the top k trees file and creating a cell of ancestry matrices
    top_k_trees_sol_cell = {};
    top_k_count_num_sol = 1;
    %Reading the _best_trees.txt.txt file and extracting the top_k adjacency list from it
    fprintf( 'Opening the top k trees file output by the GPU code\n' );

    %fileoutput_top_k_trees = ['/home/jbento/our_GPU_code/output_cpu_multithread_cost1_396_5_4_best_trees.txt'];
    fID_output_top_k_trees = fopen(fileoutput_k_best_trees,'r');

    %k = 4;
    %number_trees = k; %number of the best trees in the file, k best trees
    %fprintf( 'The number of best tree solutions in file = %d\n', number_trees );

    while (~feof(fID_output_top_k_trees) )
        fprintf( 'Tree solution number = %d\n', top_k_count_num_sol );

        %Reading the first line, cost of this tree:
        least_cost = fgets(fID_output_top_k_trees); %# read line by line
        fprintf( 'The cost of this least cost tree = %s\n', least_cost );

        top_k_trees_ancestry_mat = zeros(n,n);

        top_k_text_line = fgets(fID_output_top_k_trees); %# read line by line
        fprintf( 'text_line = %s\n', top_k_text_line );

        while ( isempty(top_k_text_line) == 0 )
            top_k_mat_parent_child = sscanf(top_k_text_line,'(%d,%d)');
            disp ( top_k_mat_parent_child );
            parent = top_k_mat_parent_child(2);
            child = top_k_mat_parent_child(1);
            top_k_trees_ancestry_mat(parent, child) = 1;
            %top_k_trees_ancestry_mat(child, parent) = 1;
            top_k_text_line = fgets(fID_output_top_k_trees); %# read line by line
            fprintf( 'text_line = %s\n', top_k_text_line );

            %If we have reached the end of this block
            block_end = strcmp( top_k_text_line, newline );
            if ( block_end == 1 )
                break
            end
        end 
        fprintf( 'The ancestry matrix for tree %d is:\n', top_k_count_num_sol );
        disp(top_k_trees_ancestry_mat);
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