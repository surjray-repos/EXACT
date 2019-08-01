%% this is a matlab wrapper for the functions of CITUP. 
% It transforms the input in a form similar to Ancestree and then calls the
% CITUP clustering.py to cluster the input mutations
% the root is node 0

% M{1} stores multiple adjcencay matrices. One per solution found. These
% solutions are AFTER clustering has been done. ALl the trees here have the
% same corresponding clustering map. The root of the trees is always node
% zero.
% M{2} cluster that takes real mutations into clustered mutations ID
% M{3} output M matrices for each of the solutions in M{1}
% M{4} output F matrices for each of the solutions in M{1}
% M{5}  F matrix for the clustered mutations. THis F matrix is a function
% of the input unclustered F matrix and the particular mutation done by the
% program
% M[6] correlation of mutations to nodes

function [M] = CITUP_wrapper(F_reduced, wrapper_dir, cluster_number, pathtoprogram, error_rate  )
    
    disp('The root is node 0');
    
    starting_directory = pwd;

    top_tree_number = cluster_number;
    all_best_trees_flag = 0; %0(false) or 1(true) flag, if 1

    n = size(F_reduced,1);
    T = size(F_reduced,2);

    %% this is a matlab wrapper for transforming F_reduced count matrix to CITUP format
    rng shuffle;    
    file_number = num2str(randi(1000000000));
    fileF = ['fileF', file_number];
    disp(fileF);
    
    %Creating a unique test_directory name for input to CITUP
    test_directory = [wrapper_dir, 'test', '_', fileF, '/'];
    command_to_exec = ['mkdir ', test_directory];
    fprintf( 'making test_directory: %s\n', command_to_exec);
    system(command_to_exec);
    
    fileF_CITUP = [test_directory, 'Frequencies', '_', fileF, '.txt'];
    pause(0.1);
    fID_CITUP = fopen(fileF_CITUP,'w');
    
    delimiterIn = '\t';
    headerlinesIn = 1;
    input_mat = F_reduced;
    
    data_table = input_mat;
    [T,n] = size(data_table);
    numberNodes = n;
    
    numberSamples = numberNodes;
    numberMutations = T;
    fprintf('For CITUP input format:\n');
    fprintf('The number of samples for CITUP = %d\n', numberSamples);
    fprintf('The number of mutations for CITUP = %d\n', numberMutations);
    
    mutation_row_array = zeros(1,numberNodes);
    
    %Printing the number of samples and mutations at beginning of file    
    fprintf(fID_CITUP,'Num_mutations: %d\n',numberMutations);
    fprintf(fID_CITUP,'Num_samples: %d\n',numberSamples);
    %error_rate = 0.03;
    fprintf(fID_CITUP,'Error_rate: %f\n',error_rate);
    
    %Printing two blank lines after Error rate in the file
    fprintf(fID_CITUP,'\n\n');
    
    for t = 1:numberMutations
        for i = 1:numberSamples
            mut_freq_CITUP = data_table(t,i);
            mutation_row_array(:,i) = mut_freq_CITUP;
        end
        
        %writing the mutation_row_array to the file
        for i = 1:numberSamples-1
            fprintf(fID_CITUP,'%f ',mutation_row_array(:,i));
        end
        fprintf(fID_CITUP,'%f ',mutation_row_array(:,numberNodes));
        
        if ( t < T )
            fprintf(fID_CITUP,'\n');
        end
    end
    
    fprintf( 'The current path is: %s\n', pwd );
    CITUP_binary_path = [wrapper_dir, 'bin/'];
    cd(CITUP_binary_path);
    fprintf( 'The current path, after changing directory to CITUP bin is: %s\n', pwd );
    
    command_to_exec = [pathtoprogram, 'runCITUP2.sh ', test_directory(1:end-1), ' ', int2str(cluster_number), ' ', fileF];
    fprintf( 'CITUP command: %s\n', command_to_exec);
    
    %Start of Try-end block right before CITUP command execution
    %try
    tic();
    system(command_to_exec);
    runtime_CITUP = toc();
    %f = msgbox(num2str(runtime_CITUP));
    
    %Reconstituting the test_directory string structure
    test_directory = ['test', '_', fileF, '/'];
    file_CITUP_cluster = [wrapper_dir, test_directory, 'cluster_membership', '_', fileF, '.txt'];
    pause(1);
    fid_CITUP_cluster = fopen(file_CITUP_cluster, 'r');
    
    genes_mat = nan(numberMutations,2); %a structure to contain row numbers, and cluster designations
    
    s1 = fgets(fid_CITUP_cluster);
    cluster_designation = textscan(s1, '%d');
    
    for k = 1:numberMutations
        genes_mat(k,1) = k;
        genes_mat(k,2) = cluster_designation{1}(k);
	end
    
    %Opening the f_hat_solution file and reading in the f_hat matrix
    file_F_hat = [wrapper_dir, test_directory, 'f_hat_solution', '_', fileF, '.txt'];
    pause(0.1);
    F_hat_fid = fopen(file_F_hat, 'r');
    
    F_hat_mat = [];
    while(1)
        line = fgetl(F_hat_fid);
        if (line == -1)
            break;
        end
        line_F_hat = str2num( line );
        F_hat_mat = [F_hat_mat; line_F_hat];
	end
    
    %Checking if we are in the CITUP root directory, for the
    %visualizeResults
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to CITUP is: %s\n', pwd );
    
    %Removing the older version of the results_summary directory
    command_to_exec = ['rm -r results_summary/'];
    fprintf( 'results_summary removal: %s\n', command_to_exec);
    system(command_to_exec);
    
    CITUP_binary_path = [wrapper_dir, 'bin/'];
    cd(CITUP_binary_path);
    fprintf( 'The current path, after changing directory to CITUP bin is: %s\n', pwd );
    
    if ( all_best_trees_flag == 0 )
        command_to_exec = ['perl ./visualizeResults.pl ', '-r ', wrapper_dir, test_directory, ' ', '-g ', wrapper_dir, 'GammaAdjMatrices ', '-o ', wrapper_dir, 'results_summary_', fileF, ' ', '-n ', int2str(top_tree_number)];
        fprintf( 'CITUP command: %s\n', command_to_exec);
    else
        command_to_exec = ['perl ./visualizeResults.pl ', '-r ', wrapper_dir, test_directory, ' ', '-g ', wrapper_dir, 'GammaAdjMatrices ', '-o ', wrapper_dir, 'results_summary_', fileF];
        fprintf( 'CITUP command: %s\n', command_to_exec);
    end
        
    system(command_to_exec);
    
    % for directory navigation in the dot subdirectory of CITUP
    % and import all trees into adjacency matrix format.
    
    fprintf( 'The current path is: %s\n', pwd );
    tree_dir_results_summary = [wrapper_dir, 'results_summary_', fileF, '/dots/'];
    cd(tree_dir_results_summary);
    fprintf( 'The path, after changing to the results_summary, is: %s\n', pwd );
    
    tree_dir_pattern = [tree_dir_results_summary, '*.dot'];
    mat_tree_files = dir(fullfile(tree_dir_pattern));
    
    [number_dot_files, temp] = size(mat_tree_files);
    
    unique_count_num_sol = 1;
    unique_count_num_sol_M = 1;
    unique_count_num_sol_F = 1;
    sol_cell = {};
    F_sol_cell = {};
    M_sol_cell = {};
    
    %List of the indices of the best trees will be stored
    list_of_best_tree_indices_in_string_format = [];
    for i = 1:number_dot_files
        
        pause(0.1);
        tree_fid = fopen(mat_tree_files(i).name);
        
        % this finds the indices of the best trees
        
        k1 = strfind(mat_tree_files(i).name,'tree');
        k2 = strfind(mat_tree_files(i).name,'.');
        k2 = k2(1);
        list_of_best_tree_indices_in_string_format = [list_of_best_tree_indices_in_string_format, str2num(mat_tree_files(i).name(k1+4:k2-1)) ];
        
        ancestery_nodes = [];
        mat = [];
        
        while (~feof(tree_fid) )
            s_nodes = fgetl(tree_fid);
            
            pattern = '\d+ -> \d+;';
            test_pattern = regexp(s_nodes, pattern, 'match');
            
            if ( isempty(test_pattern) == 0 )
                %to extract node numbers from cell, first convert to char array
                ancestry_nodes = char(test_pattern);
                parent_child_array = sscanf(ancestry_nodes, '%d -> %d');
                
                parent_node = parent_child_array(1);
                child_node = parent_child_array(2);
                
                parent_node = parent_node + 1;
                child_node = child_node + 1;
                mat(parent_node,child_node) = 1;
                mat(child_node,parent_node) = 1;
            end
		end
        
        %Identify if this matrix is the same as others before it, and then
        %store it in the M cell array if unique
        unique_flag = 0;
        unique_counter = unique_count_num_sol - 1;
        if ( unique_counter ~= 0 )
            for j = 1:unique_counter
                tf = isequal(sol_cell{j}, mat);
                if (tf == 1)
                    unique_flag = 1;
                    break;
                end
            end
        end
        
        if (unique_flag == 0)
            sol_cell{unique_count_num_sol} = mat;
            unique_count_num_sol = unique_count_num_sol + 1;
		end  
        
        %Reopen the .dot file to read the alpha and q frequency values
        pause(0.1);
        tree_fid2 = fopen(mat_tree_files(i).name);
        
        while(1)
            line = fgetl(tree_fid2);
            if (strcmp(line(1:min(5,length(line))),'aFreq'))
                break;
            end
        end
        
        alpha_M_mat = [];
        while(1)
            line = fgetl(tree_fid2);
            lineM = str2num( line(2:end) );
            alpha_M_mat = [alpha_M_mat; lineM];
            if (strcmp(line(1:min(5,length(line))),'qFreq'))
                break;
            end
		end
        %Identify if this matrix is the same as others before it, and then
        %store it in the M cell array if unique
        unique_flag_M = 0;
        unique_counter_M = unique_count_num_sol_M - 1;
        if ( unique_counter_M ~= 0 )
            for j = 1:unique_counter_M
                tf_M = isequal(M_sol_cell{j}, alpha_M_mat);
                if (tf == 1)
                    unique_flag_M = 1;
                    break;
                end
            end
        end
        
        if (unique_flag_M == 0)
            M_sol_cell{unique_count_num_sol_M} = alpha_M_mat;
            unique_count_num_sol_M = unique_count_num_sol_M + 1;
		end
        
        beta_F_mat = [];
        while(1)
            line = fgetl(tree_fid2);
            if (line == -1)
                break;
            end
            lineM = str2num( line(2:end) );
            beta_F_mat = [beta_F_mat; lineM];
		end		
        %Identify if this matrix is the same as others before it, and then
        %store it in the M cell array if unique
        unique_flag_F = 0;
        unique_counter_F = unique_count_num_sol_F - 1;
        if ( unique_counter_F ~= 0 )
            for j = 1:unique_counter_F
                tf_F = isequal(F_sol_cell{j}, beta_F_mat);
                if (tf == 1)
                    unique_flag_F = 1;
                    break;
                end
            end
        end
        
        if (unique_flag_F == 0)
            F_sol_cell{unique_count_num_sol_F} = beta_F_mat;
            unique_count_num_sol_F = unique_count_num_sol_F + 1;
		end
        
    end
    
    %Reading the corresponding top results_ files and getting correlation of mutations to
    %nodes and storing in list_of_mut_to_nodes cell.
    list_of_best_tree_indices_in_string_format = unique(list_of_best_tree_indices_in_string_format);
    number_unique_best_trees = length(list_of_best_tree_indices_in_string_format);
    
    list_of_mut_to_nodes_cell = cell(number_unique_best_trees, 1);
    list_of_objective_val_cell = cell(number_unique_best_trees, 1);
    for r = 1:number_unique_best_trees
        id_str = num2str(list_of_best_tree_indices_in_string_format(r));
        file_F_res = [wrapper_dir, test_directory, 'Results_', fileF, '__',id_str,'.txt'];
        pause(0.1);
        F_res_fid = fopen(file_F_res, 'r');
        
        %Getting the first line which contains the objective value
        tmp = fgetl(F_res_fid);
        objective_val = sscanf(tmp, ' Objective_value %f');
        list_of_objective_val_cell{r} = objective_val;
        
        while(1)
            tmp = fgetl(F_res_fid);
            if ( strcmp(tmp,'cluster_centers') )
                break;
            end
        end
        tmp = fgetl(F_res_fid);
        tmp = fgetl(F_res_fid);
        list_of_mut_to_nodes = str2num(tmp);
        fclose(F_res_fid);
        list_of_mut_to_nodes_cell{r} = list_of_mut_to_nodes;
    end

    M = {sol_cell, genes_mat, M_sol_cell, F_sol_cell, F_hat_mat, list_of_mut_to_nodes_cell, list_of_objective_val_cell};
    
    fclose(fID_CITUP);
    fclose(fid_CITUP_cluster);
    fclose(tree_fid);
    fclose(tree_fid2);
    fclose(F_hat_fid);
    
    %End of Try-end block right after CITUP command execution
    %end
    
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to CITUP root is: %s\n', pwd )
    
    %removing all the directories and temporary files generated from this
    %run
    %Cleaning up the /test directory from this very run/process
    command_to_exec = ['rm -r test_', fileF];
    system(command_to_exec);
    %Removing the newer version of the results_summary directory
    command_to_exec = ['rm -r results_summary_', fileF];
    fprintf( 'results_summary removal: %s\n', command_to_exec);
    system(command_to_exec);
    
    cd(starting_directory);
end    