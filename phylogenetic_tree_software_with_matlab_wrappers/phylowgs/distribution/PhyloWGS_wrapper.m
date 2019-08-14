%% Matlab wrapper for PhyloWGS
% calls the PhyloWGS executable to infer phylogenetic trees

% INPUTS:
% F_reduced = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% scale = multiplying factor to transform the mutation frequencies back to read counts 
% wrapper_dir = path to the folder where PhyloWGS will create temporary files and folders
% phylowgs_exec_dir = full path to the directory where PhyloWGS executables are located
% Constants for PhyloWGS, for sequencing data from https://github.com/morrislab/phylowgs
% mu_r_val = fraction of expected reference allele sampling from the reference population
% mu_v_val = fraction of expected reference allele sampling from variant population
% OUTPUTS:
% M is a matlab cell object having 3 components.
% Each component is a cell object indexed by sol_id, which lists different good solutions.
% The number of solutions that PhyloWGS outputs is given by how many different sol_id indices there are in the output 
% M{1}{sol_id} = adjacency matrix for the sol_id th output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
% M{2}{sol_id} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
% M{3}{sol_id} = clustered frequencies of mutants.

function [M] = PhyloWGS_wrapper(F_reduced, scale, wrapper_dir, phylowgs_exec_dir, mu_r_val, mu_v_val)
close all;

    disp('The virtual root is node 0');
    
    starting_directory = pwd;
    
    n = size(F_reduced',1);
    T = size(F_reduced',2);
    
	% the virtual root is node 0
	
    rng shuffle;   
    file_number = num2str(randi(1000000000));
    fileF = ['fileF', file_number];
    
    %Making sure we are in the path of PhyloWGS
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to PhyloWGS root is: %s\n', pwd )
    
    %Creating a unique test_directory name for input to PhyloWGS
    test_directory = ['data', '_', fileF, '/'];
    command_to_exec = ['mkdir ', test_directory];
    fprintf( 'making test_directory: %s\n', command_to_exec);
    system(command_to_exec);
    wrapper_dir_unique = [wrapper_dir, test_directory]
    cd(wrapper_dir_unique);
    fprintf( 'The new current path, after changing to unique directory is: %s\n', pwd )
    
    fileF_PhyloWGS = [wrapper_dir, test_directory, fileF, '.txt'];
    fID = fopen(fileF_PhyloWGS,'w');
    
    right_table = round(0.5*F_reduced*scale);
    left_table = scale - right_table;
    
    first_half = zeros(1,n);
    second_half = zeros(1,n);
    
    one = 1;
    mu_r = mu_r_val;
    mu_v = mu_v_val;
        
    for t = 0:T
        if (t == 0)
            fprintf(fID,'id\t');
            fprintf(fID,'gene\t');
            fprintf(fID,'a\t');
            fprintf(fID,'d\t');
            fprintf(fID,'mu_r\t');
            fprintf(fID,'mu_v\n');
        else
            fprintf(fID,'s%d\t',t - 1);
            fprintf(fID,'gene_%d\t',t - 1);
            for i = 1:n
                %Checking if this is the rightmost column or not since
                %rightmost column should not end with a tab /t               
                if ( i <= n )
                    %Checking if the value to be inserted is 0 and if so
                    %change it to 1
                    %left_table elements will correspond to the first set
                    if ( round(left_table(t,i)) ~= 0 )
                        first_half(:,i) = round(left_table(t,i));
                    elseif ( round(left_table(t,i)) == 0 ) 
                        first_half(:,i) = one;                  
                    end  
                    
                    %right_table elements will correspond to the second set
                    if ( round(right_table(t,i)) ~= 0 )
                        second_half(:,i) = round(right_table(t,i) + first_half(:,i));
                    elseif ( round(right_table(t,i)) == 0 ) 
                        second_half(:,i) = round(1 + first_half(:,i));                 
                    end 
                end
            end
 
            %writing first_half and second_half to file
            for i = 1:n-1
                fprintf(fID,'%d,',first_half(:,i));
            end
            fprintf(fID,'%d\t',first_half(:,n));
        
            for i = 1:n-1
                fprintf(fID,'%d,',second_half(:,i));
            end
            fprintf(fID,'%d\t',second_half(:,n));
            
            fprintf(fID,'%.3f\t',mu_r);
            fprintf(fID,'%.2f',mu_v);
            
            if ( t < T )
                fprintf(fID,'\n');
            end
        end
    end
    fclose(fID);

    disp( fileF_PhyloWGS );    
    type (fileF_PhyloWGS);
    cnv_data_file = [wrapper_dir, 'cnv_data2.txt'];
    
    command_to_exec = ['python2 ', phylowgs_exec_dir, 'evolve.py -nS 1 ', fileF_PhyloWGS, ' ', cnv_data_file];
    fprintf( 'PhyloWGS command: %s\n', command_to_exec );
    
    %Start of Try-end block right before PhyloWGS command execution
    %try
    tic();
    system(command_to_exec);
    runtime_PhyloWGS = toc();
    %f = msgbox(num2str(runtime_PhyloWGS));
    
    file_k_trees = [wrapper_dir_unique, 'top_k_trees_adjacency_matrix.txt'];
    fID_k_trees = fopen(file_k_trees,'r');
    
    %n designates the number of nodes
    number_trees = 0;
    matrix_header = ['Adjacency_matrix:'];

    while ~feof(fID_k_trees)
        line = fgetl(fID_k_trees); %# read line by line
        disp( line );
   
        if ( strcmp(line, matrix_header) == 1 )
            number_trees = number_trees + 1;
        end
	end
    file_k_trees_max_element = [wrapper_dir_unique, 'top_k_trees_max_element.txt'];
    fid_dimension = fopen(file_k_trees_max_element);
    
    Number_nodes = fscanf(fid_dimension, '%d\n');
    sol_cell = {};
    count_num_sol = 1;
    unique_count_num_sol = 1;
    unique_trees = [];
    unique_clusters = [];
    unique_phi = [];
    
    fID_k_trees2 = fopen(file_k_trees,'r');
    while (~feof(fID_k_trees2) )
    
        dim = Number_nodes(count_num_sol,:);
    
        fscanf(fID_k_trees2, '%s:\n');

        mat = [];

        for i = 1:dim

            fscanf(fID_k_trees2, ' [');
            s1 = fscanf(fID_k_trees2, ' %f ', dim);

            fscanf(fID_k_trees2, ' ]\n');
            mat = [mat; s1'];
        end

        fscanf(fID_k_trees2, ' \n');
        fscanf(fID_k_trees2, ' \n');
    
        disp ( mat );
        
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
            unique_trees(unique_count_num_sol) = 1;
            unique_count_num_sol = unique_count_num_sol + 1;
        end
        %In this version we are storing all trees
        sol_cell{count_num_sol} = mat;
        count_num_sol = count_num_sol + 1;
    end
    
    disp( sol_cell );
    max_count = count_num_sol - 1;
    
    %transforming PhyloWGS cluster genes files to Matlab cells
    file_k_trees_genes = [wrapper_dir_unique, 'top_k_trees_max_element_genes.txt'];
    fID_k_trees_genes = fopen(file_k_trees_genes,'r');
    
    number_trees = 0;
    matrix_header = ['Nodes:'];
    
    while ~feof(fID_k_trees_genes)
        line = fgetl(fID_k_trees_genes); %# read line by line
        disp( line );
        
        if ( strcmp(line, matrix_header) == 1 )
            number_trees = number_trees + 1;
        end
	end
    
    file_k_trees_max_element = [wrapper_dir_unique, 'top_k_trees_max_element.txt'];
    fclose(fid_dimension);
    fid_dimension = fopen(file_k_trees_max_element);
    
    Number_nodes = fscanf(fid_dimension, '%d\n');
    
    sol_cell_genes = {};
    count_num_sol_genes = 1;
    unique_count_num_sol_genes = 1;
    
    fID_k_trees_genes2 = fopen(file_k_trees_genes,'r');
    while (~feof(fID_k_trees_genes2) )
        
        dim = Number_nodes(count_num_sol_genes,:);
        
        s_nodes = fgetl(fID_k_trees_genes2);
        
        %Matrix to capture the clustering data for the genes
        genes_mat = {};
        j = 1; %counter for the genes_mat matrix
        
        for i = 1:dim
            
            s1 = fgets(fID_k_trees_genes2);
            
            cluster_number = textscan(s1, '%d');
            
            pattern = 's\d+';
            test_pattern = regexp(s1, pattern, 'match');
            disp(test_pattern);
            [test_pattern_rows, test_pattern_columns] = size(test_pattern);
            cluster_member_no = test_pattern_columns;
            
            for k = 1:cluster_member_no
                mutation_position = sscanf(test_pattern{k},'s%d');
                genes_mat{j,1} = mutation_position;
                genes_mat{j,2} = cluster_number{1};
                j = j + 1;
            end
        end
        disp(genes_mat);
        
        %Identify if this matrix is the same as others before it, and then
        %store it in the M cell array if unique
        unique_flag = 0;
        unique_counter = unique_count_num_sol_genes - 1;
        if ( unique_counter ~= 0 )
            for m = 1:unique_counter
                tf = isequal(genes_cluster_cell{m}, genes_mat);
                if (tf == 1)
                    unique_flag = 1;
                    break;
                end
            end
        end
        
        if (unique_flag == 0)
            %genes_cluster_cell{unique_count_num_sol_genes} = genes_mat;
            unique_clusters(unique_count_num_sol_genes) = 1;
            unique_count_num_sol_genes = unique_count_num_sol_genes + 1;
        end
        
        fscanf(fID_k_trees_genes2, ' \n')
        
        %In this version we are storing all clusters
        genes_cluster_cell{count_num_sol_genes} = genes_mat;
        count_num_sol_genes = count_num_sol_genes + 1;
    end
    
    %Extracting the phi values from the top_k_trees file:
    file_k_trees_max_element = [wrapper_dir_unique, 'top_k_trees_max_element.txt'];
    file_k_trees_phi_values = [wrapper_dir_unique, 'top_k_trees_phi_values'];
    fclose(fid_dimension);
    fid_dimension = fopen(file_k_trees_max_element);
    
    Number_nodes = fscanf(fid_dimension, '%d\n');
    sol_cell_phi = {};
    count_num_sol = 1;
    unique_count_num_sol = 1;
    
    fID_k_trees_phi = fopen(file_k_trees_phi_values,'r');
    while (~feof(fID_k_trees_phi) )
        
        %Checking if we have reached the end of possible solutions
        if ( count_num_sol > max_count )
            break;
        end
        
        dim = Number_nodes(count_num_sol,:);
        
        fscanf(fID_k_trees_phi, '%s:\n');
        
        phi_mat = [];
        
        for i = 1:dim
            
            fscanf(fID_k_trees_phi, ' [');
            s1 = fscanf(fID_k_trees_phi, ' %f ', n);
            
            fscanf(fID_k_trees_phi, ' ]\n');
            phi_mat = [phi_mat; round(s1', 4)];
        end
        
        fscanf(fID_k_trees_phi, ' \n');
        fscanf(fID_k_trees_phi, ' \n');
        
        disp( phi_mat );
        
        %Identify if this matrix is the same as others before it, and then
        %store it in the M cell array if unique
        unique_flag = 0;
        unique_counter = unique_count_num_sol - 1;
        if ( unique_counter ~= 0 )
            for j = 1:unique_counter
                tf = isequal(sol_cell_phi{j}, phi_mat);
                if (tf == 1)
                    unique_flag = 1;
                    break;
                end
            end
        end
        
        if (unique_flag == 0)
            unique_phi(unique_count_num_sol) = 1;
            unique_count_num_sol = unique_count_num_sol + 1;
        end
        
        %In this version we are storing all phi_matrices
        sol_cell_phi{count_num_sol} = phi_mat;
        count_num_sol = count_num_sol + 1;      
    end

    M = {sol_cell, genes_cluster_cell, sol_cell_phi, unique_trees, unique_clusters, unique_phi};
    
    
    fclose(fID_k_trees);
    fclose(fID_k_trees_genes);
    fclose(fID_k_trees_genes2);
    fclose(fID_k_trees_phi);
    fclose(fID_k_trees2);
    fclose(fid_dimension);
    %delete(fileF);
    %End of Try-end block right after PhyloWGS command execution
    %end

    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to PhyloWGS root is: %s\n', pwd )
    
    %removing all the directories and temporary files generated from this
    %run
    %Cleaning up the /data_ directory from this very run/process
    command_to_exec = ['rm -r data_', fileF];
    system(command_to_exec);
    
    cd(starting_directory);
	
end