%% Matlab wrapper for AncesTree
% calls the EXACT executable to infer phylogenetic trees

% INPUTS:
% F_reduced = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% scale = multiplying factor to transform the mutation frequencies back to read counts i.e. read counts = F_reduced * scale
% Desription of alpha, beta and gamma from https://github.com/raphael-group/AncesTree
% alpha = Controls the clustering of mutations in the graph clustering phase: only arcs (v_j, v_k) with 0.5 - alpha <= min_p P(X_pj < X_pk) <= 0.5 + alpha are considered
% beta = Controls the confidence in ancestral relationships in the graph: there is an arc (v_j, v_k) if min_p P(X_pj < X_pk) >= beta
% gamma = Controls the allowed pertubation of observed variant frequencies by defining (1 - gamma) confidence intervals
% wrapper_dir = path to the folder where AncesTree will create temporary files and folders
% pathtoprogram = full path to the AncesTree executable
% OUTPUTS:
% M contains structures with the output trees and cluster information, mutant frequencies and calculated F values
% M being a matlab cell object having 7 components, namely, 
%	M{1} = recovered (clean) frequencies of mutations
%	M{2} = clustered frequencies of mutants
%	M{3} = adjacency matrix for the best tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
%	M{4} = cluster membership information in form of a cell array, each row cell designates which particular cluster/node that group of mutations belongs to.
%	M{5} = input complete set of frequencies of mutations
%	M{6} = pre-clustering assignment of mutations. An array with 2 columns, the 1st column designating the cluster/node, and the 2nd column designating the mutation that belongs to that cluster
% M{1}, M{2}, M{3} and M{4} might have, inside them, just one or several cells depending on how many solutions AncesTree has inferred.	


function [M] = ancestree_wrapper(F_reduced, scale, alpha, beta, gamma, wrapper_dir, pathtoprogram )

    starting_directory = pwd;

    disp('The virtual root is node 0');

    n = size(F_reduced,1);
    T = size(F_reduced,2);

    rng shuffle;   
    file_number = num2str(randi(1000000000));
    fileF = ['fileF', file_number];
    disp(fileF);
    
    %Making sure we are in the path of Ancestree
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to AncesTree root is: %s\n', pwd )
    
    %Creating a unique test_directory name for input to AncesTree
    test_directory = ['data', '_', fileF, '/'];
    command_to_exec = ['mkdir ', test_directory];
    fprintf( 'making test_directory: %s\n', command_to_exec);
    system(command_to_exec);
    wrapper_dir_unique = [wrapper_dir, test_directory]
    cd(wrapper_dir_unique);
    fprintf( 'The new current path, after changing to unique directory is: %s\n', pwd )

    fileoutput = [wrapper_dir, test_directory, 'vec_output.sol'];   
    
    fileF_AncesTree = [wrapper_dir, test_directory, fileF, '.txt'];
    fID = fopen(fileF_AncesTree,'w');
    
    % the formula is F_according_to_them = right_counts / (right_counts + left_counts)
    % however, F_according_to_them is equal to 0.5*F_in_reality
    % so F_in_reality = 2 *  right_counts / (right_counts + left_counts)
    right = round((0.5 * F_reduced) * scale);
    left = scale - right;
    
    bigtable = zeros(n,2*T);
    
    bigtable(:,1:2:end) = left;
    bigtable(:,2:2:end) = right;
    
    one = 1;
        
    for t = 0:n
        if (t == 0)
            fprintf(fID,'gene_id\t');
            for i = 1:T
                if ( i < T )
                    fprintf(fID,'time_%d\ttime_%d\t',i,i);
                elseif ( i == T )
                    fprintf(fID,'time_%d\ttime_%d',i,i);
                end
            end
        else
            fprintf(fID,'mut_%d\t',t);
            for i = 1:T*2
                %Checking if this is the rightmost column or not since
                %rightmost column should not end with a tab /t
                if ( i < (T*2) )
                    %Checking if the value to be inserted is 0 and if so
                    %change it to 1
                    if ( round(bigtable(t,i)) ~= 0 )
                        fprintf(fID,'%d\t',round(bigtable(t,i)));
                    elseif ( round(bigtable(t,i)) == 0 ) 
                        fprintf(fID,'%d\t', 0*one );
                    end    
                elseif ( i == (T*2) )
                    if ( round(bigtable(t,i)) ~= 0 )
                        fprintf(fID,'%d',round(bigtable(t,i)));
                    elseif ( round(bigtable(t,i)) == 0 ) 
                        fprintf(fID,'%d', 0*one );
                    end 
                end
            end
        end
        fprintf(fID,'\n');
    end
    fclose(fID);
    
    fprintf('alpha = %f\n', alpha);
    whos alpha;
    fprintf('beta = %f\n', beta);
    whos beta;
    fprintf('gamma = %f\n', gamma);
    whos gamma;
    
    %Reconstructing the file name containing the F input 
    command_to_exec1 = [pathtoprogram,' --version'];
    command_to_exec2 = [pathtoprogram ,' --alpha ',num2str(alpha), ' --beta ',num2str(beta), ' --gamma ',num2str(gamma), ' ', fileF_AncesTree, ' --sol ', fileoutput];
    disp('Ancestree command: ');
    disp(command_to_exec2);
    
    %Start of Try-end block right before AncesTree command execution
    %try
    tic();
    system(command_to_exec1);
    system(command_to_exec2);
    runtime_Ancestree = toc();
    %f = msgbox(num2str(runtime_Ancestree));
    
    %This output has been temporarily set to the .sol solution file. In the
    %future this will be the CPLEX ILP matrix
    foutputID = fopen(fileoutput,'r');
    
    line = fgetl(foutputID); % read line by line

    num_sols = sscanf(line, '%d');
    
    line = fgetl(foutputID); % read line by line
    line = fgetl(foutputID); % read line by line
    num_samples = sscanf(line, '%d');
    line = fgetl(foutputID); % read line by line
    num_mutations = sscanf(line, '%d');
    input_full_F = nan(num_samples, num_mutations);
    for i = 1:num_samples
        line = fgetl(foutputID); % read line by line
        input_full_F(i,:) = str2num(line);
    end
    line = fgetl(foutputID); % read line by line
    line = fgetl(foutputID); % read line by line
    line = fgetl(foutputID); % read line by line
    line = fgetl(foutputID); % read line by line

    all_Fs = cell(num_sols,1);
    all_Ms = cell(num_sols,1);
    all_Ts = cell(num_sols,1);
    all_clusts = cell(num_sols,1);
    for i = 1:num_sols
        line = fgetl(foutputID); % read line by line
        num_samples_of_sol = sscanf(line, '%d');
        line = fgetl(foutputID); % read line by line
        num_mutations_of_sol = sscanf(line, '%d');
        M_of_sol = nan(num_samples_of_sol, num_mutations_of_sol);
        for j = 1:num_samples_of_sol
            line = fgetl(foutputID); % read line by line
            M_of_sol(j,:) = str2num(line);
        end
        all_Ms{i} = M_of_sol;
        line = fgetl(foutputID); % read line by line
        U_size_x = sscanf(line, '%d');
        line = fgetl(foutputID); % read line by line
        U_size_y = sscanf(line, '%d');
        U_of_sol = nan(U_size_x, U_size_y);
        for j = 1:U_size_x
            line = fgetl(foutputID); % read line by line
            U_of_sol(j,:) = str2num(line);
        end
        all_Ts{i} = eye(U_size_x) - inv(U_of_sol');
        line = fgetl(foutputID); % read line by line
        line = fgetl(foutputID); % read line by line
        line = fgetl(foutputID); % read line by line
        num_samples_of_sol = sscanf(line, '%d');
        line = fgetl(foutputID); % read line by line
        num_mutations_of_sol = sscanf(line, '%d');
        F_of_sol = nan(num_samples_of_sol, num_mutations_of_sol);
        for j = 1:num_samples_of_sol
            line = fgetl(foutputID); % read line by line
            F_of_sol(j,:) = str2num(line);
        end
        all_Fs{i} = 2*F_of_sol;
        line = fgetl(foutputID); % read line by line
        line = fgetl(foutputID); % read line by line
        
        cluster_string_line = fgetl(foutputID); % read line by line
        cluster_member_tokens = strsplit(cluster_string_line, ' ')';

        cluster_of_sol = cell(num_mutations_of_sol,1);
        for j=1:num_mutations_of_sol
            c_row = strsplit(cluster_member_tokens{j}, ';')
            [~, element_number] = size(c_row);

            %To loop over every element on c_row and put them into
            %genes_mat[]
            for k=1:element_number
                cluster_row_element = 1 + str2num(c_row{k}); % we do the plus one because we want the indices of the mutations to start at 1
                cluster_of_sol{j} =  [ cluster_of_sol{j}; cluster_row_element  ];
            end
        end
        line = fgetl(foutputID); % read line by line
        all_clusts{i} = cluster_of_sol;
    end
    
    fclose(foutputID);
    
    file_solution_cluster_genes = [wrapper_dir_unique, 'solution_cluster_file.txt'];
    fID_solution_cluster_genes = fopen(file_solution_cluster_genes,'r');
    pre_clust_assign = fscanf(fID_solution_cluster_genes,' value of i = %d value of it row number = %d');
    pre_clust_assign = 1 + reshape(pre_clust_assign,2,length(pre_clust_assign)/2)';
    fclose(fID_solution_cluster_genes);
    
    M = {all_Fs, all_Ms, all_Ts, all_clusts, input_full_F, pre_clust_assign};
    
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to Ancestree root is: %s\n', pwd )
    
    %removing all the directories and temporary files generated from this
    %run
    %Cleaning up the /data_ directory from this very run/process
    command_to_exec = ['rm -r data_', fileF];
    fprintf( '%s\n', command_to_exec);
    system(command_to_exec);
    
    cd(starting_directory);

end
