% Modifications Copyright (c) 2019 Surjyendu Ray
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

%% Matlab wrapper for Canopy
% calls the Canopy R executables to infer phylogenetic trees

% INPUTS:
% full_input_file_name = full path of the input file with mutation counts, from which variant allele frequencies are calculated. "The input is a tab-separated ASCII text file. The first line contains the sample headers. The first column contains gene ids. Then every consecutive pair of columns contains read counts for reference alleles and alternate alleles, respectively."
% wrapper_dir = path to the folder where Canopy R scripts will create temporary files and folders
% Description of the followinginput variables from https://github.com/yuchaojiang/Canopy and the CRAN vignette of Canopy
% burnin_val = burnin of MCMC chains
% thin_val = MCMC chains thinning
% K_min_val, K_max_val = minimum and maximum number of subclones allowed respectively
% numchains_val = number of MCMC chains with random initiations
% maxsimrun_val, minsimrun_val = maximum and minimum number of simutation iterations for each chain, respectively
% writeskip_val = interval to store sampled trees
% cluster_start, cluster_end = determines the range of the number of mutation clusters (BIC as model selection metric)
% OUTPUTS:
% M_canopy_output is a matlab cell object having 7 components.
% M_canopy_output{1} = matrix representation of a tumor's clonal composition. Zsk is the indicator of whether the sth SNA is present at the kth clone. We first sort the columns, cluster the columns and the rows of the matrix to get unique relations between the mutants, and also add a null mutation
% M_canopy_output{2} = clustered frequencies of mutants
% M_canopy_output{3} = records rows of the unique_cols_rows_sorted_Z_with_null_mut matrix that got clustered together, in the form of a cell array. The ith cell is an array that lists the rows that belong to the ith cluster.
% M_canopy_output{4} = records columns of the unique_cols_rows_sorted_Z_with_null_mut matrix that got clustered together, in the form of a cell array. The ith cell is an array that lists the columns that belong to the ith cluster.
% M_canopy_output{5} = pre-clustering membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
% M_canopy_output{6} = post MCMC final cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
% M_canopy_output{7} = run time (in seconds) for the Canopy executable to infer the most optimal tree.


function [M_canopy_output] = canopy_wrapper(full_input_file_name, wrapper_dir, burnin_val, thin_val, K_min_val, K_max_val, numchains_val, maxsimrun_val, minsimrun_val, writeskip_val, cluster_start, cluster_end)

	%Script to create R_file,X_file, WM_file, Wm_file and Y_file from
	%Ancestree data
	
	% Default values of parameters from Canopy demo test code:
	% burnin = 10;
	% thin = 5;
	% K_min = 5;
	% K_max = 5;
	% numchains = 15;
	% maxsimrun = 100000;
	% minsimrun = 10000;
	% writeskip = 200;
	
	starting_directory = pwd;
	%Making sure we are in the path of Canopy
    fprintf( 'The current path is: %s\n', pwd );
    cd(wrapper_dir);
    fprintf( 'The current path, after changing directory to Canopy root is: %s\n', pwd )
	
	burnin = burnin_val;
	thin = thin_val;
	K_min = K_min_val;
	K_max = K_max_val;
	numchains = numchains_val;
	maxsimrun = maxsimrun_val;
	minsimrun = minsimrun_val;
	writeskip = writeskip_val;
	
	fID = fopen('parameters_file.tsv','w');
	fprintf(fID, '%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n',burnin, thin, K_min, K_max, numchains, maxsimrun, minsimrun, writeskip);
	fclose(fID);
	
	fID = fopen('cluster_parameters_file.tsv','w');
	fprintf(fID, '%d\n%d\n',cluster_start, cluster_end);
	fclose(fID);
	
	%getting the data from an Ancestry format input file for
	%R_matrix_from_SampleData and X_matrix_from_SampleData
	input_file = full_input_file_name;
	
	fID = fopen(input_file,'r');
	disp(fID);
	
	delimiterIn = '\t';
	headerlinesIn = 1;
	input_mat = importdata(input_file, delimiterIn, headerlinesIn);
	
	data_table = input_mat.data;
	[n,T] = size(data_table);
	numberNodes = n;
	samples = T/2;
	
	R_matrix_from_SampleData = [];
	X_matrix_from_SampleData = [];
	
	%Left set shall be the counts of the reference allele
	%right set shall be the counts of the alternate allele
	left_set = data_table(:,1:2:end);
	right_set = data_table(:,2:2:end);
	
	R_alt_allele = nan(1,samples);
	X_total_reads = nan(1,samples);
	
	for t = 1:numberNodes
		for i = 1:samples
			%R shall get the right_set elements (alternate alleles)
			R_alt_allele(i) = right_set(t,i);
			%X shall get the sum of right_set elements (alternate alleles) and
			%left_set
			X_total_reads(i) = left_set(t,i) + right_set(t,i);
		end
		
		R_matrix_from_SampleData = [R_matrix_from_SampleData; R_alt_allele];
		X_matrix_from_SampleData = [X_matrix_from_SampleData; X_total_reads];
	end
	
	fileF_R_matrix = ['R_file.tsv'];
	fID_R = fopen(fileF_R_matrix,'w');
	fileF_X_matrix = ['X_file.tsv'];
	fID_X = fopen(fileF_X_matrix,'w');
	
	%Writing the matrices to tab separated .tsv files
	for t = 0:numberNodes
		if (t == 0)
			%the column header row will have 'samples' number of columns
			for i = 1:samples
				if (i < samples) %last column will have a \n nextline
					fprintf(fID_R,'sample%d\t',i);
					fprintf(fID_X,'sample%d\t',i);
				else
					fprintf(fID_R,'sample%d\n',i);
					fprintf(fID_X,'sample%d\n',i);
				end
			end
		else
			%printing the row header at the beginning of each row
			fprintf(fID_R,'sna%d\t',t);
			fprintf(fID_X,'sna%d\t',t);
			
			%Checking if this is the rightmost column or not since
			%rightmost column should not end with a tab /t
			%writing R_matrix and X_matrix to file
			for i = 1:samples - 1
				fprintf(fID_R,'%d\t',R_matrix_from_SampleData(t,i));
				fprintf(fID_X,'%d\t',X_matrix_from_SampleData(t,i));
			end
			fprintf(fID_R,'%d',R_matrix_from_SampleData(t,samples));
			fprintf(fID_X,'%d',X_matrix_from_SampleData(t,samples));
			
			%rightmost column should not end with a tab /t
			if ( t < numberNodes )
				fprintf(fID_R,'\n');
				fprintf(fID_X,'\n');
			end
		end
	end
	
	%Creating the Y matrix and the Y_file.tsv
	Y_matrix = [];
	fileF_Y_matrix = ['Y_file.tsv'];
	fID_Y = fopen(fileF_Y_matrix,'w');
	
	%The number of columns will be 1 (non-cna) plus 3 cna columns since we are reusing
	%cna information from Canopy data.
	%Writing the Y switches to tab separated .tsv file
	one  = 1;
	zero = 0;
	for t = 0:numberNodes
		if (t == 0)
			%the column header row will have 1 plus 3 cna number of columns
			for i = 0:3
				if (i == 0) %first column will have the non-cna switches
					fprintf(fID_Y,'non-cna\t');
				elseif (i < 3)
					fprintf(fID_Y,'cna%d\t',i);
				else
					fprintf(fID_Y,'cna%d\n',i);
				end
			end
		else
			%printing the row header at the beginning of each row
			fprintf(fID_Y,'sna%d\t',t);
			
			%Checking if this is the rightmost column or not since
			%rightmost column should not end with a tab /t
			for i = 0:3
				if (i == 0)
					fprintf(fID_Y,'%d\t',one);
				elseif (i < 3)
					fprintf(fID_Y,'%d\t',zero);
				else
					fprintf(fID_Y,'%d',zero);
				end
			end
			
			%rightmost column should not end with a tab /t
			if ( t < numberNodes )
				fprintf(fID_Y,'\n');
			end
		end
	end
	fclose(fID_R);
	fclose(fID_X);
	fclose(fID_Y);
	
	% call the R program
	tic()
	system("/home/surjray/Downloads/R-3.5.1/bin/Rscript canopy_demo_toy_clustering.R");
	canopy_runtime = toc(); % runtime variable
	%f = msgbox(num2str(canopy_runtime));
	
	% read the output
	M_values = importdata('output_P_file.tsv', '\t', 1);
	M_values = M_values.data;
	raw_Z = importdata('output_Z_file.tsv', '\t', 1);
	raw_Z = raw_Z.data;
	
	% sort columns lexicographically
	[tmp, sorted_indices] = sortrows(flip(raw_Z,1)');
	sorted_raw_Z = flip(tmp',1);
	
	% cluster the columns
	[a, b, c] = unique(sorted_raw_Z','stable', 'rows');
	unique_cols_sorted_Z = a';
	clust_col_labels = c';
	
	% cluster the rows
	[a, b, c] = unique(unique_cols_sorted_Z,'stable', 'rows');
	unique_cols_rows_sorted_Z = a;
	clust_row_labels = c';
	
	% add the null mutation
	unique_cols_rows_sorted_Z_with_null_mut = [ones(1,size(unique_cols_rows_sorted_Z,2)); unique_cols_rows_sorted_Z];
	clust_row_labels = clust_row_labels +1; % the row added, the null mut, is clust id = 1;
	% get the clusters
	col_clust = cell(max(clust_col_labels),1);
	row_clust = cell(max(clust_row_labels),1);
	for i = 1:max(clust_col_labels)
		col_clust{i} = sorted_indices(find(clust_col_labels == i));
	end
	for i = 1:max(clust_row_labels)
		row_clust{i} = (find(clust_row_labels == i));
	end
	
	% compute new M matrix, based on the clustering obtained
	P_col_compress = [];
	for i = 1:max(clust_col_labels)
		P_col_compress(:,i) = sum(M_values(col_clust{i},:),1);
	end
	P_col_compress = P_col_compress';
	
	% read the output cluster values from the output_sna_cluster_file.tsv
	% file, preliminary clustering stage
	cluster_values = importdata('output_sna_cluster_file.tsv', '\t', 1);
	cluster_designation = cluster_values.data;
	
	%Getting the number of mutation rows in the cluster matrix
	[numberMutations, b] = size(cluster_designation);
	
	genes_mat = nan(numberMutations,2); %a structure to contain row numbers, and cluster designations
	
	for k = 2:numberMutations + 1
		mutation_ID = textscan(cluster_values.textdata{k}, '%d');
		mutation_ID_row_integer = mutation_ID{1};
		genes_mat(k - 1,1) = mutation_ID_row_integer;
		genes_mat(k - 1,2) = cluster_designation(k - 1, 1);
	end
	cluster_number = max(genes_mat(:,2));
	
	% read the output cluster values from the
	% toy_config_highest_likelihood_clusters.txt
	% file, MCMC run clustering output of clonal mutations
	genes_mat_mcmc = nan(numberMutations,2); %a structure to contain row numbers, and cluster designations
	counter = 1;
	
	cluster_values_mcmc = importdata('toy_config_highest_likelihood_clusters.txt', '\t');
	
	num_mcmc_clusters = length(cluster_values_mcmc);
	for k = 1:num_mcmc_clusters
		cluster_row_text = strsplit(cluster_values_mcmc{k});
		
		cluster_row_text_length = length(cluster_row_text);		
	
		for i = 1:cluster_row_text_length
			%Since the first element, from left, is always the mut
			%number, i.e. cluster number ID
			if ( i == 1 )
				clonal_mut = sscanf(cluster_row_text{i}, 'mut%d');
				clonal_mut_cluster_ID = clonal_mut;
			else
				sna_ID = sscanf(cluster_row_text{i}, 'sna%d');
				if ( isempty(sna_ID) ~= 1 )
					genes_mat_mcmc(counter, 1) = sna_ID;
					genes_mat_mcmc(counter, 2) = clonal_mut_cluster_ID;
					counter = counter + 1;
				end
			end
		end
		
	end
		
	% returning M with ancestry like compressed matrix, clustered M values and the clustering information
	M_canopy_output = {unique_cols_rows_sorted_Z_with_null_mut, P_col_compress, row_clust, col_clust, genes_mat, genes_mat_mcmc, canopy_runtime};
	
	%removing all the directories and temporary files generated from this
    %R Canopy run/process
    command_to_exec = ['rm parameters_file.tsv cluster_parameters_file.tsv R_file.tsv X_file.tsv Y_file.tsv output_P_file.tsv output_Z_file.tsv toy_config_highest_likelihood_clusters.txt'];
    fprintf( '%s\n', command_to_exec);
    system(command_to_exec);
    
    cd(starting_directory);
	
end