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

%% Function to generate tree diagrams, donut plot, Muller plot and error plots

% INPUTS:
% U = ancestry matrix showing parent-child relationships between mutants, U appears in the PPM model as F = UM.
% clust = clustering information to show which mutations belong to which cluster/node.
% Mutant_Frequencies_M = clustered frequencies of mutants
% Ugt = ancestry matrix, for the ground truth tree, showing parent-child relationships between mutants, U appears in the PPM model as F = UM.
% clustgt = clustering information, for the ground truth tree, to show which mutations belong to which cluster/node.

function generate_joint_plot_of_tree_donut_muller_and_errors(U, clust, Mutant_Frequencies_M, Ugt, clustgt)

	% build tree labels from the cluster information
	nodelbs = cell(1,length(U));
	for i = 1:length(U)
		nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
	end

	% compute colors for each mutant
	node_col = cell(1,length(U));
	cum_nodelbs = cell(1,length(U));
	for i = 1:length(U)
		for j = find(U(:,i))'
			cum_nodelbs{i} = [cum_nodelbs{i}, ' ', nodelbs{j}];
		end
		array_rep = [0, sort(str2num(cum_nodelbs{i}))];
		colormap_array = colormap(colorcube);

		rng(mod(  sum(((length(clust)).^[1:length(array_rep)]).*array_rep )    ,(2^32)-1)); % use a simple hash to generate a seed to then generate a random color for each mutant
		clrix = randi(64);
		node_col{i} = colormap_array(clrix,:); 
		cum_nodelbs{i} = num2str(sort(str2num(cum_nodelbs{i})));
	end
	%% Tree diagram with information about which nodes contain which mutations
	%f = figure;
	subplot(1,4,1);
	h = plot(digraph( eye(length(U)) - inv(U) ));
	title('First output tree');
	labelnode(h,[1:length(U)],nodelbs);
	highlight(h,[1:length(U)],'MarkerSize',20);
	for i = 1:length(U)
		highlight(h,i,'NodeColor',node_col{i});
	end
	set(gca,'visible','off');
	set(findall(gca, 'type', 'text'), 'visible', 'on');
	%% draw donut plot with mutant mixing ratios
	subplot(1,4,2);
	Mutant_Frequencies_M_num_size = size(Mutant_Frequencies_M', 1);
	[~, lgd] = donut(Mutant_Frequencies_M(:,1:Mutant_Frequencies_M_num_size)', cum_nodelbs, node_col, 1);
	title('Donut plot of mutants mixing ratios');
	set(lgd,'visible','off');
	set(gca,'visible','off');
	set(findall(gca, 'type', 'text'), 'visible', 'on');
	%% assuming that the different samples are obtained in time, we can draw a muller plot
	subplot(1,4,3);
	
	% checking if the M matrix has a similar number of rows as U
	% U is a square matrix with equal number of rows and columns
	ncol_U = size(U, 2);
	nrow_M = size(Mutant_Frequencies_M, 1);
	if (nrow_M < ncol_U)
		new_M = zeros(size(Mutant_Frequencies_M));
		new_M(size(Mutant_Frequencies_M) + 1,:) = 0;
		new_M(2:end, 1:end) = Mutant_Frequencies_M;
		Mutant_Frequencies_M = new_M;
	end
	[~, lgd] = generate_simple_muller_plots(U, Mutant_Frequencies_M, cum_nodelbs, node_col);
	title('Muller plot of mutants mixing ratios');
	set(lgd,'visible','off');
	set(gca,'visible','off');
	set(findall(gca, 'type', 'text'), 'visible', 'on');
	%% compute four different error types between synthetic data ground truth and EXACT inferred tree
	
	%Calculate the 4 different error types 
	[error_rates_EXACT] = compare_trees_using_U_matrices_and_clustering(Ugt, clustgt, U, clust);
	subplot(1,4,4);
	bar(error_rates_EXACT);
	title('Four error types comparison with ground truth');
	set(gca, 'xticklabel', {'Type I', 'Type 2', 'Type 3', 'Type 4'});
	set(findall(gca, 'type', 'text'), 'visible', 'on');
	%% maximize window
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

end