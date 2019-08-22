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

%% Function to plot the inferred trees from CITUP output along with their respective fitting cost

% INPUTS:
% all_inter_sol contains all the best tree found by CITUP for each tree size in the range min_clust_no and max_clust_no, indexed by tree_size_id here, each cell storing 
%	all_inter_sol{tree_size_id}{1}{1} = adjacency matrix for the optimal output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
%	all_inter_sol{tree_size_id}{2} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
%	all_inter_sol{tree_size_id}{3}{1} =  clustered frequencies of mutants. Rows are associated with samples, and columns with mutants
%	all_inter_sol{tree_size_id}{4}{1} = recovered (clean) frequencies of clustered mutations. Rows are associated with samples, and columns with mutations
%	all_inter_sol{tree_size_id}{5} =  frequencies of clustered-mutations after the kmeans preclustering	
% min_cluster_no = minimum tree size that CITUP will explore during inference


function generate_list_of_CITUP_trees_with_fitting_costs(all_inter_sol, min_cluster_no)

	figure;

	num_sizes_tested = size(all_inter_sol, 2);

	for tree_size_id = 1:num_sizes_tested % this will show us all output trees

		U = all_inter_sol{tree_size_id}{1}{1};
		root_citup = 1; % the root of the treee is always node 1.
		n_virt_nodes_citup = size(U,1);
		AdjLrecon = {};
		for j =1:n_virt_nodes_citup
			AdjLrecon{j} = find(U(:,j));
		end
		Treerecon = BFS(AdjLrecon,root_citup); %root is the node 0 in unlabled tree
		% get adj mat of tree
		AdjTrecon = zeros(n_virt_nodes_citup);
		for j =1:n_virt_nodes_citup
			for r = Treerecon{j}
				AdjTrecon(j,r) = 1;
			end
		end
		U = inv(eye(n_virt_nodes_citup) - AdjTrecon);
		clust = all_inter_sol{tree_size_id}{2};
		clust(:,2) = 1 + clust(:,2);

		cost = norm(all_inter_sol{tree_size_id}{5}' - all_inter_sol{tree_size_id}{4}{1} ,'fro')^2; 

		nodelbs = cell(1,length(U));
		for i = 1:length(U)
			nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
		end

		ax = subplot(2, num_sizes_tested/2, tree_size_id);
		h = plot(digraph(  eye(length(U)) - inv(U)  ));
		labelnode(h, [1:length(U)],nodelbs);
		if (tree_size_id == 1)
			title(ax, ['CITUP (fitting cost, tree size) for the best tree for each size : (' , num2str(cost), ',', num2str(min_cluster_no + tree_size_id - 1),')']);
		else
			title(ax,['(' , num2str(cost), ',', num2str(min_cluster_no + tree_size_id - 1),')']);
		end
		set(ax,'visible','off');
		set(findall(ax, 'type', 'text'), 'visible', 'on');
	end

	% maximize window
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end