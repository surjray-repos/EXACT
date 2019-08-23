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

%% Function to plot the inferred trees from PhyloWGS output along with their respective fitting cost

% INPUTS:
% phylowgs_output is a matlab cell object having 3 components.
% Each component is a cell object indexed by sol_id, which lists different good solutions.
% The number of solutions that PhyloWGS outputs is given by how many different sol_id indices there are in the output 
% phylowgs_output{1}{sol_id} = adjacency matrix for the sol_id th output tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
% phylowgs_output{2}{sol_id} = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
% phylowgs_output{3}{sol_id} = frequencies of the input, after mutations have been clustered, helps in obtaining clustered frequencies of mutants.


function generate_list_of_PhyloWGS_trees(phylowgs_output)

	figure;
	
	num_solutions = size(phylowgs_output{1}, 2);

	for sol_id = 1:num_solutions % this will show us all output trees
		Tree_Matrix_T = phylowgs_output{1}{sol_id};

		U = phylowgs_output{1}{sol_id};
		root_phylowgs = 1; % the root of the treee is always node 1.
		n_virt_nodes_phylowgs = size(U,1);
		AdjLrecon = {};
		for j =1:n_virt_nodes_phylowgs
			AdjLrecon{j} = find(U(:,j));
		end
		Treerecon = BFS(AdjLrecon,root_phylowgs); %root is the node 0 in unlabled tree
		% get adj mat of tree
		AdjTrecon = zeros(n_virt_nodes_phylowgs);
		for j =1:n_virt_nodes_phylowgs
			for r = Treerecon{j}
				AdjTrecon(j,r) = 1;
			end
		end
		U = inv(eye(n_virt_nodes_phylowgs) - AdjTrecon);
		num_mutants = size(phylowgs_output{2}{sol_id},1);

		clust = zeros(num_mutants,2);
		for i = 1:size(phylowgs_output{2}{sol_id},1)
			clust(1 + phylowgs_output{2}{sol_id}{i,1},:) = [1 + phylowgs_output{2}{sol_id}{i,1},1 + phylowgs_output{2}{sol_id}{i,2}];
		end

		Mutant_Frequencies_M = inv(U)*phylowgs_output{3}{sol_id};

		nodelbs = cell(1,length(U));
		for i = 1:length(U)
			nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
		end

		ax = subplot(2, round(num_solutions/2), sol_id);
		h = plot(digraph(  eye(length(U)) - inv(U)  ));
		labelnode(h, [1:length(U)],nodelbs);
		if (sol_id == 1)
			title(ax, ['PhyloWGS list of different trees:']);
		end
		set(ax,'visible','off');
		set(findall(ax, 'type', 'text'), 'visible', 'on');
	end

	% maximize window
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end