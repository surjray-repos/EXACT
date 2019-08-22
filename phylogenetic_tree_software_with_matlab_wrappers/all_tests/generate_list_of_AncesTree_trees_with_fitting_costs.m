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

%% Function to plot the inferred trees from Ancestree output along with their respective fitting cost

% INPUTS:
% ancestree_output is a matlab cell object having 6 components.
% The first four components are in turn a cell object indexed by sol_id, which lists different good solutions.
% The number of solutions that AncesTree outputs is given by how many different sol_id indices there are in the output
%	M{1}{sol_id} = recovered (clean) frequencies of clustered mutations.
%	Each row is associated to a different sample, and each column to a different cluster of mutations
%	M{2}{sol_id} = clustered frequencies of mutants
%	M{3}{sol_id} = adjacency matrix for the optimal tree. This is a directed tree. If we have this matrix T, then U = inv(I - T), where U appears in the PPM model as F = UM.
%	M{4}{sol_id} = cluster membership information for the clustering
%	associated to M{2}{sol_id} in the form of a cell array. The ith cell is
%	an array that lists the nodes that belong to the ith cluster. These clusters are a subset of the clusters in M{6}
%	
% The last two components are not indexed by sol_id, and are the same for all of the solutions that AncesTree outputs
%	M{5} = F_reduced'/2
%	M{6} = pre-clustering of mutations. An array with 2 columns, the 1st column designating the cluster ID, and the 2nd column designating the mutation that belongs to that cluster
% F_from_SampleData = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 


function generate_list_of_AncesTree_trees_with_fitting_costs(ancestree_output, F_from_SampleData)

	figure;

	num_solutions = size(ancestree_output{1}, 1);

	for sol_id = 1:num_solutions % this will show us all output trees
		Mutant_Frequencies_M = ancestree_output{2}{sol_id}';

		U = ancestree_output{3}{sol_id};
		U = inv(eye(length(U)) - U);
		clust = [];
		for i = 1:length(ancestree_output{4}{sol_id})
			for j = ancestree_output{4}{sol_id}{i}'
				clust = [clust; [j,i]];
			end
		end

		cost = compute_ancestree_fitting_cost(F_from_SampleData, clust, ancestree_output{1}{sol_id});

		nodelbs = cell(1,length(U));
		for i = 1:length(U)
			nodelbs{i} = num2str(clust(clust(:,2)==i,1)');
		end

		ax = subplot(2, num_solutions/2, sol_id);
		h = plot(digraph(  eye(length(U)) - inv(U)  ));
		labelnode(h, [1:length(U)],nodelbs);
		if (sol_id == 1)
			title(ax, ['AncesTree fitting cost of different trees: ' , num2str(cost)]);
		else
			title(ax, [num2str(cost)]);
		end
		set(ax,'visible','off');
		set(findall(ax, 'type', 'text'), 'visible', 'on');
	end
	
	% maximize window
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end