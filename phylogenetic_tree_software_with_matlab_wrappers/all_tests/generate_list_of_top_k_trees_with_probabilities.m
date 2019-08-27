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

%% Function to plot the inferred trees from EXACT output along with their respective probabilities

% INPUT
% display_top_k = how many top scoring inferred trees we want as output in the plot
% top_k_value = how many top k best trees we want as output
% max_tree_size, min_tree_size = maximum and minimum tree size to explore during inference
% error_rate = 

function generate_list_of_top_k_trees_with_probabilities(display_top_k, top_k_value, max_tree_size, min_tree_size, error_rate, ourcode_all_Ms, ourcode_output)

	figure;

	tree_probs = [];
	sorted_ids = nan(top_k_value, 1 + max_tree_size - min_tree_size);
	sorted_vals = nan(top_k_value, 1 + max_tree_size - min_tree_size);

	for j = 1:1 + max_tree_size - min_tree_size
		for i = 1:top_k_value
			sorted_vals(i,j) = ourcode_all_Ms{j}{i}{1};
		end
		[sorted_vals(:,j), sorted_ids(:,j)] = sort(sorted_vals(:,j));
	end

	for t = 1:display_top_k % this will show us the top 10 trees
		best_size = size(ourcode_output{3}, 1); % we can choose which size to work with here
		best_tree_ix = sorted_ids(t, best_size - min_tree_size);
		best_ourcode_output = ourcode_all_Ms{best_size - min_tree_size}{best_tree_ix};

		% the probability decays fast, so that we can normalized by the sum over the first few trees
		tree_probs(t) = (exp(-best_ourcode_output{1}./(2*error_rate*error_rate)))/sum(exp(-sorted_vals(:,best_size - min_tree_size)./(2*error_rate*error_rate)));

		M_target = best_ourcode_output{4};

		U2 = inv(eye(size(best_ourcode_output{3})) - best_ourcode_output{3});
		clust2 = [ [1: length(best_ourcode_output{6})]' , best_ourcode_output{6}];

		nodelbs = cell(1,length(U2));
		for i = 1:length(U2)
			nodelbs{i} = num2str(clust2(clust2(:,2)==i,1)');
		end

		ax = subplot(2,5,t);
		h = plot(digraph(  eye(length(U2)) - inv(U2)      ));
		labelnode(h, [1:length(U2)],nodelbs);
		if (t==1)
			title(ax, ['List of probabilities: ' num2str(tree_probs(t))]);
		else
			title(ax, num2str(tree_probs(t)));
		end
		set(ax,'visible','off');
		set(findall(ax, 'type', 'text'), 'visible', 'on');

	end

	% maximize window
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end