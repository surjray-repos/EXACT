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