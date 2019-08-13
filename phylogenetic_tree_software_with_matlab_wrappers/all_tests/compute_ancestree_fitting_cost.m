function cost = compute_ancestree_fitting_cost(F_input, clust, clean_F)

	cost = 0;
	
	for t = 1:size(F_input,2)
		for i = 1:size(F_input,1)
			if (~isempty(    clust(clust(:,1) == i,2)    ))
				cost = cost + abs(F_input(i,t) - clean_F( t,  clust(clust(:,1) == i,2)  ));
			end
		end
	end

	cost = cost / (size(F_input,2) * size(F_input,1));

end