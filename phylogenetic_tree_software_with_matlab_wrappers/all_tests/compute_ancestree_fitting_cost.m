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

%% Matlab wrapper to compute the fiting cost for inferred trees by AncesTree

% INPUTS: 
% F_input = matrix with frequency of mutation values, each row is associated with a mutated position, each column is associated with a sample. 
% clust = cluster membership information for the clustering. An array with 2 columns, the 2nd column designating the cluster ID, and the 1st column designating the mutation that belongs to that cluster
% clean_F = recovered (clean) frequencies of clustered mutations.
% OUTPUTS:
% cost = fitting cost for the respective inferred tree 


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