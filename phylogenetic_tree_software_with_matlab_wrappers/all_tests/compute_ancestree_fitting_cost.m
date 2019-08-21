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