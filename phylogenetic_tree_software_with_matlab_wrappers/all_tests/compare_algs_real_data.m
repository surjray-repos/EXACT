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

%% Calculates of standard deviation among the results of the four error types for the various tools
% The tools have been run on the the simulated data provided in the
% AncesTree paper. The results have been contrasted against ground truth
% provided in the same paper to calculate the errors.

function compare_algs_real_data[phylowgs_exact_mean_error, phylowgs_exact_mean_error_std, phylowgs_citup_mean_error, phylowgs_citup_mean_error_std, phylowgs_ancestree_mean_error, phylowgs_ancestree_mean_error_std, exact_citup_mean_error, exact_citup_mean_error_std, exact_ancestree_mean_error, exact_ancestree_mean_error_std, citup_ancestree_mean_error, citup_ancestree_mean_error_std, phylowgs_canopy_mean_error, phylowgs_canopy_mean_error_std, exact_canopy_mean_error, exact_canopy_mean_error_std, citup_canopy_mean_error, citup_canopy_mean_error_std, ancestree_canopy_mean_error, ancestree_canopy_mean_error_std] = compute_means_and_std_of_error_for_all_tools_vs_all_tools(all_four_errors_phylowgs_exact, all_four_errors_phylowgs_citup, all_four_errors_phylowgs_ancestree, all_four_errors_exact_citup, all_four_errors_exact_ancestree, all_four_errors_citup_ancestree, all_four_errors_phylowgs_canopy, all_four_errors_exact_canopy, all_four_errors_citup_canopy, all_four_errors_ancestree_canopy)

	phylowgs_exact_mean_error = mean(all_four_errors_phylowgs_exact,2);
	phylowgs_exact_mean_error_std = zeros(4,1);
	for i = 1:4
	phylowgs_exact_mean_error_std(i,1) = std(all_four_errors_phylowgs_exact(i,:))/sqrt(size(all_four_errors_phylowgs_exact,2));
	end

	phylowgs_citup_mean_error = mean(all_four_errors_phylowgs_citup,2);
	phylowgs_citup_mean_error_std = zeros(4,1);
	for i = 1:4
	phylowgs_citup_mean_error_std(i,1) = std(all_four_errors_phylowgs_citup(i,:))/sqrt(size(all_four_errors_phylowgs_citup,2));
	end

	phylowgs_ancestree_mean_error = mean(all_four_errors_phylowgs_ancestree,2);
	phylowgs_ancestree_mean_error_std = zeros(4,1);
	for i = 1:4
	phylowgs_ancestree_mean_error_std(i,1) = std(all_four_errors_phylowgs_ancestree(i,:))/sqrt(size(all_four_errors_phylowgs_ancestree,2));
	end

	exact_citup_mean_error = mean(all_four_errors_exact_citup,2);
	exact_citup_mean_error_std = zeros(4,1);
	for i = 1:4
	exact_citup_mean_error_std(i,1) = std(all_four_errors_exact_citup(i,:))/sqrt(size(all_four_errors_exact_citup,2));
	end

	exact_ancestree_mean_error = mean(all_four_errors_exact_ancestree,2);
	exact_ancestree_mean_error_std = zeros(4,1);
	for i = 1:4
	exact_ancestree_mean_error_std(i,1) = std(all_four_errors_exact_ancestree(i,:))/sqrt(size(all_four_errors_exact_ancestree,2));
	end

	citup_ancestree_mean_error = mean(all_four_errors_citup_ancestree,2);
	citup_ancestree_mean_error_std = zeros(4,1);
	for i = 1:4
	citup_ancestree_mean_error_std(i,1) = std(all_four_errors_citup_ancestree(i,:))/sqrt(size(all_four_errors_citup_ancestree,2));
	end

	phylowgs_canopy_mean_error = mean(all_four_errors_phylowgs_canopy,2);
	phylowgs_canopy_mean_error_std = zeros(4,1);
	for i = 1:4
	phylowgs_canopy_mean_error_std(i,1) = std(all_four_errors_phylowgs_canopy(i,:))/sqrt(size(all_four_errors_phylowgs_canopy,2));
	end

	exact_canopy_mean_error = mean(all_four_errors_exact_canopy,2);
	exact_canopy_mean_error_std = zeros(4,1);
	for i = 1:4
	exact_canopy_mean_error_std(i,1) = std(all_four_errors_exact_canopy(i,:))/sqrt(size(all_four_errors_exact_canopy,2));
	end

	citup_canopy_mean_error = mean(all_four_errors_citup_canopy,2);
	citup_canopy_mean_error_std = zeros(4,1);
	for i = 1:4
	citup_canopy_mean_error_std(i,1) = std(all_four_errors_citup_canopy(i,:))/sqrt(size(all_four_errors_citup_canopy,2));
	end

	ancestree_canopy_mean_error = mean(all_four_errors_ancestree_canopy,2);
	ancestree_canopy_mean_error_std = zeros(4,1);
	for i = 1:4
	ancestree_canopy_mean_error_std(i,1) = std(all_four_errors_ancestree_canopy(i,:))/sqrt(size(all_four_errors_ancestree_canopy,2));
	end

end