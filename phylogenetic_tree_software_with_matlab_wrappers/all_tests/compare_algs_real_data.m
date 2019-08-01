%load('/home/sam_safavi/Desktop/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/real_data_first_four_algs_errors_April_19.mat')


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