%% this code tries to generate all the labeled trees with N vertices and a root node using an indexing scheme
% to check if the code is correct we just need to generate all the trees
% and see if we get n^(n-1) different trees for n nodes

n = 4;

num_trees = n^(n-1);

list_trees = [];

for tree_ix = 1 : num_trees
   
   [Adj , ~, ~ , ~ , ~] = index_all_rooted_and_labeled_trees(tree_ix , n);
   % Adj is a matrix stored in vector format
    
    list_trees = [list_trees ; Adj(:)'];
    
end



if (size(unique(list_trees,'rows'),1) == num_trees)
    disp('SUCCESS');
else
    disp('SOMETHING IS WRONG');
end

