function [Guessed_Adj, Guessed_root, Guessed_AdjL, Guessed_T, Guessed_AdjT] = get_rooted_tree_from_index(tree_ix,n)

    [Guessed_Adj , ~, ~ , ~ , ~] = index_all_rooted_and_labeled_trees(tree_ix , n);

    Guessed_root = find(diag(   Guessed_Adj    ));

    Guessed_Adj = Guessed_Adj - diag(diag(Guessed_Adj));

    % get the adj list from the adj mat of the tree
    Guessed_AdjL = {};
    for i =1:n
        Guessed_AdjL{i} = find(Guessed_Adj(:,i));
    end

    % get adj list of tree
    Guessed_T = BFS(Guessed_AdjL,Guessed_root);

    % get adj mat of tree
    Guessed_AdjT = zeros(n);
    for i =1:n
        for j = Guessed_T{i}
            Guessed_AdjT(i,j) = 1;
        end
    end


end