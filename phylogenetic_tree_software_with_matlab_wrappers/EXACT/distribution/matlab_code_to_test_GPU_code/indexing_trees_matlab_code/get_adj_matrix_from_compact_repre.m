% this function builds an adjacency matrix from
% the compact representation we formed before
% to distinguish the root node we put a self loop at the root

function [Adj , root] = get_adj_matrix_from_compact_repre(n, depth, balls_levels_config, groups , links)

    Adj = zeros(n,n); % this is a matrix that will have to be stored in vector format

    root = groups(1);
    Adj(root,root) = 1;

    
    ixs = cumsum([0 , balls_levels_config]);   %this cumsum can do on the fly inside the for loop. You do not need to allocate memory to compute it.
    for l = 2:depth

        for i = 1+ixs(l):ixs(l+1)

            node1 = groups(i);

            node2 = groups( links(i - 1) + ixs(l-1)  );

            Adj(node1,node2) = 1;
            Adj(node2,node1) = 1;

        end

    end



end

