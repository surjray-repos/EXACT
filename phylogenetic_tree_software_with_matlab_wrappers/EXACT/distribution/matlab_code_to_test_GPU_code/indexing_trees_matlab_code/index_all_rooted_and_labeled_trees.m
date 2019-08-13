%% this code indexes all rooted and labeled trees on n nodes
% it outputs an adjencency matrix where root has a self-loop
% it also outputs a compact representation of the tree per level
% note that, for n nodes, there are n^(n-1) such trees
% the tree_ix takes values 1, 2, 3, ...

% Adj is an adjancecy matrix stored in vector format
% balls_level_config is a vector describing how many nodes per level we have
% groups is a vector describing which nodes per level we have
% links is a vector describing how the elements of each level link to the elements of the upper level


function [Adj , depth, balls_levels_config , groups , links] = index_all_rooted_and_labeled_trees(tree_ix , n)


    flag = 0;
    r = 0;
    for depth = 2:n
        for partition_ix = 1:nchoosek(  n-2  ,   depth - 2  )
       
            [balls_levels_config] = index_n_balls_d_levels_one_ball_1st_atleastone_later(partition_ix , n , depth);
            size_interval_balls_levels_config = count_n_multi_a_with_k_groups( n , balls_levels_config) * count_n_links_a_with_k_groups( n , balls_levels_config);

            r = r + size_interval_balls_levels_config;

            if (tree_ix <= r)
                flag = 1;
                break;
            end
        end 
        if (flag == 1)
            break;
        end
    end

    sub_ix = tree_ix - r + size_interval_balls_levels_config;
    
    sub_ix_1 = 1 + mod( sub_ix - 1  ,   count_n_multi_a_with_k_groups( n , balls_levels_config)  );
    [groups] = index_n_multi_a_with_k_groups(sub_ix_1 , n , balls_levels_config);
    
    sub_ix_2 = 1 + floor( (sub_ix - 1)  /   count_n_multi_a_with_k_groups( n , balls_levels_config)  );
    [links] = index_a_links_with_k_groups(sub_ix_2 , balls_levels_config);
 
    [Adj, ~] = get_adj_matrix_from_compact_repre(n , depth, balls_levels_config , groups , links  );
end