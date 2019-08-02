% this function indexes all ways to link elements backward into groups of size a(1), a(2), ..., a(k)
% there are prod_i a(i+1)^(a(i)) values for ix
% the ix variable can take values 1, 2, 3, 4, 5, ...
% the output groups will be concatenated. we now how to read them based on the size of a

function [groups] = index_a_links_with_k_groups(ix , a)
    
    % a is a vector for size k

    groups = nan(1 , sum(  a(2:end)  ) );

    % groups is a vector of size equal to the sum of the all the elements in a except the first

    k = length(a);

    ix = ix  - 1;

    r = 0;
    for i = 1:k-1
    
        sub_ix = 1 + mod( ix  , a(i)^(a(i+1)) );
        
        [labels] = index_n_labels_k_elements(sub_ix , a(i) , a(i+1) );
        
        groups(1+r: r + a(i+1) ) = labels;
        
        ix = floor(ix / (   a(i)^(a(i+1))   ) );  % in C, when working with ints, you do not need the floor function. It rounds down by default.

        r = r + a(i+1);
    end



end