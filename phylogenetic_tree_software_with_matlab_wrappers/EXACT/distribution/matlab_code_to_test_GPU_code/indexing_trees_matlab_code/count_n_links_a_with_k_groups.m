% this function computes the number of ways link the elements for a(i+1)
% with the elements of a(i) for every i > 1

function [count] = count_n_links_a_with_k_groups( n , a)
    k = length(a);
    count = 1;
    for i = 1:k-1
        count = count *  ( a(i)^(a(i+1)) ) ;
    end
end