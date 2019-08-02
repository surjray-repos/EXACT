% this function computes the number of ways to create groups of sizes a(1), a(2) , a(3) , ... from n elements.

function [count] = count_n_multi_a_with_k_groups( n , a)
    % a is a vector of size k
    % there needs to be a function to access nchoosek from memory or something
    k = length(a);
    count = 1;
    r = 0;
    for i = 1:k
        count = count *  nchoosek( n-r   ,  a(i) ) ;
        r = r + a(i);
    end
end

