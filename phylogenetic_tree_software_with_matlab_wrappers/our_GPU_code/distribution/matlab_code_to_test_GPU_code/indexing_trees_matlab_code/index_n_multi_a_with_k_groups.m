% this function indexes all ways to break n elements into groups of size a(1), a(2), ..., a(k)
% there are n choose (a(1) , a(2) , ..., a(k) ) ix
% the ix variable can take values 1, 2, 3, 4, 5, ...
% the output groups will be concatenated. we now how to read them based on the size of a

% a is a vector 

function [groups] = index_n_multi_a_with_k_groups(ix , n , a)
    
    groups = nan(1 , sum(a) ); %this is a vector of size equal to the sum of all elements in a

    free_elem = 1:n;  % this is a vector with all elements from 1 to n. Maybe it can be hard coded into   int free_elem[13] = {1, 2, 3, …, 13};

    k = length(a);

    ix = ix  - 1;

    r = 0;
    for i = 1:k
    
        sub_ix = 1 + mod( ix  , nchoosek( n-r   ,  a(i) ) );
        
        [comb] = index_n_choose_k(sub_ix , n-r , a(i)  );
        
        groups(1+r: r + a(i) ) = free_elem(comb);
        
        free_elem(comb) = []; %this part, when implemented in C, needs to be done differently. Right now is it is not very efficient. 
        
        ix = floor(ix / nchoosek( n-r   ,  a(i) ) );
        r = r + a(i);
    end



end