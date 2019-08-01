% this function indexes all ways to assign n labels to k elements. Lables can be repeated
% there are n^k many possible options for ix
% the ix variable can take values 1, 2, 3, 4, 5, ...

function [labels] = index_n_labels_k_elements(ix , n , k)
    
    labels = nan(1,k); %this is a vector of size k
    ix = ix - 1;
    
    for i = 1:k
        labels(i) = 1 + mod(ix , n);
        ix = floor( ix  / n); %you do not need the floor function here
    end

end