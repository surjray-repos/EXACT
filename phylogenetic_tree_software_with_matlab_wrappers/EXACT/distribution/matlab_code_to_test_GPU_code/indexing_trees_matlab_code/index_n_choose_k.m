% this function indexes all combinations of n elements choose k
% there are nchoosek(n,k) many possible options for ix
% the ix variable can take values 1, 2, 3, 4, 5, ...
% maybe the function   index_n_choose_k    should go on a table to make the program faster

function [comb] = index_n_choose_k(ix , n , k)
    comb = nan(1,k); %vector
    eff_ix = ix;
    lb = 1;
    for l = 1:k
        r = 0;
        for a = lb : n - k + l 
            r = r + nchoosek( n - a , k - l );
            if(r >= eff_ix)
                break;
            end
        end
        eff_ix = eff_ix - r + nchoosek( n - a , k - l );
        comb(l) = a;
        lb = a + 1;
    end
end