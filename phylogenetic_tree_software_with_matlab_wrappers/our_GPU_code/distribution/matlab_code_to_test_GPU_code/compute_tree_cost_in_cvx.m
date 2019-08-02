function [value , M] = compute_tree_cost_in_cvx(F, U, n, T,costfunctionid)

    warning('off','all');


    switch costfunctionid
            case 1
                cvx_begin quiet
                    variable M( n , T )
                    minimize ( norm(( F  - U*M),'fro') )
                    subject to
                        sum(M,1)==1
                        M >= 0
                cvx_end
            case 2
                cvx_begin quiet
                    variable M(n,T)
                    minimize norm((F - U*M),'fro')
                    subject to
                        sum(M,1) <= 1
                        M >= 0
                cvx_end
            case 3
                cvx_begin quiet
                    variable M(n,T)
                    minimize norm((inv(U)*F - M),'fro')
                    subject to
                        sum(M,1) == 1
                        M >= 0
                cvx_end
            case 4
                cvx_begin quiet
                    variable M(n,T)
                    minimize norm((inv(U)*F - M),'fro')
                    subject to
                        sum(M,1) <= 1
                        M >= 0
                cvx_end
    end
    
    value = cvx_optval;
        
end