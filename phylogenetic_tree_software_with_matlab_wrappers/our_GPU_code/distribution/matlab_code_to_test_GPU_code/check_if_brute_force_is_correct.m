%% all data. This is Vaughn's data

all_data = [1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000,1.0000
    , 0, 0, 0, 0, 0, 0, 0, 0, 0,1.0000, 0
    , 0, 0,0.0500,0.2000,0.2600,0.1000,0.2600,0.6000,0.6600,0.8800,0.9700
    , 0,0.2800, 0, 0,0.7100,0.9100,0.7600,0.4600,0.4000,0.1300,0.0400
    , 0, 0, 0, 0, 0,0.2500,0.3900,0.2200,0.4000,0.1400,0.0300
    , 0, 0, 0, 0, 0, 0,0.0800,0.1200,0.1200,0.3400,0.4000
    , 0, 0, 0, 0,0.1500,0.3800,0.3400,0.2900,0.3300,0.1500,0.0300
    , 0, 0, 0, 0, 0, 0, 0,0.0800,0.0200,0.2400,0.3100
    ,0.2800,0.2800,0.2200,0.2000, 0, 0, 0, 0, 0, 0, 0
    ,0.2100,0.2700,0.1800, 0, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0, 0, 0,0.1900,0.2600
    , 0, 0, 0, 0, 0, 0, 0, 0,0.2400,0.1200,0.0300
    , 0, 0, 0, 0, 0, 0, 0, 0,0.0300,0.2100, 0
    , 0,0.0600,0.1000,0.1700, 0, 0, 0, 0, 0, 0, 0 
    , 0,0.0800,0.0600,0.1700, 0, 0, 0, 0, 0, 0, 0
    ,0.0900,0.1100,0.1600, 0, 0, 0, 0, 0, 0, 0, 0
    ,0.1100, 0,0.1600,0.1300,0.1600, 0, 0, 0, 0,0.1300,0.1100
    , 0, 0, 0, 0, 0, 0, 0,0.0500,0.1000,0.1500,0.1400
    , 0, 0, 0, 0, 0,0.1500,0.1000,0.0600, 0, 0, 0
    , 0, 0,0.1400, 0, 0,0.1500, 0, 0,0.1400,0.0700, 0
    , 0, 0, 0, 0,0.1400, 0,0.0300, 0, 0, 0, 0
    , 0, 0,0.1400, 0, 0, 0, 0, 0, 0, 0,0.0500
    , 0, 0, 0, 0, 0, 0, 0,0.0400,0.1000,0.1400,0.1200
    , 0, 0,0.0300,0.1300, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0,0.1300, 0, 0,0.0200
    ,0.1200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0, 0,0.0600,0.0600,0.1100
    , 0, 0, 0, 0, 0, 0, 0, 0,0.1100,0.0600,0.0600
    , 0, 0,0.1100, 0, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0,0.1100, 0, 0,0.0200
    , 0, 0,0.1000, 0, 0, 0, 0, 0, 0, 0, 0
    ,0.0500, 0,0.1000, 0, 0, 0, 0,0.0400, 0, 0, 0
    , 0, 0,0.0900,0.0300, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0,0.0300,0.0900, 0,0.0400,0.0200
    ,0.0800,0.0900,0.0700, 0, 0, 0, 0, 0, 0, 0, 0
    , 0,0.0900,0.0400,0.0500, 0,0.0600,0.0500, 0,0.0700, 0, 0
    , 0,0.0200,0.0200, 0,0.0800, 0,0.0300, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0, 0,0.0600,0.0800, 0
    ,0.0200,0.0800,0.0200, 0, 0, 0, 0, 0, 0, 0, 0
    , 0,0.0400, 0, 0,0.0800, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0,0.0800, 0, 0, 0, 0, 0, 0
    ,0.0700, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    , 0, 0,0.0700,0.0400, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0, 0,0.0600,0.0700,0.0400
    ,0.0400, 0, 0, 0,0.0700,0.0500,0.0600, 0,0.0500, 0, 0
    , 0, 0, 0, 0, 0, 0, 0, 0,0.0600,0.0600,0.0700
    , 0, 0,0.0600, 0, 0, 0, 0, 0, 0, 0,0.0100
    ,0.0600, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    , 0, 0, 0, 0, 0, 0, 0,0.0600, 0, 0, 0
    ,0.0300, 0, 0, 0, 0, 0,0.0400, 0,0.0600, 0, 0];


%% select data amount
T = 11;
n = 3;
F = all_data(1:n,1:T);
topk = 13;

%%
% device = 'cpu';
% costfunction = 'cost1';
% pathtoprogram = '/home/jose/Dropbox/NIH\ Bio/C\ code\ to\ index\ trees/bei_code_to_score_trees_F_minus_UM/c\ code\ and\ cuda\ code\ with\ command\ line\ caller/a.out';
% 
% [cost, treeix, adj] = wrapper_brute_force_CPU_GPU_MP_cost_1234(F, device, costfunction, pathtoprogram);

%% generate a tree
% we are going to assume that the root is 1 because that is what our brute
% force software assumes

costfunctionid = 4; 
all_tree_vals  = nan(n^(n-2),2);

parpool(25);
parfor ii = 1:n^(n-2)
    
    
    tree_ix = 1 + n*(ii-1); % we jump every 6 tree indices because we want the root to always be zero.
    

    [Adj , depth, balls_levels_config , groups , links] = index_all_rooted_and_labeled_trees(tree_ix , n);
    root = find(diag(   Adj    ));
%     if (root ~= 1)
%         disp('Error');
%         return;
%     end
    Adj = Adj - diag(diag(Adj));
    % get the adj list from the adj mat of the tree
    AdjL = {};
    for i =1:n
        AdjL{i} = find(Adj(:,i));
    end
    % get adj list of tree
    Tree = BFS(AdjL,root);
    % get adj mat of tree
    AdjT = zeros(n);
    for i =1:n
        for j = Tree{i}
            AdjT(i,j) = 1;
        end
    end
    
    disp(AdjT);
    
    U = inv(eye(n) - AdjT);
    
  
    val = compute_tree_cost_in_cvx(F, U, n, T,costfunctionid);
        
    
    disp( val*val );
    
    %%
    all_tree_vals(ii,:) =  [val*val, tree_ix]; % square this value because the C code produces an objective squared

    
end
delete(gcp('nocreate'));
%% finally we extract the tree adj for the optimal tree
[Cost , b] = min(all_tree_vals(:,1));
[Adj , ~, ~ , ~ , ~] = index_all_rooted_and_labeled_trees(1 + n*(b-1) , n);
Adj = Adj - diag(diag(Adj));
r = sort(all_tree_vals,1);
for i = 1:topk
disp(r(i,1));
[Adj , ~, ~ , ~ , ~] = index_all_rooted_and_labeled_trees(1 + n*(b-1) , n);
Adj = Adj - diag(diag(Adj));
%disp(Adj);
end

%% This code checks if the brute force computation is correct

T = 11;
n = 6;
U = all_fs(b(1:n),:) - T*all_fs(b(1:n),:);

cost = 0;

for i = 1:T

    cvx_begin
        variable x(n)
        minimize norm(x-U(:,i))
        subject to
        sum(x)==1
        x >= 0
    cvx_end
    
    cost = cost + cvx_optval*cvx_optval;
    
end