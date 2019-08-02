function h = generate_simple_muller_plots(U,M,cum_nodelbs,node_col)

%% how easy, or hard, is it to generate Muller plots?

% start with a simple tree and frequency matrices

% U = [ 1     1     1     1     1
%      0     1     0     0     1
%      0     0     1     0     0
%      0     0     1     1     0
%      0     0     0     0     1];
% 
% 
% M = [ 0.968384536881600   0.465248933683037   0.080651055801976   0.000000017366520   0.000000001857589
%    0.000003623239530   0.211425199822591   0.029129860186950   0.021889807919922   0.000000000935296
%    0.031609196671340   0.105263158424963   0.037422039296875   0.192468616817904   0.039682539395339
%    0.000000000332238   0.215389015850170   0.837577959234908   0.025773140117940   0.041624816489589
%    0.000002642875286   0.002673692219226   0.015219085479284   0.759868417777693   0.918692641322175];

%node_col = {[0 0   0.50],[ 1 0 0],[ 0  0.33  0.50],[0    0   0.66],[0    0   0.50]};

% get sizes
num_clust = size(M,1);
num_samples = size(M,2);

% get the frequency matrix

F = U*M;

% get the ancestors matrix
T = eye(length(U)) - inv(U);

% compute the lines

cvx_begin
	
	
	variable Up(size(M))
	variable Lw(size(M))
	
	minimize (	 norm(Lw) + norm(Up) 	); %force spacing between plots

	
	0 <= Up <= 1
	0 <= Lw <= 1
	Up >= Lw
	Up - Lw == F;
	for t = 1:num_samples	
		for i = 1:num_clust
			child = sort(find(T(i,:)));
			for j_ix = 1:length(child)
				j = child(j_ix);
				Lw(j,t) >= Lw(i,t);
				Up(j,t) <= Up(i,t);
				if (j_ix < length(child))
					Up(child(j_ix),t) <= Lw(child(j_ix+1),t);
				end
			end
		end
	end
	
cvx_end

%figure;
% do the filling
order_of_visit = bfsearch(digraph(T),1); %
hold on;
for curve_ix = order_of_visit'
	x = 1:num_samples;
	h = shade(x,Lw(curve_ix,:),x,Up(curve_ix,:),'FillType',[1 2;2 1],'FillColor',node_col{curve_ix},'LineStyle','none','FillAlpha',1);
	r(curve_ix) = plot(NaN,NaN,'Color',node_col{curve_ix});
end
	legend(r, cum_nodelbs);
hold off;

