% Copyright (c) 2019 Surjyendu Ray
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%% Generates Muller plots from U ancestry matrix, mutant frequency and
% related information
function [h,lgd] = generate_simple_muller_plots(U, M, cum_nodelbs, node_col)

% get sizes
num_clust = size(M,1);
num_samples = size(M,2);

% get the frequency matrix
F = U*M;

% get the ancestors matrix
T = eye(length(U)) - inv(U);

% compute the lines

cvx_begin quiet
	
	
	variable Up(size(M))
	variable Lw(size(M))
	variable slackUp(size(M,1),size(M,1),size(M,2))
	variable slackLw(size(M,1),size(M,1),size(M,2))
	variable slackNx(size(M))

	% a nice Muller plot just depends on the correct choice of optimization objective
	minimize (	 sum(sum(sum(slackUp.^2))) + sum(sum(sum(slackLw.^2))) + sum(sum(slackNx.^2))+  0*norm(Lw) + 0*norm(Up) + 0.*norm( Lw(:,2:end) - Lw(:,1:end-1)) + 0.*norm( Up(:,2:end) - Up(:,1:end-1))	);
	
	subject to
	
		slackUp >= 0
		slackLw >= 0
		slackNx >= 0

		0 <= Up <= 1
		0 <= Lw <= 1
		Up >= Lw
		Up - Lw == F;
		for t = 1:num_samples	
			for i = 1:num_clust
				child = sort(find(T(i,:)));
				for j_ix = 1:length(child)
					j = child(j_ix);
					Lw(j,t) == Lw(i,t) + slackLw(i,j,t);
					Up(j,t) + slackUp(i,j,t) == Up(i,t);
					if (j_ix < length(child))
						Up(child(j_ix),t) + slackNx(child(j_ix),t) == Lw(child(j_ix+1),t);
					end
				end
			end
		end
	
cvx_end

order_of_visit = bfsearch(digraph(T),1); %
hold on;
for curve_ix = order_of_visit'
	x = 1:num_samples;
	h = shade(x,Lw(curve_ix,:),x,Up(curve_ix,:),'FillType',[1 2;2 1],'FillColor',node_col{curve_ix},'LineStyle','none','FillAlpha',1);
	r(curve_ix) = plot(NaN,NaN,'Color',node_col{curve_ix});
end
	lgd = legend(r, cum_nodelbs);
hold off;




