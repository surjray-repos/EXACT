load('/home/surjray/Phylogeny_repo/phylogenetic_tree_software_with_matlab_wrappers/all_tests/all_results/our_code_on_real_data_tree_sizes_6_to_10_top_100trees.mat','all_our_code_outputs');

%%
c = newline;
clust_names = {['BCL2L13;',c],['COL24A1;',c],['DAZAP1;',c],['EXOC6B;',c],['GHDC;',c],['GPR158;',c],['HMCN1;',c],['KLHDC2;',c],['LRRC16A;',c],['MAP2K1;',c],['NAMPTL;',c],['NOD1;',c],['OCA2;',c],['PLA2G16;',c],['SAMHD1;',c],['SLC12A1;',c]};	
%clust_names = {'BCL2L13;','COL24A1;','DAZAP1;','EXOC6B;','GHDC;','GPR158;','HMCN1;','KLHDC2;','LRRC16A;','MAP2K1;','NAMPTL;',
               %'NOD1;','OCA2;','PLA2G16;','SAMHD1;','SLC12A1;'};	

%%
fileId = 27;
tree_size = 7;

all_costs = [];
for i = 1:100
	all_costs = [all_costs, all_our_code_outputs{fileId}{3 + tree_size}{i}{1}];
end

[~, sorted_order] = sort(all_costs);

%%
for i = 1:4
	T = all_our_code_outputs{fileId}{3 + tree_size}{sorted_order(i)}{3};
	clustIDs = all_our_code_outputs{fileId}{3 + tree_size}{sorted_order(i)}{6};
	s = cell(tree_size + 1, 1);
	for j = 1:tree_size + 1
		s{j} = cell2mat(clust_names(find(clustIDs == j)));
	end
	s{1} = '';

	figure;
	h = plot(graph(T + T'));
	labelnode(h, [1:tree_size + 1], s')
end
%%

all_probs = exp(-all_costs./(2*0.03*0.03));
all_probs = all_probs./(sum(all_probs));

prob_of_connection = zeros(16);
for i = 1:100
	for r = 1:16
		for l = 1:16
			if (r ~= l)
				node1 = all_our_code_outputs{fileId}{tree_size + 3}{i}{6}(r);	
				node2 = all_our_code_outputs{fileId}{tree_size + 3}{i}{6}(l);

				prob_of_connection(r,l) = prob_of_connection(r,l) + ( node1 == node2 || 1 == all_our_code_outputs{fileId}{tree_size + 3}{i}{3}(node1, node2)) * all_probs(i);
			end
		end
	end
	
end
%%

%plotting the probability of connection
prob_of_connection = prob_of_connection.*(prob_of_connection>0.013);

G = digraph(prob_of_connection);

Loc=rand(length(prob_of_connection),2);
weight=nonzeros(prob_of_connection);
figure;
h=plot(G);
labelnode(h,1:length(clust_names), clust_names)
%h.NodeFontSize = 14;
colormap jet           % select color palette 
h.EdgeCData=weight;    % define edge colors
%h.MarkerSize=10*rand(1,length(prob_of_connection)); % define node size
%h.XData=Loc(:,1);      % place node locations on plot
%h.YData=Loc(:,2);
colorbar
% % hide axes:
% set(gca,'XTickLabel',{' '})
% set(gca,'YTickLabel',{' '})
% set(gca,'YTick',[])
% set(gca,'XTick',[])

