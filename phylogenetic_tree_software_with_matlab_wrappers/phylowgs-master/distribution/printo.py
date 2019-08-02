from numpy	  import *
from numpy.random import *
from tssb		import *
from util		import *
from util2 import *

from ete2 import *

from subprocess import call
import sys

ctr=0

ancestry_Matrix = []
#ancestry_Matrix = [0] * 100

def print_top_trees(tree_archive,fout,k=5):
	global ctr;
	fout = open(fout,'w')
	fout_adj = open( 'top_k_trees_adjacency_matrix.txt', 'w' )
	fout_max_element = open( 'top_k_trees_max_element.txt', 'w' )
	fout_node_genes = open( 'top_k_trees_max_element_genes.txt', 'w' )
	fout_phi = open( 'top_k_trees_phi_values', 'w')
	tree_reader = TreeReader(tree_archive)

	for idx, (tidx, llh, tree) in enumerate(tree_reader.load_trees_and_metadata(k)):
			ctr=0
			del ancestry_Matrix[:]
			remove_empty_nodes(tree.root, None)
			# print top K trees in ascii
			fout_node_genes.write( "Nodes:\n" )
			print_best_tree(tree,fout,fout_node_genes,fout_phi)
			fout_node_genes.write( "\n" )
			
			#printing the ancestry_Matrix from the print_top_trees function
			print ('ancestry_Matrix for print_top_trees function: ')
			print( ancestry_Matrix ) 
			print( '\n' )
			#adding elements to the adjacency matrix
			#Finding the maximum element in the ancestry_Matrix array
			max_element = max( ancestry_Matrix )
			print ('Maximum node number = ')
			print( max_element ) 
			print( '\n' )

			adjacency_Matrix = numpy.zeros( (max_element + 1, max_element + 1), order='C' )
			for i in range( 0, len(ancestry_Matrix), 2 ):
				node_parent = ancestry_Matrix[i]
				node_child = ancestry_Matrix[i + 1]
				print('****node_parent' + str(node_parent) + '****node_child' + str(node_child) + '\n')
				adjacency_Matrix[node_parent][node_child] = 1
				adjacency_Matrix[node_child][node_parent] = 1			

			print( adjacency_Matrix )
			print( '\n' )

			for item in adjacency_Matrix:
				fout.write( "%s\n" % item )
			fout.write( "\n" )

			fout_adj.write( "Adjacency_matrix:\n" )	
			for item2 in adjacency_Matrix:
				fout_adj.write( "%s\n" % item2 )
			fout_adj.write( "\n\n" )	

			fout_max_element.write( "%d\n" % (max_element + 1) )


	tree_reader.close()
	fout.close()	
	fout_adj.close()
	fout_max_element.close()
	fout_node_genes.close()
	fout_phi.close()

def print_best_tree(tssb,fout,fout_node_genes,fout_phi):
	print ('Inside print_best_tree function!\n')
	nodes = tssb.get_nodes()
	nnodes = sum( [ 1 for node in nodes if len(node.get_data()) ] )
	
	#fout=open(fout,'w')
	t = Tree();t.name='0'
	fout.write('id, \t phi, \t nChildren, \t nGenes, \t genes \n')
	
	fout_phi.write('Phi_F_values\n')
	print_node2(tssb.root,None,t,fout,fout_node_genes,fout_phi)
	fout_phi.write('\n\n')
	fout.write('\n\n')
	fout.write(t.get_ascii(show_internal=True))
	fout.write('\n\n')	
	fout.write('Number of non-empty nodes in the tree: ' + repr(nnodes))	
	
	fout.write( '\nTransforming tree into adjacency matrix!\n' )
	#currnode = t.root;
	#fout.write( 'Example of a current node: ' )
	#fout.write( str(currnode) )


	fout.write('\n\n\n')	
	#currnode = t.root;
	#for 	
	#t.root['children']



	#fout.close()

def print_node2(node, parent,tree,fout,fout_node_genes,fout_phi):
	global ctr;
	global ancestry_Matrix
	print ('Inside print_node2 function!\n')

	num_data = node['node'].num_data()
	node_name  = ctr ; ctr+=1;
	fout.write( 'Current parent node is: ' + str(parent) + ' *******current counter is = ' + str(node_name) )
	fout.write( '\n' )
	parent_coord = parent
	child_coord = node_name
	print( parent_coord )
	print( child_coord )
	print( '\n' )

	#ancestry_Matrix = []
	if parent is not None:
		ancestry_Matrix.append(parent)
		ancestry_Matrix.append(node_name)
	print ('ancestry_Matrix for this recursive stage: ')
	print( ancestry_Matrix ) 
	print( '\n' )
	
	#if parent is not None:
	#	ancestry_Matrix[node_name] = 1

	genes = node['node'].get_data()
	gnames = ''
	if len(genes)>0:
		gnames = genes[0].id#name
		for g in arange(1,len(genes)):
			gnames = gnames + '; ' + genes[g].id#name;
	out_str = str(node_name) + ',\t' + str(around(node['node'].params,3)) + ',\t' + str(len(node['children'])) + ',\t' + str(len(genes)) + ',\t' + gnames +  '\n'
	genes_str = str(node_name) + '\t' + gnames + ']' + '\n'
	fout.write(out_str)
	phi_str = str(around(node['node'].params,3)) + '\n'
	fout_phi.write(phi_str)
	fout_node_genes.write(genes_str)

	for child in node['children']:
		name_string = str(ctr)#+'('+str(len(child['node'].get_data()))+')'
		fout.write( 'Example of a current ancestor to this node: ')
		#ancestors = []
		print( tree.get_ancestors()) 
		print( '\n' )
		fout.write( str(ctr) +'('+ str(child['node'].children()) + ')')
		fout.write( '\n' )
		print_node2(child, node_name,tree.add_child(name=name_string),fout,fout_node_genes,fout_phi)
	
### printing stuff #################
def print_best_tree_pdf(tssb,fout,score=0):
	#wts, nodes = tssb.get_mixture()
	#w = dict([(n[1], n[0]) for n in zip(wts,nodes)])
	print_tree_latex(tssb,fout,score)	


################ LATEX PRINTING ######################
global count
# writes code for tree
# root: root of the tree
# tree_file: string with latex code
def write_tree(root, tree_file):
	global count
	count+=1
	tree_file+='node {{{0}}}'.format(count)
	for child in root.children():
		tree_file+='child {'
		tree_file=write_tree(child, tree_file)
		tree_file+='}'
	return tree_file

# writes code for index
# root: root of the tree
# tree_file: string with latex code
def print_index(root, tree_file):
	global count
	count+=1
	tree_file+='{0} & '.format(count)
	ssms=''
	for datum in root.get_data():
		ssms+='{0}, '.format(datum.name)
	tree_file+=ssms.strip().strip(',')
	if root.get_data()==[]:
		tree_file+='-- '
	#else:
	#	tree_file=tree_file[:-1]
		#tree_file+=' & '
	tree_file+=' & '
	for i in range(len(root.params)):
		tree_file+='{0}, '.format(str(around(root.params[i],3)))
	tree_file=tree_file[:-2]
	tree_file+='\\\\\n'
	for child in root.children():
		tree_file=print_index(child, tree_file)
	return tree_file

# writes the latex code
# tssb: tssb structure of the tree
# fout: output file for latex
def print_tree_latex(tssb,fout,score):
	global count

	fout = open(fout,'w')
	count=-1
	#tree_file='\documentclass{article}\n'
	tree_file='\documentclass{standalone}\n'	
	tree_file+='\usepackage{tikz}\n'
	tree_file+='\usepackage{multicol}\n'
	tree_file+='\usetikzlibrary{fit,positioning}\n'
	tree_file+='\\begin{document}\n'
	tree_file+='\\begin{tikzpicture}\n'
	tree_file+='\\node (a) at (0,0){\n'
	tree_file+='\\begin{tikzpicture}\n'	
	tree_file+='[grow=east, ->, level distance=20mm,\
	every node/.style={circle, minimum size = 8mm, thick, draw =black,inner sep=2mm},\
	every label/.append style={shape=rectangle, yshift=-1mm},\
	level 2/.style={sibling distance=50mm},\
	level 3/.style={sibling distance=20mm},\
	level 4/.style={sibling distance=20mm},\
	every edge/.style={-latex, thick}]\n'
	tree_file+='\n\\'
	tree_file=write_tree(tssb.root['node'], tree_file)
	tree_file+=';\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='};\n'	
	count=-1
	tree_file+='\\node (b) at (a.south)[anchor=north,yshift=-.5cm]{\n'
	tree_file+='\\begin{tikzpicture}\n'	
	tree_file+='\\node (table){\n'
	tree_file+='\\begin{tabular}{|c|p{5cm}|p{5cm}|'
	for i in range(len(tssb.root['node'].params)):
		tree_file+='l|'
	tree_file+='}\n'
	tree_file+='\\hline\n'
	tree_file+='Node & \multicolumn{{1}}{{|c|}}{{Mutations}} & \multicolumn{{1}}{{|c|}}{{Frequencies}} \\\\\n'.format(len(tssb.root['node'].params))
	tree_file+='\\hline\n'
	tree_file=print_index(tssb.root['node'], tree_file)
	tree_file+='\\hline\n'
	tree_file+='\\end{tabular}\n'
	tree_file+='};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='};\n'	
	#tree_file+='\\node at (b.south) [anchor=north,yshift=-.5cm]{Posterior probability: ' + str(score) + '};\n'
	tree_file+='\\end{tikzpicture}\n'
	tree_file+='\end{document}\n'
	fout.write(tree_file)
	fout.close()	

