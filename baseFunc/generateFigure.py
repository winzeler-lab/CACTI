# Write txt files
# Karla Godinez

# DECLARATIONS
import networkx as nx
from networkx.readwrite import cytoscape_data
import numpy as np

# import figure packages
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Create txt from dictionary
def similNet(prefix,path_output,net):
	#net = ['Source', 'Target', 'Simil']
	#nodes = ['Cmp','SMILES','Canonical smiles','Type']

	# Network visualization properties
	node_color = []
	edge_color = []
	
	# Draw network
	pos = nx.spring_layout(net)
	#edges,weights = zip(*nx.get_edge_attributes(net,'Simil').items())
	#nx.draw(net,pos, edgelist=edges, edge_color=weights, edge_cmap=plt.cm.Blues)
	nx.draw(net,pos,alpha=0.6)
	nx.draw_networkx_labels(net,pos)

	# Save figure
	#plt.show()
	plt.plot()
	plt.savefig(path_output+'/'+prefix+'_similarityNetwork.svg')
	plt.savefig(path_output+'/'+prefix+'_similarityNetwork.png')
	plt.close()

def gridSimilNet(prefix,path_output,net):
	pos = nx.spring_layout(net)
	
	# Get node shape
	node_db = [n for (n,ty) in nx.get_node_attributes(net,'Type').items() if ty == 'database']
	node_ds = [n for (n,ty) in nx.get_node_attributes(net,'Type').items() if ty == 'dataset'] 


	# Count node relationships
	val = []
	shape = []
	for n in net.nodes():
		val.append(len(net.edges([n])))
		net.node[n]['val'] = len(net.edges([n]))
		shape.append('o') if net.node[n]['Type'] == 'dataset' else shape.append('d')
		if net.node[n]['Type'] == 'dataset':
			net.node[n]['shape'] = 'o'
		else:
			net.node[n]['shape'] = 'd'
		
	
	# Plot network
	#nx.draw_networkx_nodes(net, pos, nodelist=node_ds, node_shape='o')
	#nx.draw_networkx_nodes(net, pos, nodelist=node_db, node_shape='d')
	nx.draw_networkx_edges(net,pos)
	nx.draw_networkx_labels(net,pos)


	nc_color = list(nx.get_node_attributes(net,'val').values())
	nc_shape = list(nx.get_node_attributes(net,'shape').values()) 
	nc = nx.draw_networkx_nodes(net,pos,node_color=nc_color, cmap=plt.cm.Blues)
	nc.set_edgecolor('black')

	plt.plot()
	plt.colorbar(nc)
	plt.savefig(path_output+'/'+prefix+'_similarityNetwork.svg')
	#plt.savefig(path_output+'/'+prefix+'_similarityNetwork.png')
	plt.close()

def simHeatmap(prefix,path_output,cmp_name,simmat):
	import matplotlib.cm as cm

	fig, ax = plt.subplots(figsize=(len(simmat)*2, len(simmat)*2))
	mask = np.zeros_like(simmat, dtype='bool')
	mask[np.triu_indices_from(simmat, k=1)] = True
	sns.heatmap(simmat, cmap=cm.viridis, mask=mask,linewidths=0.5,
				xticklabels=cmp_name, yticklabels=cmp_name, ax=ax)
	plt.savefig(path_output + '/' + prefix + '_similarityHeatmap.png')
	plt.close()

# Network cytoscape based
def similCytNet(prefix,path_output,net):
	from networkx.readwrite import cytoscape_data
	from cyjupyter import Cytoscape
	plt.plot()
	cy_g = cytoscape_data(net)
	cyobj=Cytoscape(data=cy_g, layout_name='cose')
	plt.savefig()
	#cyobj.savefig(path_output+'/'+prefix+'_fig')
	




