import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

name_list = ["E", "P", "S"]

def draw_static_graph(G, data, cluster_name, pos_param, name_list):
    print(G)
    plt.clf()
    dot_n = 1
    edges_with_weight = [(name_list[i], name_list[j], data[j,i])  for i in range(len(name_list)) for j in range(len(name_list)) if i!=j]
    G.add_weighted_edges_from(edges_with_weight)
    edgelist = list(G.edges(data=True))
   
    edge_weight = [edge[2].get('weight',1) for edge in edgelist]
    edge_colors = [(0.8, 0.1, 0.1, 0.8) if x > 0 else (0.1, 0.1, 0.8, 0.8) if x<0 else (0.1, 0.1, 0.8, 0) for x in edge_weight]  # 所有边都使用相同的颜色
   
    edge_width = [max([5*abs(x),0]) for x in edge_weight] 
  
    dot_n = len(list(G.nodes()))
    pos = dict()
    for n in range(dot_n):
        node = list(G.nodes())[n]
        pos[node] = pos_param[n]

    nx.draw_networkx_nodes(G, pos, node_color=[(0.8, 0.1, 0.1, 1) if x > 0 else (0.1, 0.1, 0.8, 0.8) for x in [data[i,i] for i in range(3)]], node_size=[50*abs(data[i,i]) for i in range(3)])

    nx.draw_networkx_edges(G, pos, edgelist=edgelist,edge_color=edge_colors, width=edge_width, alpha=0.5, arrows=True)
  
    nx.draw_networkx_labels(G, pos, font_size=28, font_family='sans-serif')
   
    plt.axis("off")
  
    plt.pause(0.0007)
   

  
    plt.savefig(f"{cluster_name}_sn_new.tiff")
    
    
G = nx.DiGraph()


draw_static_graph(G, np.array(yysdes.Total_weight_["64"])[-1000,:,:], f'microbial 63 interactions', pos_param, name_list)
