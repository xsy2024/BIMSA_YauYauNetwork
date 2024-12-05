import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

name_list = ["E", "P", "S"]

def draw_graph(t,G, data, cluster_name, pos_param, name_list):
    print(G)
    plt.clf()
   
    dot_n = 1

    G.add_weighted_edges_from(edges_with_weight)
    edgelist = list(G.edges(data=True))
   
    edge_weight = [edge[2].get('weight',1) for edge in edgelist]
    edge_colors = [(0.8, 0.1, 0.1, 0.8) if x > 0 else (0.1, 0.1, 0.8, 0.8) if x<0 else (0.1, 0.1, 0.8, 0) for x in edge_weight]  # 所有边都使用相同的颜色
    
    edge_width = [max([2.5*abs(x),0]) for x in edge_weight] 
    
    dot_n = len(list(G.nodes()))
    pos = dict()
    for n in range(dot_n):
        node = list(G.nodes())[n]
        pos[node] = pos_param[n]
   
    nx.draw_networkx_nodes(G, pos, node_color=[(0.8, 0.1, 0.1, 1) if x > 0 else (0.1, 0.1, 0.8, 0.8) for x in [data[10*t,i,i] for i in range(3)]], node_size=[50*abs(data[10*t,i,i]) for i in range(3)])

    nx.draw_networkx_edges(G, pos, edgelist=edgelist,edge_color=edge_colors, width=edge_width, alpha=0.5, arrows=True)
   
    nx.draw_networkx_labels(G, pos, font_size=32, font_family="monospace")
  
    plt.title('Video 2: Yau-Yau stochastic network among three bacterial species (E, P, S)')
    plt.pause(0.0007)
    plt.axis('off')  # 关闭坐标轴
    if t == 0:
        plt.savefig(f"{edge_width[0]}.png")
   
    
G = nx.DiGraph()
dot_num = len(name_list)
time_step = list(range(400))
ani = FuncAnimation(plt.gcf(), draw_graph, frames=time_step, interval = 50,fargs=(G,yysdes.Total_weight["64"],f'microbial 64 interactions',pos_param,name_list))
ani.save("microbial_sample_ESP_3_63_new_yy.mp4",fps=20,dpi=500)    
