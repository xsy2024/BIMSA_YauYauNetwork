import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define a list of node names (e.g., representing three bacterial species)
name_list = ["E", "P", "S"]

def draw_static_graph(G, data, cluster_name, pos_param, name_list):
    """
    Draws and saves a static directed graph with weighted edges.

    Parameters:
        G            : A NetworkX graph instance.
        data         : A 2D NumPy array containing weight values between nodes.
        cluster_name : A string used for the output filename.
        pos_param    : A list or dictionary of positions for each node.
        name_list    : A list of node names.
    """
    print(G)
    plt.clf()  # Clear the current figure

    dot_n = 1  # Initialize dot_n (will be overwritten below)

    # Generate weighted edges between nodes.
    # Create an edge from node i to node j with weight data[j,i] if i != j.
    edges_with_weight = [
        (name_list[i], name_list[j], data[j, i])
        for i in range(len(name_list))
        for j in range(len(name_list))
        if i != j
    ]
    
    # Add weighted edges to graph G.
    G.add_weighted_edges_from(edges_with_weight)
    edgelist = list(G.edges(data=True))  # Get list of edges along with their attribute dictionaries

    # Extract weights from the edge data; default weight is 1 if not specified.
    edge_weight = [edge[2].get('weight', 1) for edge in edgelist]
    
    # Define edge colors:
    # Use a reddish color for positive weights, bluish for negative weights, 
    # and nearly transparent blue for zero weight.
    edge_colors = [
        (0.8, 0.1, 0.1, 0.8) if x > 0 else (0.1, 0.1, 0.8, 0.8) if x < 0 else (0.1, 0.1, 0.8, 0)
        for x in edge_weight
    ]
    
    # Set edge widths proportional to the absolute weight (scaled by a factor of 5)
    edge_width = [max([5 * abs(x), 0]) for x in edge_weight]
    
    # Update dot_n to the number of nodes in graph G.
    dot_n = len(list(G.nodes()))
    pos = dict()  # Initialize node positions dictionary
    for n in range(dot_n):
        node = list(G.nodes())[n]
        pos[node] = pos_param[n]  # Assign positions from pos_param

    # Draw nodes:
    # - The node color is determined by the diagonal value data[i, i] for i in range(3)
    # - If the value is positive, use a reddish color; if negative, use a bluish color.
    # - Node size is scaled by 50 times the absolute value of data[i, i].
    nx.draw_networkx_nodes(
        G, pos,
        node_color=[
            (0.8, 0.1, 0.1, 1) if x > 0 else (0.1, 0.1, 0.8, 0.8)
            for x in [data[i, i] for i in range(3)]
        ],
        node_size=[50 * abs(data[i, i]) for i in range(3)]
    )

    # Draw edges with the calculated edge colors and widths.
    # Set transparency (alpha=0.5) and draw arrows on directed edges.
    nx.draw_networkx_edges(
        G, pos,
        edgelist=edgelist,
        edge_color=edge_colors,
        width=edge_width,
        alpha=0.5,
        arrows=True
    )

    # Draw labels for each node with a font size of 28 and a sans-serif font.
    nx.draw_networkx_labels(G, pos, font_size=28, font_family='sans-serif')
   
    plt.axis("off")  # Turn off the axis
    plt.pause(0.0007)  # Pause briefly to update the plot
  
    # Save the static graph as a TIFF image with a filename based on cluster_name.
    plt.savefig(f"{cluster_name}_sn_new.tiff")
    
# Create a directed graph instance.
G = nx.DiGraph()

# Call draw_static_graph with:
# - Graph G.
# - Data from yysdes.Total_weight_["64"] converted to a NumPy array and taking the last 1000 rows.
# - A cluster name string.
# - pos_param: positions for nodes (assumed to be defined elsewhere).
# - name_list: the list of node names.
draw_static_graph(G, np.array(yysdes.Total_weight_["64"])[-1000, :, :], f'microbial 63 interactions', pos_param, name_list)

