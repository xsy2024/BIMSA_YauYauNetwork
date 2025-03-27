import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define a list of node names (e.g., representing three bacterial species)
name_list = ["E", "P", "S"]

def draw_graph(t, G, data, cluster_name, pos_param, name_list):
    """
    Draw the network graph for a given frame.
    
    Parameters:
        t           : Current time/frame index (used to index the data).
        G           : A NetworkX graph instance.
        data        : Data used to determine node sizes and colors (e.g., interaction strengths).
        cluster_name: A string for cluster naming (not directly used in drawing).
        pos_param   : List or dictionary of node positions.
        name_list   : List of node names.
    """
    print(G)
    plt.clf()  # Clear the current figure
    
    # The variable dot_n will hold the number of nodes in the graph
    dot_n = 1

    # Add weighted edges to the graph from a global variable (edges_with_weight should be defined elsewhere)
    G.add_weighted_edges_from(edges_with_weight)
    # Get list of edges along with their attributes
    edgelist = list(G.edges(data=True))
   
    # Extract the weight for each edge; if no weight is found, default to 1
    edge_weight = [edge[2].get('weight', 1) for edge in edgelist]
    # Define edge colors: if weight > 0 use a reddish color, if weight < 0 use a bluish color,
    # otherwise nearly transparent blue.
    edge_colors = [
        (0.8, 0.1, 0.1, 0.8) if x > 0 else (0.1, 0.1, 0.8, 0.8) if x < 0 else (0.1, 0.1, 0.8, 0)
        for x in edge_weight
    ]
    
    # Set edge width proportional to the absolute weight (scaled by 2.5)
    edge_width = [max([2.5 * abs(x), 0]) for x in edge_weight]
    
    # Get the actual number of nodes in the graph
    dot_n = len(list(G.nodes()))
    # Create a dictionary to store node positions using pos_param
    pos = dict()
    for n in range(dot_n):
        node = list(G.nodes())[n]
        pos[node] = pos_param[n]
   
    # Draw nodes with color and size based on diagonal values from 'data'
    # Assumes that for each node, the value is data[10*t, i, i] for i in range(3)
    nx.draw_networkx_nodes(
        G, pos,
        node_color=[
            (0.8, 0.1, 0.1, 1) if x > 0 else (0.1, 0.1, 0.8, 0.8)
            for x in [data[10*t, i, i] for i in range(3)]
        ],
        node_size=[50 * abs(data[10*t, i, i]) for i in range(3)]
    )

    # Draw edges with the computed colors and widths, set transparency to 0.5 and draw arrows
    nx.draw_networkx_edges(
        G, pos, edgelist=edgelist,
        edge_color=edge_colors,
        width=edge_width,
        alpha=0.5, arrows=True
    )
   
    # Draw node labels with a large font size and monospace font
    nx.draw_networkx_labels(G, pos, font_size=32, font_family="monospace")
  
    # Set the title of the plot
    plt.title('Video 2: Yau-Yau stochastic network among three bacterial species (E, P, S)')
    plt.pause(0.0007)  # Pause briefly to update the plot
    plt.axis('off')    # Turn off the axis
    # For the first frame, save the figure to a PNG file using the first edge width as filename
    if t == 0:
        plt.savefig(f"{edge_width[0]}.png")
   
# Create a directed graph object
G = nx.DiGraph()
# Define the number of nodes based on name_list
dot_num = len(name_list)
# Define the time steps (here from 0 to 399)
time_step = list(range(400))

# Create an animation using FuncAnimation. The draw_graph function is called for each frame.
ani = FuncAnimation(
    plt.gcf(),              # Get current figure
    draw_graph,             # Function to call for each frame
    frames=time_step,       # List of frame indices
    interval=50,            # Interval between frames in milliseconds
    fargs=(G, yysdes.Total_weight["64"], f'microbial 64 interactions', pos_param, name_list)
)

# Save the animation as an mp4 video file with specified FPS and DPI
ani.save("microbial_sample_ESP_3_63_new_yy.mp4", fps=20, dpi=500)
 
