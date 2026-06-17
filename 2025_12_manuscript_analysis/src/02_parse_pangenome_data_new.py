import igraph as ig
from tqdm import tqdm
import pandas as pd
import csv
import plotly.graph_objects as go
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import argparse

# ===========================================================================================================
# Chapter 1 : Input
# ===========================================================================================================

parser = argparse.ArgumentParser(description="Pangenome graph analysis and plots")
parser.add_argument(
    "--input_path",
    required=True,
help = "Path to input directory containing the pangenome graph files"
)
parser.add_argument(
    "--output_path",
    required=True,
    help="Path to output directory where results will be written"
)
args = parser.parse_args()

# input
input_path = args.input_path
output_path = args.output_path

# make sure output directory exists
os.makedirs(output_path, exist_ok=True)

degree_output_file_path = os.path.join(input_path, 'degree_results.csv')
panaroo_file_path = os.path.join(input_path, 'panaroo_final_graph.gml')
PPanGGOLiN_file_path = os.path.join(input_path, 'pangenomeGraph.gexf')
ggcaller_file_path = os.path.join(input_path, 'ggcaller_final_graph.gml')
bifrost_file_path = os.path.join(input_path, 'bifrost.gfa')
cuttlefish_file_path = os.path.join(input_path, 'cuttlefish_2.gfa1')
pangraph_file_path = os.path.join(input_path, 'pangraph.gfa')

# ===========================================================================================================
# Chapter 2 : Import pangenome graph from file to igraph and query the number of nodes, edges, and subgraphs
# ===========================================================================================================

# Panaroo
panaroo_graph = ig.load(panaroo_file_path, format='gml')

# Query the number of nodes and edges
panaroo_num_nodes = panaroo_graph.vcount()
panaroo_num_edges = panaroo_graph.ecount()
panaroo_num_subgraphs = panaroo_graph.decompose()
print("Finish query nodes and edges panaroo")


# PPanGGOLiN
def PPanGGOLiN_to_igraph(input_file):
    # Creating an igraph
    graph = ig.Graph()

    # Used to keep track of edges that have been added to avoid duplicate additions
    added_edges = set()

    # Used to keep track of nodes that have been added to avoid duplicate additions
    added_nodes = set()
    total_lines = sum(1 for line in open(input_file, 'r'))

    with open(input_file, 'r') as file:
        for line in tqdm(file, desc="Processing PPanGGOLiN nodes and edges", unit=" lines", total=total_lines):
            # Process nodes
            if '<node' in line and 'id=' in line:
                node_id = line.split('id="')[1].split('"')[0]
                node_label = line.split('label="')[1].split('"')[0]
                graph.add_vertex(name=node_id, label=node_label)
                added_nodes.add(node_id)

            # Process edges
            elif '<edge' in line and 'id=' in line:
                edge_id = line.split('id="')[1].split('"')[0]
                source = line.split('source="')[1].split('"')[0]
                target = line.split('target="')[1].split('"')[0]

                # Ensure source and target nodes exist before adding the edge
                if source in added_nodes and target in added_nodes:
                    edge_name = f"{edge_id}"  # You can customize the edge name as needed
                    if edge_name not in added_edges:
                        graph.add_edge(source, target, name=edge_name)
                        added_edges.add(edge_name)

    return graph


PPanGGOLiN_output_igraph = PPanGGOLiN_to_igraph(PPanGGOLiN_file_path)

# ggcaller
ggcaller_graph = ig.load(ggcaller_file_path, format='gml')

# Query the number of nodes and edges
ggcaller_num_nodes = ggcaller_graph.vcount()
ggcaller_num_edges = ggcaller_graph.ecount()
ggcaller_num_subgraphs = ggcaller_graph.decompose()
print("Finish query nodes and edges ggcaller")


# bifrost
def bifrost_to_igraph(gfa_file):
    # Creating a igraph
    graph = ig.Graph()

    # Used to keep track of edges that have been added to avoid duplicate additions
    added_edges = set()

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Adding a progress bar with tqdm
        for line in tqdm(file, desc="Processing bifrost nodes and edges", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                _, node_id, sequence = line.strip().split('\t')
                graph.add_vertex(name=node_id, sequence=sequence)

            elif line.startswith('L'):  # Edge
                _, source, source_orient, target, target_orient, length = line.strip().split('\t')

                # Remove same but opposite edges
                forward_edge = (source, target, source_orient, target_orient)
                reverse_edge = (target, source, '-' if target_orient == '+' else '+',
                                '-' if source_orient == '+' else '+')

                if forward_edge not in added_edges and reverse_edge not in added_edges:
                    graph.add_edge(source, target, source_orient=source_orient, target_orient=target_orient)
                    added_edges.add(forward_edge)

    return graph


bifrost_output_igraph = bifrost_to_igraph(bifrost_file_path)


# cuttlefish
def cuttlefish_to_igraph(gfa_file):
    # Creating a igraph
    graph = ig.Graph()

    # Used to keep track of edges that have been added to avoid duplicate additions
    added_edges = set()

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Process lines starting with 'S'
        for line in tqdm(file, desc="Processing cuttlefish nodes", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                _, node_id, sequence, _ = line.strip().split('\t')
                graph.add_vertex(name=node_id, sequence=sequence)

        # Reset file pointer to the beginning of the file for processing 'L' lines
        file.seek(0)

        # Process lines starting with 'L'
        for line in tqdm(file, desc="Processing cuttlefish edges", unit=" lines", total=total_lines):
            if line.startswith('L'):  # Edge
                _, source, source_orient, target, target_orient, length = line.strip().split('\t')

                # Remove same but opposite edges
                forward_edge = (source, target, source_orient, target_orient)
                reverse_edge = (target, source, '-' if target_orient == '+' else '+',
                                '-' if source_orient == '+' else '+')

                if forward_edge not in added_edges and reverse_edge not in added_edges:
                    graph.add_edge(source, target, source_orient=source_orient, target_orient=target_orient)
                    added_edges.add(forward_edge)

    return graph


cuttlefish_output_igraph = cuttlefish_to_igraph(cuttlefish_file_path)


# pangraph
def pangraph_to_igraph(gfa_file):
    # Creating a igraph
    pangraph_graph = ig.Graph()

    # Used to keep track of edges that have been added to avoid duplicate additions
    added_edges = set()

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Adding a progress bar with tqdm
        for line in tqdm(file, desc="Processing Pangraph nodes and edges", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                _, node_id, _, _, _, = line.strip().split('\t')
                pangraph_graph.add_vertex(name=node_id)
            elif line.startswith('L'):  # Edge
                _, source, source_orient, target, target_orient, _, _, = line.strip().split('\t')

                # Remove same but opposite edges
                forward_edge = (source, target, source_orient, target_orient)
                reverse_edge = (target, source, '-' if target_orient == '+' else '+',
                                '-' if source_orient == '+' else '+')

                if forward_edge not in added_edges and reverse_edge not in added_edges:
                    pangraph_graph.add_edge(source, target, source_orient=source_orient, target_orient=target_orient)
                    added_edges.add(forward_edge)

    return pangraph_graph


pangraph_output_igraph = pangraph_to_igraph(pangraph_file_path)

# record nodes, edge, and subgraphs
results = []

# record nodes, edge, and subgraphs of panaroo
panaroo_result = ['panaroo', panaroo_num_nodes, panaroo_num_edges, len(panaroo_num_subgraphs)]
results.append(panaroo_result)

# record nodes, edge, and subgraphs of PPanGGOLiN
PPanGGOLiN_subgraphs = PPanGGOLiN_output_igraph.decompose()
PPanGGOLiN_result = ['PPanGGOLiN', PPanGGOLiN_output_igraph.vcount(), PPanGGOLiN_output_igraph.ecount(),
                     len(PPanGGOLiN_subgraphs)]
results.append(PPanGGOLiN_result)

# record nodes, edge, and subgraphs of ggcaller
ggcaller_result = ['ggcaller', ggcaller_num_nodes, ggcaller_num_edges, len(ggcaller_num_subgraphs)]
results.append(ggcaller_result)

# record nodes, edge, and subgraphs of bifrost
bifrost_subgraphs = bifrost_output_igraph.decompose()
bifrost_result = ['bifrost', bifrost_output_igraph.vcount(), bifrost_output_igraph.ecount(), len(bifrost_subgraphs)]
results.append(bifrost_result)

# record nodes, edge, and subgraphs of cuttlefish
cuttlefish_subgraphs = cuttlefish_output_igraph.decompose()
cuttlefish_result = ['cuttlefish', cuttlefish_output_igraph.vcount(), cuttlefish_output_igraph.ecount(),
                     len(cuttlefish_subgraphs)]
results.append(cuttlefish_result)

# record nodes, edge, and subgraphs of pangraph
pangraph_subgraphs = pangraph_output_igraph.decompose()
pangraph_result = ['pangraph', pangraph_output_igraph.vcount(), pangraph_output_igraph.ecount(),
                   len(pangraph_subgraphs)]
results.append(pangraph_result)

print("Finish query nodes and degrees")

# create DataFrame
columns = ['Method', 'Node Count', 'Edge Count', 'Subgraph Count']
df = pd.DataFrame(results, columns=columns)

# result
# nodes, edge, and subgraphs
output_file_path = os.path.join(output_path, 'Tool_basic_informations.tsv')
df.to_csv(output_file_path, sep='\t', index=False)

# ===========================================================================================================
# Chapter 3 : Create ridgeline plot for node degrees
# ===========================================================================================================

# create input format for ridgeline graph
# query degree
ggcaller_degrees = ggcaller_graph.degree()
panaroo_degrees = panaroo_graph.degree()
pangraph_degrees = pangraph_output_igraph.degree()
PPanGGOLiN_degree = PPanGGOLiN_output_igraph.degree()
bifrost_degree = bifrost_output_igraph.degree()
cuttlefish_degree = cuttlefish_output_igraph.degree()


# calculate degree count
def count_nodes_by_degree(degrees):
    degree_counts = {}
    for degree in degrees:
        if degree in degree_counts:
            degree_counts[degree] += 1
        else:
            degree_counts[degree] = 1
    return degree_counts


# get degree information
panaroo_degree_counts = count_nodes_by_degree(panaroo_degrees)
PPanGGOLiN_degree_counts = count_nodes_by_degree(PPanGGOLiN_degree)
ggcaller_degree_counts = count_nodes_by_degree(ggcaller_degrees)
bifrost_degree_counts = count_nodes_by_degree(bifrost_degree)
cuttlefish_degree_counts = count_nodes_by_degree(cuttlefish_degree)
pangraph_degree_counts = count_nodes_by_degree(pangraph_degrees)

panaroo_degree_counts = dict(sorted(panaroo_degree_counts.items()))
PPanGGOLiN_degree_counts = dict(sorted(PPanGGOLiN_degree_counts.items()))
ggcaller_degree_counts = dict(sorted(ggcaller_degree_counts.items()))
bifrost_degree_counts = dict(sorted(bifrost_degree_counts.items()))
cuttlefish_degree_counts = dict(sorted(cuttlefish_degree_counts.items()))
pangraph_degree_counts = dict(sorted(pangraph_degree_counts.items()))

# find max degree and mini degree
all_degrees = set(list(panaroo_degree_counts.keys()) + list(PPanGGOLiN_degree_counts.keys()) + list(
    ggcaller_degree_counts.keys()) + list(bifrost_degree_counts.keys()) + list(cuttlefish_degree_counts.keys()) + list(
    pangraph_degree_counts.keys()))
min_degree = min(all_degrees)
max_degree = max(all_degrees)

# create dataframe
degree_columns = ['Tool', 'Degree', 'Node_Count']
degree_results = []

# supplementi non-existent integer values
for tool, degree_counts in [('panaroo', panaroo_degree_counts), ('PPanGGOLiN', PPanGGOLiN_degree_counts),
                            ('ggcaller', ggcaller_degree_counts), ('bifrost', bifrost_degree_counts),
                            ('cuttlefish', cuttlefish_degree_counts), ('pangraph', pangraph_degree_counts)]:
    for degree in range(min_degree, max_degree + 1):
        count = degree_counts.get(degree, 0)
        degree_results.append([tool, degree, count])
degree_df = pd.DataFrame(degree_results, columns=degree_columns)

# record degree information to a csv file

print("Finish query degrees")
temp = degree_df

# Set tools list
Tool_list = ['panaroo', 'PPanGGOLiN', 'ggcaller', 'bifrost', 'cuttlefish', 'pangraph']

# set color
pal = sns.color_palette(palette='coolwarm', n_colors=6)

# calculate total node count for each tool
total_node_counts = {}
for Tool in Tool_list:
    total_node_counts[Tool] = temp[temp['Tool'] == Tool]['Node_Count'].sum()

# normalized degree and set x,y axis
array_dict = {}  # instantiating an empty dictionnary
for Tool in Tool_list:
    array_dict[f'x_{Tool}'] = temp[temp['Tool'] == Tool]['Degree']  # storing the temperature data for each Tool
    array_dict[f'y_{Tool}'] = temp[temp['Tool'] == Tool]['Node_Count'] / total_node_counts[
        Tool]  # storing the normalized node count for each Tool

# convert Seaborn colour schemes to Plotly colour format (RGB)
plotly_colors = [f'rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})' for r, g, b in pal]

fig = go.Figure()
# replace the 'line_color' parameter with the new Plotly colour scheme
for index, Tool in enumerate(Tool_list):
    fig.add_trace(go.Scatter(
        x=[0, max_degree], y=np.full(2, len(Tool_list) - 1 - index),
        mode='lines',
        line=dict(color=plotly_colors[index], width=1)  # use Plotly colour scheme
    ))

    fig.add_trace(go.Scatter(
        x=array_dict[f'x_{Tool}'],
        y=array_dict[f'y_{Tool}'] + (len(Tool_list) - 1 - index),
        fill='tonexty',
        name=f'{Tool}',
        line=dict(color=plotly_colors[index], width=0)  # use Plotly colour scheme
    ))

    # plotly.graph_objects' way of adding text to a figure
    fig.add_annotation(
        x=max_degree,
        y=len(Tool_list) - 1 - index,
        text=f'{Tool}',
        showarrow=False,
        font=dict(size=15, color=plotly_colors[index], family="Arial Black"),
        yshift=20
    )

    # Adding horizontal lines
    for i in range(1, 61):
        y_val = len(Tool_list) - 1 - index + i * 0.2
        if y_val <= 6:  # Ensure the line is within the y-axis range
            fig.add_trace(go.Scatter(
                x=[0, max_degree],
                y=[y_val, y_val],
                mode='lines',
                line=dict(color='grey', width=0.3, dash='dash')  # You can adjust the color and style as needed
            ))

fig.update_layout(
    showlegend=False,
    xaxis=dict(
        title=dict(
            text='Degree',  # set x axis title
            font=dict(
                size=15,  # set x axis font size
                family='Arial Black'
            )
        ),
        tickfont=dict(
            size=15,  # set y axis font size
            color='black',
            family='Arial Black'
        )
    ),
    yaxis=dict(
        title=dict(
            text='Proportion of Nodes',  # set y axis title
            font=dict(
                size=15,  # set y axis font size
                family='Arial Black'
            )
        ),
        tickfont=dict(
            size=15,  # set y axis font size
            color='black',
            family='Arial Black'
        ),
        range=[0, 6]  # set y axis range
    ),
    plot_bgcolor='white',
    paper_bgcolor='white'
)

# result
# ridgeline graph result
degree_output_filename = 'ridgeline_degree.png'
degree_output_path = os.path.join(output_path, degree_output_filename)

fig.write_image(degree_output_path)

# ============================================================================================================================
temp = temp[temp['Degree'] <= 10]
# Set tools list
Tool_list = ['panaroo', 'PPanGGOLiN', 'ggcaller', 'bifrost', 'cuttlefish', 'pangraph']

# set color
pal = sns.color_palette(palette='coolwarm', n_colors=6)

# calculate total node count for each tool
total_node_counts = {}
for Tool in Tool_list:
    total_node_counts[Tool] = temp[temp['Tool'] == Tool]['Node_Count'].sum()

# normalized degree and set x,y axis
array_dict = {}  # instantiating an empty dictionnary
for Tool in Tool_list:
    array_dict[f'x_{Tool}'] = temp[temp['Tool'] == Tool]['Degree']  # storing the temperature data for each Tool
    array_dict[f'y_{Tool}'] = temp[temp['Tool'] == Tool]['Node_Count'] / total_node_counts[
        Tool]  # storing the normalized node count for each Tool

# convert Seaborn colour schemes to Plotly colour format (RGB)
plotly_colors = [f'rgb({int(r * 255)}, {int(g * 255)}, {int(b * 255)})' for r, g, b in pal]

fig = go.Figure()
# replace the 'line_color' parameter with the new Plotly colour scheme
for index, Tool in enumerate(Tool_list):
    fig.add_trace(go.Scatter(
        x=[0, 10], y=np.full(2, len(Tool_list) - 1 - index),
        mode='lines',
        line=dict(color=plotly_colors[index], width=1)  # use Plotly colour scheme
    ))

    fig.add_trace(go.Scatter(
        x=array_dict[f'x_{Tool}'],
        y=array_dict[f'y_{Tool}'] + (len(Tool_list) - 1 - index),
        fill='tonexty',
        name=f'{Tool}',
        line=dict(color=plotly_colors[index], width=0)  # use Plotly colour scheme
    ))

    # plotly.graph_objects' way of adding text to a figure
    fig.add_annotation(
        x=10,
        y=len(Tool_list) - 1 - index,
        text=f'{Tool}',
        showarrow=False,
        font=dict(size=15, color=plotly_colors[index], family="Arial Black"),
        yshift=20
    )

    # Adding horizontal lines
    for i in range(1, 61):
        y_val = len(Tool_list) - 1 - index + i * 0.2
        if y_val <= 6:  # Ensure the line is within the y-axis range
            fig.add_trace(go.Scatter(
                x=[0, 10],
                y=[y_val, y_val],
                mode='lines',
                line=dict(color='grey', width=0.3, dash='dash')  # You can adjust the color and style as needed
            ))

fig.update_layout(
    showlegend=False,
    xaxis=dict(
        title=dict(
            text='Degree',  # set x axis title
            font=dict(
                size=15,  # set x axis font size
                family='Arial Black'
            )
        ),
        tickfont=dict(
            size=15,  # set y axis font size
            color='black',
            family='Arial Black'
        )
    ),
    yaxis=dict(
        title=dict(
            text='Proportion of Nodes',  # set y axis title
            font=dict(
                size=15,  # set y axis font size
                family='Arial Black'
            )
        ),
        tickfont=dict(
            size=15,  # set y axis font size
            color='black',
            family='Arial Black'
        ),
        range=[0, 6]  # set y axis range
    ),
    plot_bgcolor='white',
    paper_bgcolor='white'
)

# result
# ridgeline graph result
degree_lower_output_filename = 'ridgeline_degree_lower_10.png'
degree_output_path_10 = os.path.join(output_path, degree_lower_output_filename)

fig.write_image(degree_output_path_10)


# ============================================================================================================================

def get_nodes_and_degrees(graph):
    nodes = graph.vs["name"]
    degrees = graph.degree()
    data = {"Node": nodes, "Degree": degrees}
    df = pd.DataFrame(data)
    return df


# create a new input format for boxplot
panaroo_df = get_nodes_and_degrees(panaroo_graph)
PPanGGOLiN_df = get_nodes_and_degrees(PPanGGOLiN_output_igraph)
ggcaller_df = get_nodes_and_degrees(ggcaller_graph)
bifrost_df = get_nodes_and_degrees(bifrost_output_igraph)
cuttlefish_df = get_nodes_and_degrees(cuttlefish_output_igraph)
pangraph_df = get_nodes_and_degrees(pangraph_output_igraph)

# set tools
panaroo_df['Tool'] = 'panaroo'
PPanGGOLiN_df['Tool'] = 'PPanGGOLiN'
ggcaller_df['Tool'] = 'ggcaller'
bifrost_df['Tool'] = 'bifrost'
cuttlefish_df['Tool'] = 'cuttlefish'
pangraph_df['Tool'] = 'pangraph'

# merge dataframes
df = pd.concat([panaroo_df, PPanGGOLiN_df, ggcaller_df, bifrost_df, cuttlefish_df, pangraph_df])

# record degree information to a csv file
degree_output_file_path = os.path.join(output_path, 'degree_results.csv')
df.to_csv(degree_output_file_path, sep=',', index=False)

# input for querying sequence degree
df_degree_merge_input = pd.read_csv(degree_output_file_path)

# Calculate the median and IQR of each tool (length)
summary_stats_degree = df_degree_merge_input.groupby('Tool')['Degree'].agg(
    ['median', lambda x: np.percentile(x, 25), lambda x: np.percentile(x, 75)])
summary_stats_degree.columns = ['Median', 'Q1', 'Q3']
degree_summary_file_path = os.path.join(output_path, 'degree_summary_stats.csv')
# output to a graph
summary_stats_degree.to_csv(degree_summary_file_path)


# ===========================================================================================================
# Chapter 4 : Create beewarm plot and boxplot for node degrees
# ===========================================================================================================

# Functions to create beeswarm plot
def simple_beeswarm2_degree(y, nbins=None, width=1.):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """

    # Convert y to a numpy array to ensure it is compatible with numpy functions
    y = np.asarray(y)

    # If nbins is not provided, calculate a suitable number of bins based on data length
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)

    # Get the histogram of y and the corresponding bin edges
    nn, ybins = np.histogram(y, bins=nbins)

    # Find the maximum count in any bin to be used in calculating the x positions
    nmax = nn.max()

    # Create an array of zeros with the same length as y, to store x-coordinates
    x = np.zeros(len(y))

    # Divide indices of y-values into corresponding bins
    ibs = []
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        # Find the indices where y falls within the current bin
        i = np.nonzero((y > ymin) * (y <= ymax))[0]
        ibs.append(i)

    # Assign x-coordinates to the points in each bin
    dx = width / (nmax // 2)

    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            # Determine the starting index (j) based on the number of elements in the bin
            j = len(i) % 2

            # Sort the indices based on their corresponding y-values
            i = i[np.argsort(yy)]

            # Separate the indices into two halves (a and b) for arranging the points
            a = i[j::2]
            b = i[j + 1::2]

            # Assign x-coordinates to points in each half of the bin
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x


# set color beeswarm plot and boxplot
coolwarm_palette = sns.color_palette(palette='coolwarm', n_colors=len(df['Tool'].unique()))

# use setted color
fig, ax = plt.subplots(1, 1, figsize=(14, 12))
ax.get_xaxis().set_visible(False)  # Mask the x-axis

# Display each group:
count = 1
boxplot_data = []
for i, group in enumerate(df['Tool'].unique()):
    # Subset only observation from the group
    y = df[df['Tool'] == group]['Degree']

    # Get position of the observations
    x = simple_beeswarm2_degree(y, width=0.40)
    # Plot
    ax.plot(x + count, y, 'o', color=coolwarm_palette[i])

    count += 1  # Moves each group 1 unit to the right for the next iteration (avoid overlapping)
    boxplot_data.append(y)  # Add the values of the group to the `boxplot_data` variable

# Add the boxplots
ax.boxplot(boxplot_data,
           widths=0.5,  # Boxplots width
           )

# Add title and axis name
plt.title('Beeswarm and box plots of node degree', fontsize=16)
ax.set_ylabel("Degree", fontsize=16)

# Add a legend with the name of each distinct label
ax.legend(df['Tool'].unique(), fontsize=16)

# result
# beeswarm plot and boxplot graph
degree_boxplot_file_path = os.path.join(output_path, 'boxplot_degree.png')
plt.savefig(degree_boxplot_file_path)

# ===========================================================================================================
# Chapter 5 : Query node lengths
# ===========================================================================================================

# length
# query panaroo nodes sequence length
panaroo_lengths = panaroo_graph.vs["lengths"]
panaroo_node_ids = panaroo_graph.vs["name"]
panaroo_components = panaroo_graph.components()
panaroo_membership = panaroo_components.membership

# create dataframe with node_id and subgraph_id
panaroo_sequence_length = pd.DataFrame({
    "Tool": "panaroo",
    "length": panaroo_lengths,
    "node_id": panaroo_node_ids,
    "subgraph_id": panaroo_membership
})

print("Finish query lengths panaroo")


# PPanGGOLiN
def PPanGGOLiN_length(input_file):
    # Create a DataFrame to store the results
    result_df = pd.DataFrame(columns=['Tool', 'length', 'node_id'])

    total_lines = sum(1 for line in open(input_file, 'r'))

    current_node_id = None

    with open(input_file, 'r') as file:
        for line in tqdm(file, desc="Processing PPanGGOLiN sequence length", unit=" lines", total=total_lines):
            line = line.strip()

            # Track current node
            if '<node' in line and 'id="' in line:
                try:
                    current_node_id = line.split('id="')[1].split('"')[0]
                except IndexError:
                    current_node_id = None

            elif '</node' in line:
                current_node_id = None

            # Process length attribute within a node
            if '<attvalue' in line and 'for="9"' in line and current_node_id is not None:
                parts = line.split('"')
                # Expected format: <attvalue for="9" value="1234"/>
                if len(parts) >= 4:
                    node_length_str = parts[3]
                    try:
                        node_length = int(node_length_str)
                    except ValueError:
                        try:
                            node_length = float(node_length_str)
                        except ValueError:
                            continue

                    result_df = pd.concat(
                        [
                            result_df,
                            pd.DataFrame(
                                {
                                    'Tool': ['PPanGGOLiN'],
                                    'length': [node_length],
                                    'node_id': [current_node_id],
                                }
                            )
                        ],
                        ignore_index=True
                    )

    return result_df


# query PPanGGOLiN nodes sequence length and create dataframe
PPanGGOLiN_sequence_length = PPanGGOLiN_length(PPanGGOLiN_file_path)

# add subgraph_id using membership from PPanGGOLiN_output_igraph
PPanGGOLiN_components = PPanGGOLiN_output_igraph.components()
PPanGGOLiN_membership = PPanGGOLiN_components.membership
PPanGGOLiN_names = PPanGGOLiN_output_igraph.vs["name"]
PPanGGOLiN_map = {PPanGGOLiN_names[i]: PPanGGOLiN_membership[i] for i in range(len(PPanGGOLiN_names))}
PPanGGOLiN_sequence_length["subgraph_id"] = PPanGGOLiN_sequence_length["node_id"].map(PPanGGOLiN_map)

# query ggcaller nodes sequence length
ggcaller_lengths = ggcaller_graph.vs["lengths"]
ggcaller_node_ids = ggcaller_graph.vs["name"]
ggcaller_components = ggcaller_graph.components()
ggcaller_membership = ggcaller_components.membership

# create dataframe
ggcaller_sequence_length = pd.DataFrame({
    "Tool": "ggcaller",
    "length": ggcaller_lengths,
    "node_id": ggcaller_node_ids,
    "subgraph_id": ggcaller_membership
})

print("Finish query lengths ggcaller")


# bifrost
def bifrost_to_length(gfa_file):
    result_df = pd.DataFrame(columns=['Tool', 'length', 'node_id'])

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Adding a progress bar with tqdm
        for line in tqdm(file, desc="Processing bifrost sequence length", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                _, node_id, sequence = line.strip().split('\t')

                node_length = len(sequence)

                result_df = pd.concat(
                    [
                        result_df,
                        pd.DataFrame(
                            {
                                'Tool': ['bifrost'],
                                'length': [node_length],
                                'node_id': [node_id],
                            }
                        )
                    ],
                    ignore_index=True
                )

    return result_df


# query bifrost nodes sequence length and create dataframe
bifrost_sequence_length = bifrost_to_length(bifrost_file_path)

# add subgraph_id for bifrost
bifrost_components = bifrost_output_igraph.components()
bifrost_membership = bifrost_components.membership
bifrost_names = bifrost_output_igraph.vs["name"]
bifrost_map = {bifrost_names[i]: bifrost_membership[i] for i in range(len(bifrost_names))}
bifrost_sequence_length["subgraph_id"] = bifrost_sequence_length["node_id"].map(bifrost_map)


# cuttlefish
def cuttlefish_to_length(gfa_file):
    result_df = pd.DataFrame(columns=['Tool', 'length', 'node_id'])

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Adding a progress bar with tqdm
        for line in tqdm(file, desc="Processing cuttlefish sequence length", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                _, node_id, sequence, _, = line.strip().split('\t')

                node_length = len(sequence)

                result_df = pd.concat(
                    [
                        result_df,
                        pd.DataFrame(
                            {
                                'Tool': ['cuttlefish'],
                                'length': [node_length],
                                'node_id': [node_id],
                            }
                        )
                    ],
                    ignore_index=True
                )

    return result_df


# query cuttlefish nodes sequence length and create dataframe
cuttlefish_sequence_length = cuttlefish_to_length(cuttlefish_file_path)

# add subgraph_id for cuttlefish
cuttlefish_components = cuttlefish_output_igraph.components()
cuttlefish_membership = cuttlefish_components.membership
cuttlefish_names = cuttlefish_output_igraph.vs["name"]
cuttlefish_map = {cuttlefish_names[i]: cuttlefish_membership[i] for i in range(len(cuttlefish_names))}
cuttlefish_sequence_length["subgraph_id"] = cuttlefish_sequence_length["node_id"].map(cuttlefish_map)


# Pangraph
def pangraph_length(gfa_file):
    # Create a DataFrame to store the results
    result_df = pd.DataFrame(columns=['Tool', 'length', 'node_id'])

    # Get total number of lines in this file
    total_lines = sum(1 for line in open(gfa_file, 'r'))

    with open(gfa_file, 'r') as file:
        # Adding a progress bar with tqdm
        for line in tqdm(file, desc="Processing Pangraph sequence length", unit=" lines", total=total_lines):
            if line.startswith('S'):  # Node
                # Expected: S <node_id> * <length_field> ...
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                _, node_id, _, length_field, *_ = parts

                # Splitting the fourth part by ':'
                length_parts = length_field.split(':')

                # Check if there is a third part after splitting by ':'
                if len(length_parts) > 2:
                    length_number = length_parts[2]

                    try:
                        node_length = int(length_number)
                    except ValueError:
                        try:
                            node_length = float(length_number)
                        except ValueError:
                            continue

                    # Add the results to the DataFrame
                    result_df = pd.concat(
                        [
                            result_df,
                            pd.DataFrame(
                                {
                                    'Tool': ['pangraph'],
                                    'length': [node_length],
                                    'node_id': [node_id],
                                }
                            )
                        ],
                        ignore_index=True
                    )

    return result_df


# query pangraph nodes sequence length and create dataframe
pangraph_sequence_length = pangraph_length(pangraph_file_path)

# add subgraph_id for pangraph
pangraph_components = pangraph_output_igraph.components()
pangraph_membership = pangraph_components.membership
pangraph_names = pangraph_output_igraph.vs["name"]
pangraph_map = {pangraph_names[i]: pangraph_membership[i] for i in range(len(pangraph_names))}
pangraph_sequence_length["subgraph_id"] = pangraph_sequence_length["node_id"].map(pangraph_map)

# merge dataframes of 6 tools
df_length_merge = pd.concat([
    panaroo_sequence_length,
    PPanGGOLiN_sequence_length,
    ggcaller_sequence_length,
    bifrost_sequence_length,
    cuttlefish_sequence_length,
    pangraph_sequence_length
])

# output merged dataframe (Input file for querying sequence length)
output_file_path = os.path.join(output_path, 'length_results.csv')
df_length_merge.to_csv(output_file_path, sep=',', index=False)

print("Finish query lengths")

# input for querying sequence length
df_length_merge_input = pd.read_csv(output_file_path)

# ensure length is numeric
df_length_merge_input['length'] = pd.to_numeric(df_length_merge_input['length'], errors='coerce')

# Calculate the median and IQR of each tool (length)
summary_stats_length = df_length_merge_input.groupby('Tool')['length'].agg(
    ['median', lambda x: np.percentile(x.dropna(), 25), lambda x: np.percentile(x.dropna(), 75)])
summary_stats_length.columns = ['Median', 'Q1', 'Q3']

# output to a graph
length_summary_file_path = os.path.join(output_path, 'length_summary_stats.csv')
summary_stats_length.to_csv(length_summary_file_path)


# ===========================================================================================================
# Chapter 6 : Create beewarm plot and boxplot for node lengths
# ===========================================================================================================

# Functions to create beeswarm plot
def simple_beeswarm2_length(y, nbins=None, width=1):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """

    # Convert y to a numpy array to ensure it is compatible with numpy functions
    y = np.asarray(y)

    # If nbins is not provided, calculate a suitable number of bins based on data length
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)

    # Get the histogram of y and the corresponding bin edges
    nn, ybins = np.histogram(y, bins=nbins)

    # Find the maximum count in any bin to be used in calculating the x positions
    nmax = nn.max()

    # Create an array of zeros with the same length as y, to store x-coordinates
    x = np.zeros(len(y))

    # Divide indices of y-values into corresponding bins
    ibs = []
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        # Find the indices where y falls within the current bin
        i = np.nonzero((y > ymin) * (y <= ymax))[0]
        ibs.append(i)

    # Assign x-coordinates to the points in each bin
    dx = width / (nmax // 2)

    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            # Determine the starting index (j) based on the number of elements in the bin
            j = len(i) % 2

            # Sort the indices based on their corresponding y-values
            i = i[np.argsort(yy)]

            # Separate the indices into two halves (a and b) for arranging the points
            a = i[j::2]
            b = i[j + 1::2]

            # Assign x-coordinates to points in each half of the bin
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x


# set color scheme
coolwarm_palette = sns.color_palette(palette='coolwarm', n_colors=len(df_length_merge_input['Tool'].unique()))

# use coolwarm scheme
fig, ax = plt.subplots(1, 1, figsize=(14, 12))
ax.get_xaxis().set_visible(False)  # Mask the x-axis

# Display each group:
count = 1
boxplot_data = []
for i, group in enumerate(df_length_merge_input['Tool'].unique()):
    # Subset only observation from the group
    y = df_length_merge_input[df_length_merge_input['Tool'] == group]['length']

    # Get position of the observations
    x = simple_beeswarm2_length(y, width=0.40)
    y_log = np.log10(y)
    # Plot
    ax.plot(x + count, y_log, 'o', color=coolwarm_palette[i])

    count += 1  # Moves each group 1 unit to the right for the next iteration (avoid overlapping)
    boxplot_data.append(y_log)  # Add the values of the group to the `boxplot_data` variable

# Add the boxplots
ax.boxplot(boxplot_data,
           widths=0.5,  # Boxplots width
           )

# Add title and axis name
plt.title('Beeswarm and box plots of node length', fontsize=16)
ax.set_ylabel("Length (log10)", fontsize=16)

# Add a legend with the name of each distinct label
ax.legend(df_length_merge_input['Tool'].unique(), fontsize=16)

# Display the chart
length_boxplot_file_path = os.path.join(output_path, 'boxplot_length.png')
plt.savefig(length_boxplot_file_path)
plt.close()
# ===========================================================================================================
# Chapter 7 : Create histogram for the No. of nodes in each subgraph
# ===========================================================================================================

# calculate each node counts of each subgraph
panaroo_node_counts = [subgraph.vcount() for subgraph in panaroo_num_subgraphs]
PPanGGOLiN_node_counts = [subgraph.vcount() for subgraph in PPanGGOLiN_subgraphs]
ggcaller_node_counts = [subgraph.vcount() for subgraph in ggcaller_num_subgraphs]
bifrost_node_counts = [subgraph.vcount() for subgraph in bifrost_subgraphs]
cuttlefish_node_counts = [subgraph.vcount() for subgraph in cuttlefish_subgraphs]
pangraph_node_counts = [subgraph.vcount() for subgraph in pangraph_subgraphs]

panaroo_node_counts_df = pd.DataFrame(panaroo_node_counts, columns=['No. of nodes'])
panaroo_node_counts_df['Tool'] = 'panaroo'
PPanGGOLiN_node_counts_df = pd.DataFrame(PPanGGOLiN_node_counts, columns=['No. of nodes'])
PPanGGOLiN_node_counts_df['Tool'] = 'PPanGGOLiN'
ggcaller_node_counts_df = pd.DataFrame(ggcaller_node_counts, columns=['No. of nodes'])
ggcaller_node_counts_df['Tool'] = 'ggcaller'
bifrost_node_counts_df = pd.DataFrame(bifrost_node_counts, columns=['No. of nodes'])
bifrost_node_counts_df['Tool'] = 'bifrost'
cuttlefish_node_counts_df = pd.DataFrame(cuttlefish_node_counts, columns=['No. of nodes'])
cuttlefish_node_counts_df['Tool'] = 'cuttlefish'
pangraph_node_counts_df = pd.DataFrame(pangraph_node_counts, columns=['No. of nodes'])
pangraph_node_counts_df['Tool'] = 'pangraph'

subgraph_df = pd.concat(
    [panaroo_node_counts_df, PPanGGOLiN_node_counts_df, ggcaller_node_counts_df, bifrost_node_counts_df,
     cuttlefish_node_counts_df, pangraph_node_counts_df])

output_file_path = os.path.join(output_path, 'subgraph_results.csv')
subgraph_df.to_csv(output_file_path, sep=',', index=False)


# Functions to create beeswarm plot
def simple_beeswarm2_subgraph(y, nbins=None, width=1.):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    """

    # Convert y to a numpy array to ensure it is compatible with numpy functions
    y = np.asarray(y)

    # If nbins is not provided, calculate a suitable number of bins based on data length
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)

    # Get the histogram of y and the corresponding bin edges
    nn, ybins = np.histogram(y, bins=nbins)

    # Find the maximum count in any bin to be used in calculating the x positions
    nmax = nn.max()

    # Create an array of zeros with the same length as y, to store x-coordinates
    x = np.zeros(len(y))

    # Divide indices of y-values into corresponding bins
    ibs = []
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        # Find the indices where y falls within the current bin
        i = np.nonzero((y > ymin) * (y <= ymax))[0]
        ibs.append(i)

    # Assign x-coordinates to the points in each bin
    dx = width / (nmax // 2)

    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            # Determine the starting index (j) based on the number of elements in the bin
            j = len(i) % 2

            # Sort the indices based on their corresponding y-values
            i = i[np.argsort(yy)]

            # Separate the indices into two halves (a and b) for arranging the points
            a = i[j::2]
            b = i[j + 1::2]

            # Assign x-coordinates to points in each half of the bin
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx

    return x


# set color beeswarm plot and boxplot
coolwarm_palette = sns.color_palette(palette='coolwarm', n_colors=len(subgraph_df['Tool'].unique()))

# use setted color
fig, ax = plt.subplots(1, 1, figsize=(14, 12))
ax.get_xaxis().set_visible(False)  # Mask the x-axis

# Display each group:
count = 1
boxplot_data = []
for i, group in enumerate(subgraph_df['Tool'].unique()):
    # Subset only observation from the group
    y = subgraph_df[subgraph_df['Tool'] == group]['No. of nodes']

    # Get position of the observations
    x = simple_beeswarm2_subgraph(y, width=0.40)
    y_log = np.log10(y)

    # Plot
    ax.plot(x + count, y_log, 'o', color=coolwarm_palette[i])

    count += 1  # Moves each group 1 unit to the right for the next iteration (avoid overlapping)
    boxplot_data.append(y_log)  # Add the values of the group to the `boxplot_data` variable

# Add the boxplots
ax.boxplot(boxplot_data,
           widths=0.5,  # Boxplots width
           )

# Add title and axis name
# 设置 y 轴为对数刻度
plt.title('Beeswarm and box plots of node counts in each subgraph', fontsize=16)
ax.set_ylabel("No. of nodes (log10)", fontsize=16)

# Add a legend with the name of each distinct label
ax.legend(subgraph_df['Tool'].unique(), fontsize=16)

# result
# beeswarm plot and boxplot graph
subgraph_plot_output_file_path = os.path.join(output_path, 'beeswarm_plot_subgraph.png')
plt.savefig(subgraph_plot_output_file_path)
plt.close()
