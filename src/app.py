import pandas as pd
import numpy as np
import matplotlib as mpl 
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import mean_squared_error
import plotly.graph_objs as go
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import os
import base64
import random


# Read data
ge_files = os.listdir('./ge_csv/')
ge_data_list = {}

# Read and melt the gene expression data matrix 
for _file in ge_files: 
    ge_data_list['_'.join(_file.split('GSE')[0].split('_')[0:3])] = pd.melt(pd.read_csv(os.path.join('./ge_csv', _file), 
                                                                                        index_col = 0), 
                                                            id_vars="index", var_name="gene", value_name="expression")
    

    
def continuous_colors(label_list, colormap='viridis', custom_color=None, orders=None):
    """
    Generate continuous colors for a list of labels.
    """
    color_dict = {}

    # Choose colormap
    if isinstance(colormap, str):
        cmap = colormaps.get_cmap(colormap)
    else:
        cmap = colormap

    # Generate custom colormap
    if custom_color is not None:
        custom_cmap = LinearSegmentedColormap.from_list('custom_cmap', custom_color)
    else:
        custom_cmap = None

    # Generate colors
    num_labels = len(label_list)
    for i, label in enumerate(label_list):
        if custom_cmap is not None:
            norm_color = i / (num_labels - 1)  # Normalize color index
            color = mpl.colors.rgb2hex(custom_cmap(norm_color))
        else:
            color = mpl.colors.rgb2hex(cmap(i / (num_labels - 1)))  # Normalize color index
        color_dict[label] = color

    # Reorder color_dict based on orders if provided
    if orders is not None:
        color_dict = {label: color_dict[label] for label in orders if label in color_dict}
    return color_dict

def distinct_colors(label_list, category=None, custom_color=None, random_state=0):
    """
    Generate distinct colors for a list of labels.

    Parameters:
    label_list (list): A list of labels for which you want to generate distinct colors.
    category (str): Category of distinct colors. Options are 'warm', 'floral', 'rainbow', or None for random. Default is None.

    Returns:
    dict: A dictionary where labels are keys and distinct colors (in hexadecimal format) are values.

    Example:
    >>> labels = ['A', 'B', 'C']
    >>> color_mapping = distinct_colors(labels, category='warm')
    >>> print(color_mapping)
    {'A': '#fabebe', 'B': '#ffd8b1', 'C': '#fffac8'}
    """
    random.seed(random_state)
    
    warm_colors = ['#fabebe', '#ffd8b1', '#fffac8', '#ffe119', '#ff7f00', '#e6194B']
    floral_colors = ['#bfef45', '#fabed4', '#aaffc3', '#ffd8b1', '#dcbeff', '#a9a9a9']
    rainbow_colors = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4']
    pastel_colors = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', 
                     '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928', 
                     '#8DD3C7', '#BEBADA', '#FFED6F']
    
    color_dict = {}

    if custom_color is not None: 
        assert len(custom_color) >= len(label_list), "Provided label_list needs to be shorter than provided custom_color"
        for i, _label in enumerate(label_list): 
            color_dict[_label] = custom_color[i]
        return color_dict

    color_palette = None
    if category == 'warm':
        color_palette = random.sample(warm_colors, len(warm_colors))
    elif category == 'floral':
        color_palette = random.sample(floral_colors, len(floral_colors))
    elif category == 'rainbow':
        color_palette = random.sample(rainbow_colors, len(rainbow_colors))
    elif category == 'pastel': 
        color_palette = random.sample(pastel_colors, len(pastel_colors))
    else:
        color_palette = random.sample(warm_colors + floral_colors + rainbow_colors + pastel_colors, len(label_list))
    
    for i, label in enumerate(label_list):
        color_dict[label] = color_palette[i % len(color_palette)]
    
    return color_dict

def find_nearest_key(dictionary, value):
    return min(dictionary.keys(), key=lambda x: abs(x - value))

def get_rmse_table(ge_data, subject_gene = 'Rtp1'):
    
    ge_table = pd.pivot_table(ge_data, 
                              index = 'index',
                              columns = 'gene', 
                              values = 'expression')    
    genes = [col for col in ge_table.columns if not any(exc_str in col for exc_str in ['dpt_average', '_sum'])]
    # Calculating mean_square_error for all the columns in a dictionary form 
    mse_values = (ge_table[genes].sub(ge_table[subject_gene], axis=0)**2).mean()
    return mse_values

def find_enriched_genes(ge_data, query_dpt, n_enriched_genes = 10): 
    
    ge_table = pd.pivot_table(ge_data, 
                              index = 'index',
                              columns = 'gene', 
                              values = 'expression')    
    
    genes = [col for col in ge_table.columns if not any(exc_str in col for exc_str in ['dpt_average', '_sum'])]
    
    dataframe = ge_table[genes].copy()
    
    # Extract the numerical data and convert it to a numpy array
    data_array = dataframe.values

    # Extract the index column and convert it to a numpy array
    index_array = dataframe.index.values

    # Find the indices of rows corresponding to the given index list
    selected_indices = np.where(np.isin(index_array, query_dpt))[0]

    # Calculate mean and standard deviation for each column
    column_mean = np.mean(data_array, axis=0)
    column_std = np.std(data_array, axis=0)

    # Calculate the threshold for increase
    threshold = 0.5 * column_std

    # Find columns where the values for the selected indexes are higher than the rest
    # is_enriched = (data_array > (column_mean + threshold))
    is_enriched = (data_array > column_mean)

    # Filter columns where maximum value corresponds to the enriched indexes
    filtered_column_indices = np.where(np.all(is_enriched[selected_indices], axis=0))[0]

    mean_diff = np.mean(data_array[selected_indices] - column_mean, axis=0)

    # Get the column names sorted by mean error (from most enriched to least enriched)
    sorted_indices = np.argsort(mean_diff[filtered_column_indices])[::-1]
    filtered_column_names = dataframe.columns[filtered_column_indices[sorted_indices]].tolist()


    # print('DEBUG 1: ', query_dpt)
    # print('DEBUG 2: ', filtered_column_names[0:n_enriched_genes])
    return filtered_column_names[0:n_enriched_genes] if n_enriched_genes > 0 else []

def expression_plot(subject_gene, query_genes, n_top_genes, n_bot_genes, enriched_celltype, n_enriched_genes,ge_data):
    """
    Generate the plot for gene expression trajectory.
    """
    
    if subject_gene is not None: 
        rmse_data = get_rmse_table(ge_data, subject_gene)
    else: 
        rmse_data = None 
    
    
    # Create traces for query genes
    traces = []
    colorscale = None 
    
    # Plot n of genes associated with subject gene 
    if (any(n_top_genes) or any(n_bot_genes)) and (rmse_data is not None): 
        # Set colorscale based on the top n of genes plotting 
        top_genes_to_plot = rmse_data.sort_values(ascending=True)[n_top_genes[0]+1:n_top_genes[1]+1].index.to_list()
        bot_genes_to_plot = rmse_data.sort_values(ascending=False)[n_bot_genes[0]:n_bot_genes[1]].index.to_list()

        genes_to_plot = top_genes_to_plot + bot_genes_to_plot
        
        genes_to_plot += query_genes
        rmse_data = rmse_data[rmse_data.index.isin(genes_to_plot)]
        colorscale = continuous_colors(list(np.round(np.linspace(min(rmse_data),
                                                                 max(rmse_data),
                                                                 100), 3)))
    
        for _gene in genes_to_plot:
            subset = ge_data[ge_data['gene'] == _gene]
            color_idx = find_nearest_key(colorscale, np.round(rmse_data[_gene], 3))
            
            # print('if pass xxxxxx', _gene, color_idx, colorscale.keys()) # for debugging

            trace = go.Scatter(x=subset['index'], y=subset['expression'],
                            mode='lines', line_shape='spline',
                            name=_gene,
                            line=dict(color=colorscale[color_idx], width=5),
                            hoverinfo='text',
                            hovertext=f"{_gene} <br> RMSE: {rmse_data[_gene]:.3f}")
            traces.append(trace)
        
    # Plot individually input genes from dropdown 
    elif query_genes: 
        # Plot genes in continuous color based on rmse refrence to subject_gene 
        if rmse_data is not None: 
            rmse_data = rmse_data[rmse_data.index.isin(query_genes)]
            colorscale = continuous_colors(list(np.round(np.linspace(min(rmse_data),
                                                                    max(rmse_data),
                                                                    100), 3)))
            for _gene in query_genes:
                query_data_q = ge_data[ge_data['gene'] == _gene]
                color_idx = find_nearest_key(colorscale, np.round(rmse_data[_gene], 3))
                
                # print('else yyyyyyyy',_gene, color_idx, colorscale.keys()) # for debugging 
                
                trace = go.Scatter(x=query_data_q['index'], y=query_data_q['expression'],
                                mode='lines', line_shape='spline',
                                name=_gene,
                                line=dict(color=colorscale[color_idx], width=5),
                                hoverinfo='text',
                                hovertext=f"{_gene} <br> RMSE: {rmse_data[_gene]:.3f}")
                traces.append(trace)
        # Plot genes in distinct color when there is no subject_gene to compare to
        else: 
            d_colors = distinct_colors(query_genes, category='pastel')
            for _gene in query_genes:
                query_data_q = ge_data[ge_data['gene'] == _gene]        
                trace = go.Scatter(x=query_data_q['index'], y=query_data_q['expression'],
                                mode='lines', line_shape='spline',
                                name=_gene,
                                line=dict(color=d_colors[_gene], width=5),
                                hoverinfo='text',
                                hovertext=f"{_gene}")
                traces.append(trace)
    
    # Plot top 10 genes enriched in defined celltype combination
    if np.any([_celltype in ['GBC', 'INP', 'early iOSN', 'late iOSN', 'mOSN'] for _celltype in enriched_celltype]): 
        cell_type_dpt = {'GBC': ['0-0.3'],
                         'INP': ['0.3-0.6'], 
                         'early iOSN': ['0.6-0.7', '0.7-0.8'], 
                         'late iOSN': ['0.8-0.9', '0.9-0.95'],
                         'mOSN': ['0.95-1']}
        query_dpt = [item for sublist in [cell_type_dpt[_celltype] for _celltype in enriched_celltype] for item in sublist]
        genes = find_enriched_genes(ge_data, query_dpt, n_enriched_genes)
        # print('DEBUG 2: ', query_dpt)
        if genes: 
            if (rmse_data is not None) & (colorscale is None): 
                rmse_data = rmse_data[rmse_data.index.isin(genes)]
                colorscale = continuous_colors(list(np.round(np.linspace(min(rmse_data),
                                                                        max(rmse_data),
                                                                        100), 3)))
                for _gene in genes:
                    subset = ge_data[ge_data['gene'] == _gene]
                    color_idx = find_nearest_key(colorscale, np.round(rmse_data[_gene], 3))
                    trace = go.Scatter(x=subset['index'], y=subset['expression'],
                                    mode='lines', line_shape='spline',
                                    name=_gene,
                                    line=dict(color=colorscale[color_idx], width=5),
                                    hoverinfo='text',
                                    hovertext=f"{_gene} <br> RMSE: {rmse_data[_gene]:.3f}")
                    traces.append(trace)
            else: 
                d_colors = distinct_colors(genes, category='pastel')
                for _gene in genes:
                    query_data_q = ge_data[ge_data['gene'] == _gene]        
                    trace = go.Scatter(x=query_data_q['index'], y=query_data_q['expression'],
                                    mode='lines', line_shape='spline',
                                    name=_gene,
                                    line=dict(color=d_colors[_gene], width=5),
                                    hoverinfo='text',
                                    hovertext=f"{_gene}")
                    traces.append(trace)
           
        
    # Create trace for subject gene
    if subject_gene is not None: 
        subject_data = ge_data[ge_data['gene'] == subject_gene]
        sub_trace = go.Scatter(x=subject_data['index'], y=subject_data['expression'],
                            mode='lines', line_shape='spline',
                            name=subject_gene,
                            line=dict(color='#F4D1FF', width=7),
                            hoverinfo='text',
                            hovertext=subject_gene)
        traces.append(sub_trace)
        
    # Define layout and create plot 
    layout = go.Layout(xaxis=dict(title="Pseudotime across OSN lineage", 
                                showticklabels=True, showgrid=False),
                    yaxis=dict(title="Normalized expression", 
                                showline=True, showgrid=False), 
                    template='simple_white')
    fig = go.Figure(data=traces, layout=layout)

    # # Create color bar based on colorscale
    if colorscale is not None: 
        c_min = np.round(list(colorscale.items())[0][0], 2)
        c_max = np.round(list(colorscale.items())[-1][0], 2)
        colorbar_trace = go.Scatter(x=[None], y=[None],
                                    showlegend=False, 
                                    mode='markers',
                                    hoverinfo='none', 
                                    marker=dict( 
                                        colorbar_x = -0.24,
                                        colorscale='Viridis', 
                                        showscale=True,
                                        cmin=c_min, cmax=c_max,
                                        colorbar=dict(thickness=10, tickvals=[c_min, c_max], ticktext=[str(c_min), str(c_max)])
                                        ),
                                    )
        fig.add_trace(colorbar_trace)
        fig.update_layout(annotations = [dict(xref="paper",
                                              yref="paper",
                                        text = 'Expression RMSE', 
                                        textangle=-90, 
                                        showarrow = False,
                                        x = -0.27, y = 0.5
                                        )])

    
    # Inserts x axis cell lineage bar. 
    image_filename = './img/geLine_cellbar.png'
    plotly_logo = base64.b64encode(open(image_filename, 'rb').read())
    fig.update_layout(images= [dict(
                        source='data:image/png;base64,{}'.format(plotly_logo.decode()),
                        xref="paper", yref="paper",
                        x=0, y=-0.2,
                        sizex=1, sizey=0.3,
                        xanchor="left",
                        yanchor="top",
                        sizing="stretch",
                        layer="above")], 
                        margin=dict(l=150, r=0, t=20, b=150)
                        )
    
    return fig

# Define Dash app
app = Dash(__name__, suppress_callback_exceptions=True)
server = app.server 

# Define dropdown options placeholders
genes = [{'label': gene, 'value': gene} for gene in ge_data_list[list(ge_data_list.keys())[0]]['gene'].unique()]
sub_gene_dropdown_options = genes
que_genes_dropdown_options = genes

# Define app layout
app.layout = html.Div(children=[
    html.Div([
        html.H1("Gene Expression trajectory in Olfactory Sensory Neurons")
    ], style={"display": "inline-block", "align-items": "center", 
              "margin-left": "50px", 'background-color': 'white', 'font-size': '16px'}),
    
    # Control panel 
    html.Div([
        html.Div([
            html.Label("Plot type: ", htmlFor="plot-type", 
                   style={'font-weight': 'bold', 'margin-bottom': '20px'}),
            dcc.Checklist(id='plot-type',
                                options=[{'label': 'Logmaritize counts', 'value': 'logmaritize'},
                                        {'label': 'Normalize maximum expression', 'value': 'norm'}],
                                value=['norm'])
        ], style={'margin-bottom': '30px'
                #   'border': '2px grey solid'
                  }),

        html.Div([
            html.Label("Input gene: ", style={'font-weight': 'bold', 'margin-bottom': '20px'}),
            html.Label('Subject gene', style={"text-align": "center", 
                                              'display': 'block', 
                                              'padding': '0px 0px 10px 10px'}),
            dcc.Dropdown(id='sub_gene_dropdown',
                           options=sub_gene_dropdown_options,
                           value='Omp',
                           multi=False),
            html.Label('Query genes', style={"text-align": "center", 
                                             'display': 'block', 
                                             'padding': '0px 0px 10px 10px'}),
            dcc.Dropdown(id='que_genes_dropdown',
                           options=que_genes_dropdown_options,
                        #    value=[],
                           value=['Hmgb2', 'Tubb5', 'Gap43','Rtp1'], 
                           multi=True),
            
            html.Div([
                html.Label("Most associated genes: ", htmlFor="top_associated_genes"),
                dcc.RangeSlider(0,30,
                                step=1,
                                id='top_associated_genes',
                                value=[0,0],
                                marks={str(num_gene): str(num_gene) for num_gene in range(0,31,10)}), 
                html.Label("Least associated genes: ", htmlFor="bottom_associated_genes", 
                                 style={'margin-bottom': '20px'}),
                dcc.RangeSlider(0,30,
                                step=1,
                                id='bottom_associated_genes',
                                value=[0,0],
                                marks={str(num_gene): str(num_gene) for num_gene in range(0,31,10)})
                        ])
            ], style={'margin-bottom': '30px'
                    #   'border': '2px grey solid'
                    }),
        
        html.Div([
            html.Label("Genes enriched in cell types: ", htmlFor="enriched_celltype", 
                    style={'font-weight': 'bold', 'margin-top': '30px'}),
            dcc.Checklist(id='enriched_celltype',
                            options=['GBC', 'INP', 'early iOSN', 'late iOSN', 'mOSN'],
                            value=[]),
            html.Label("Number of genes: ", htmlFor="enriched_celltype", 
                    style={'margin-top': '20px', 'margin-bottom': '20px'}),
            dcc.Slider(0,
                        30,
                        step=1,
                        id='n_enriched_genes',
                        value=10,
                        marks={str(num_gene): str(num_gene) for num_gene in range(0,31,10)})
        ], style={'margin-bottom': '30px'
                #   'border': '2px grey solid'
                  })
    ], style={'width': '25%', 
              'float': 'right', 
              'backgroundColor': '#FFFFFF'}),
    
    # Expression plot 
    html.Div([
        dcc.Graph(id='expression_plot', config={'displayModeBar': False})
    ],
            style={'width': '75%', 'height': '80%', 
                'float': 'left', 
                'display': 'block'}),
    
    html.Div([
        html.P("Discover gene expression patterns in Olfactory Sensory Neurons using our interactive tool. Choose a subject gene as your reference point, then explore gene trajectories in real-time. Highlight genes enriched in specific cell types and adjust the top and bottom associated genes to refine your analysis. Explore with uncover the dynamics of gene expression in OSNs developemnet.", 
                style={"margin": "50px 30px 30px 10px", 'font-size': '16px'}),
        html.H5("For source code visit ", style={'display': 'inline-block'}),
        html.A("geLine github", href="https://github.com/Justice-Lu/geLine_dash")
    ], style={"display": "inline-block", "float": 'left', 
              'background-color': 'white', 'font-size': '16px', 
              'width': '75%'})
], style={'backgroundColor': '#FFFFFF'} )

# Callback to update dropdown options and graph
@app.callback(
    Output('expression_plot', 'figure'),
    Input('sub_gene_dropdown', 'value'),
    Input('que_genes_dropdown', 'value'),
    Input('plot-type', 'value'),
    Input('top_associated_genes', 'value'),
    Input('bottom_associated_genes', 'value'), 
    Input('enriched_celltype', 'value'), 
    Input('n_enriched_genes', 'value')
)
def update_expression_plot(subject_gene, query_genes, plot_type, n_top_genes, n_bot_genes, enriched_celltype, n_enriched_genes):
    selected_dataset = 'ge_log1p' if 'logmaritize' in plot_type else 'ge_normalized'
    selected_dataset += '_normExp' if 'norm' in plot_type else '_'

    # print("DEBUG: ", enriched_celltype) 
    
    return expression_plot(subject_gene, query_genes, n_top_genes, n_bot_genes, enriched_celltype, n_enriched_genes, ge_data_list[selected_dataset])

if __name__ == '__main__':
    app.run_server(debug=True)
