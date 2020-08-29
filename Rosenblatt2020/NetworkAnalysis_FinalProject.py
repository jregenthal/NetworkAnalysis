#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: johannaregenthal

Date: August 2020
"""

# 1. Presettings ################################################################

import os
# Set working directory to path where the script is
os.chdir('/Users/johannaregenthal/Documents/GitHub/NetworkAnalysis/Rosenblatt2020')

# Import required packages 
import pandas as pd
import numpy as np
import random
import networkx as nx
from spreading_CR import SpreadingProcess #https://github.com/gstonge/spreading_CR/blob/master/README.md
from EstimatingCentralityRobustness import robustness_calculator_builder #https://github.com/crsqq/EstimatingCentralityRobustness


# Define user-specific functions
def create_random_networks(original_graph, number_of_random_graphs, networksize):
    """
    Function to generate a defined number of random graphs with a specified networksize n based on one original graph. 
    More precisely, the function performs a configuration model which creates random graphs from a given degree sequence.
    The function returns a list of random graphs.
    """
    # Get degree sequence from original graph
    degree_sequence = [d for n, d in original_graph.degree()]

    # Create true configuration networks with the same degree sequence without self-loops and parallel edges
    N = 0
    graphs = []
    while N < number_of_random_graphs:
        # Sample n nodes from original degree sequence at random
        deg_sequence_sample = random.sample(degree_sequence, networksize)
        # Check for valid degree sequence 
        if nx.is_graphical(deg_sequence_sample):
            # Generate random graph based on the sampled nodes and append to list
            G_random = nx.configuration_model(deg_sequence_sample, create_using = nx.Graph())
            graphs.append(G_random)
            N += 1
    return graphs


def create_networkcopy_with_missing_nodes(graph, list_missing_percentages):
    """
    Function to generate copies of a graph with varying proportion p of missing nodes given by a list P.
    The function returns a list of graphs, where each graph is a network copy with p% of the nodes missing.
        graph = true network
        list_missing_percentages = list that contains different values for p
    """
    lst = []
    for p in list_missing_percentages:
        # Generate copy of the graph, randomly remove p% of the nodes from the set and append to list
        G = graph.copy()
        G.remove_nodes_from(random.sample(list(G.nodes()), int(len(graph.nodes()) * p)))
        lst.append((G, p))
    return lst


def sim_infection(graph, immune_nodes, numOfSimulations, beta, gamma):
    """ 
    Function to perform a defined number of outbreak simulations based on the SIR model in a given graph under 
    consideration of a defined set of immunized nodes.
    The function returns a tuple with mean and standard deviation (std) of the outbreak size over all simulations, 
    where the outbreak size is given by the proportion of infected nodes in the graph. 
        graph = graph in which the outbreak is simulated
        immune_nodes = list of nodes that are immune
        numOfSimulations = number of simulations
        beta, gamma = outbreak parameters for the SIR model
    """
    final_size_list = [] # List containing final size of infected nodes per simulation
    
    # Create list of unimmune nodes
    unimmune_nodes = list(set(graph.nodes) - set(immune_nodes))
    
    # Identify potential starting nodes of the outbreak simulation (= patient0s) by checking for nodes with at least
    # one neighbor. Otherwise, the function sp-initialize throws an error
    nodes_with_neighbors = [n for n in unimmune_nodes if len(list(graph.neighbors(n))) > 0]
    
    # Draw patient0s from the nodes_with_neighbors for each simulation at random
    array_patient0s = np.random.choice(nodes_with_neighbors, size=numOfSimulations)
    
    # Convert each element to an own list due to the form of the sp.initialize function
    list_patient0s = [[patient0] for patient0 in array_patient0s] 
    
    # Create simulation object
    sp = SpreadingProcess(list(graph.edges()), beta, gamma, 0)
    
    # Start simulations for numOfSimulations runs
    for i in range(numOfSimulations): 
        seed = np.random.randint(MAX_CPLUSPLUS_INT+1)
        
        # Initialize simulation with infected patient0 and immune_nodes
        sp.initialize(list_patient0s[i], immune_nodes, seed)
        
        # Run outbreak simulation infinitly (without time constraint)
        sp.evolve(np.inf)
        
        # Check the number of recovered nodes (including recovered and immune nodes in terms of the report) and 
        # store in finalSizeInclTargets
        finalSizeInclTargets = sp.get_Rnode_number_vector()[-1]
        
        # Subtract the number of immune nodes from finalSizeInclTargets to obtain the number of infected nodes
        adjustedFinalSize = finalSizeInclTargets -  len(immune_nodes)
        
        # Append everything to final_size_list containing the number of infected nodes per simulation
        final_size_list.append(adjustedFinalSize) 
        
        # Reset spreading process
        sp.reset()
    
    # Calculate mean and standard deviation of outbreak size over all simulations by using final_size_list
    pctTotalPopInf_mean = np.mean(final_size_list) / graph.number_of_nodes() # mean of outbreak size in % 
    pctTotalPopInf_std = np.std(final_size_list) / graph.number_of_nodes()   # std of outbreak size
    
    return (pctTotalPopInf_mean, pctTotalPopInf_std)


# 2. Definition of parameters ###################################################

MAX_CPLUSPLUS_INT = 4294967295 # Seed similar to original 

P = np.arange(0, 0.9, 0.1) # Part of missing nodes (nine different error levels between 0 and 80%)
q = 0.1                    # Fraction of immunized nodes = 10%
n = 1000                   # Number of nodes = 1000
N = 5000                   # Number of random networks = 5000
numOfSimulations = 2000    # Number of simulations per random network = 2000
beta = 0.95; gamma = 1     # Outbreak parameters of the SIR model

# List of examined centralities to loop over
centralities = [('degree', nx.degree_centrality), 
                ('between', nx.betweenness_centrality), 
                ('eigen', nx.eigenvector_centrality_numpy), 
                ('page', nx.pagerank)]


# 3. Creation of random networks ################################################

# Import network data
data = pd.read_csv('edgelist.truecolsprings.csv', usecols = ['V1', 'V2'])
nodes = list(set(data.V1).union(set(data.V2)))
edges = list(zip(data.V1, data.V2))

# Create Colorado Springs network (= Original network) 
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)

# Create N random graphs of size n from the Colorado Springs network G (= True networks) 
graphs = create_random_networks(G, N, n)

# Delete unnecessary variables from network data
del data, nodes, edges, G


# 4. Simulations ################################################################

# Initialize list to store data for each centrality measure
data_degree, data_between, data_eigen, data_page = [], [], [], []

# Define number of nodes to be immunized (q% of the n nodes in the true network)
num_immune_nodes = int(n * q)

# For each true network G_true 
for i, G_true in enumerate(graphs[:5]):
    
    # For each observed network G_observed with error level p (p% of nodes missing)
    for G_observed, p in create_networkcopy_with_missing_nodes(G_true, P):

        # Print network ID and error level p
        print('G_true: ', i, 'p: ', p)
        
        # Random immunization: Perform outbreak simulation in G_true with randomly immunized nodes sampled from G_observed
        # and obtain mean and std of the outbreak size
        immune_nodes = random.sample(list(G_observed.nodes), num_immune_nodes)
        mean_ran, std_ran = sim_infection(G_true, immune_nodes, numOfSimulations, beta, gamma)
        
        # For each centrality measure 
        for name, centrality in centralities:
            
            # Strategic immunization: Perform outbreak simulation in G_true with strategically immunized nodes from 
            # G_observed determined by the respective centrality measure and obtain mean and std of the outbreak size
            dic = centrality(G_observed)
            immune_nodes = sorted(dic, key=dic.get, reverse=True)[:num_immune_nodes]
            mean_strat, std_strat = sim_infection(G_true, immune_nodes, numOfSimulations, beta, gamma)
            
            # Measure by Rosenblatt et al. (2020): Use the concept proposed by Rosenblatt et al. (2020) to compute 
            # robustness of the centrality measure with respect to the corresponding immunization strategy 
            robustness_Ros = mean_ran - mean_strat
            
            # Measure by Martin and Niemeyer (2019): Use the function proposed by Martin and Niemeyer (2019) to compute 
            # robustness of the centrality measure
            robustness_calculator = robustness_calculator_builder(centrality)
            robustness_Nie = robustness_calculator(G_true, G_observed)            
            
            # Save relevant values in list 
            exec(f'data_{name}.append([i, p, mean_ran, std_ran, mean_strat, std_strat, robustness_Nie, robustness_Ros])')
            

# Create dataframe for the simulation results of each centrality measure
for name, c in centralities:
    exec(f"data_{name} = pd.DataFrame(data_{name}, columns = ['graph_ID', 'p', 'mean_ran', 'std_ran', 'mean_strat', 'std_strat', 'robustness_Nie', 'robustness_Ros'])")
    exec(f"data_{name}.name = name")

# Save results as csv
for df in [data_degree, data_between, data_eigen, data_page]:
    df.to_csv('data_' + df.name + '.csv', sep = ';')


# 5. Analysis ###################################################################

# see R skript "NetworkAnalysis_FinalProject_Plotting.R"
