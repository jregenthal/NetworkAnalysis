#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: johannaregenthal

Date: August 2020
"""

# Presettings ################################################################
# Import packages
import pandas as pd
import numpy as np
import random
import networkx as nx
from collections import Counter
import os
from spreading_CR import SpreadingProcess #https://github.com/gstonge/spreading_CR/blob/master/README.md

# User-defined functions
def create_random_networks(original_graph, number_of_random_graphs, networksize):
    """
    Function to create several random graphs based on one original graph.
    Returns list of random graphs
    """
    # Get degree sequence from original graph
    degree_sequence = [d for n, d in original_graph.degree()]

    # Create true configuration networks without self-loops and parallel edges
    N = 0
    graphs = []
    while N < number_of_random_graphs:
        #Sample from original degree_sequence
        deg_sequence_sample = random.sample(degree_sequence, networksize)
        #Check for valid degree sequence
        if nx.is_graphical(deg_sequence_sample):
            G_random = nx.configuration_model(deg_sequence_sample, create_using = nx.Graph())
            graphs.append(G_random)
            N += 1
    return graphs


def create_networkcopy_with_missing_nodes(graph, list_missing_percentages):
    """
    Function to return a list of graphs which have p% missing nodes based on a given list P.
        graph = true network
        list_missing_percentages = list that contains the percentages of which the observed networks should be created
    """
    lst = []
    for p in list_missing_percentages:
        G = graph.copy()
        G.remove_nodes_from(random.sample(list(G.nodes()), int(len(graph.nodes()) * p)))
        lst.append((G, p))
    return lst



def sim_infection(graph, immune_nodes, unimmune_nodes, numOfSimulations, beta, gamma, seed):
    """ 
    Function to simulate an several outbreaks based on the SIR model.
    Returns a tuple with (mean, std) of outbreak size
        graph = graph in which the outbreak is simulated
        immune_nodes = list of nodes that are immune
        unimmune_nodes = list of nodes that are unimmune
        numOfSimulations = number of simulations
        beta, gamma = parameters for the SIR model outbreak
        seed = similar seed for each simulation
    """
    final_size_list = []
    
    # patient0s should be nodes which have at least one neighbor, otherwise sp.initialize will throw an error, so we first check can be used
    nodes_with_neighbors = [n for n in unimmune_nodes if len(list(graph.neighbors(n))) > 0]
    array_patient0s = np.random.choice(nodes_with_neighbors, size=numOfSimulations)
    
    list_patient0s = [[patient0] for patient0 in array_patient0s] 
    del nodes_with_neighbors, array_patient0s
    
    sp = SpreadingProcess(list(graph.edges()), beta, gamma, 0)
    # Start simulations
    for i in range(numOfSimulations): 
        sp.initialize(list_patient0s[i], immune_nodes, seed)
        sp.evolve(np.inf)
        
        # Check how many nodes are infected and append to final_size_list
        finalSizeInclTargets = sp.get_Rnode_number_vector()[-1]
        adjustedFinalSize = finalSizeInclTargets -  len(immune_nodes)
        final_size_list.append(adjustedFinalSize) # final_size_list = Number of infected nodes per simulation
        
        sp.reset()
    
    # Outbreak size = % of nodes that are recovered/immune and were infected
    pctTotalPopInf_mean = np.mean(final_size_list) / graph.number_of_nodes() # outbreak size in % 
    pctTotalPopInf_std = np.std(final_size_list) / graph.number_of_nodes() # std of outbreak size
    
    return (pctTotalPopInf_mean, pctTotalPopInf_std)


# Define parameters ##########################################################

# Set working directory
os.chdir('/Users/johannaregenthal/Documents/GitHub/NetworkAnalysis/Rosenblatt2020')
MAX_CPLUSPLUS_INT = 4294967295

P = np.arange(0, 0.9, 0.1) #Part of missing nodes (error level 0 - 80%)
q = 0.1 #Fraction of immunized nodes = 10%
n = 1000 #Number of nodes = 1000
N = 50 #5000 # Number of random networks
num_SIR = 2000 #Number of simulations per network

# Outbreak parameters
beta = 0.95
gamma = 1

centralities = [('degree', nx.degree_centrality), 
                ('between', nx.betweenness_centrality), 
                ('eigen', nx.eigenvector_centrality_numpy), 
                ('page', nx.pagerank)]


# Creating random networks ###################################################

# Import network data
data = pd.read_csv('edgelist.truecolsprings.csv', usecols = ['V1', 'V2'])
nodes = list(set(data.V1).union(set(data.V2)))
edges = list(zip(data.V1, data.V2))

# Create Colorado Springs network
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)
#print(nx.info(G))

graphs = create_random_networks(G, N, n)

# delete unnecessary variables
del data, nodes, edges, G


# Simulations ################################################################
data_degree, data_between, data_eigen, data_page = [], [], [], []
# For each true network
for i, G_true in enumerate(graphs[:5]):
    # For each observed network with p% missing nodes
    for G_observed, p in create_networkcopy_with_missing_nodes(G_true, P):
        seed = np.random.randint(MAX_CPLUSPLUS_INT+1)
        num_immune_nodes = int(len(G_observed.nodes) * q)
        print('Graph: ', i, G_observed.number_of_nodes(), p)
        
        #Simulation with no immunization
        mean_noImm, std_noImm = sim_infection(G_observed, [], list(G_observed.nodes), num_SIR, beta, gamma, seed)
        
        #Random immunization
        immune_nodes = random.sample(list(G_observed.nodes), num_immune_nodes) # Select q% of nodes for immunization
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_ran, std_ran = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma, seed)
        
        # For each centrality
        for name, c in centralities:
            dic = c(G_observed)
            immune_nodes = sorted(dic, key=dic.get, reverse=True)[:num_immune_nodes]
            unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
            mean, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma, seed)
            exec(f'data_{name}.append([i, p, mean_ran, std_ran, mean, std])')
            

# Create dataframe for results of each centrality
for name, c in centralities:
    exec(f"data_{name} = pd.DataFrame(data_{name}, columns = ['graph_ID', 'p', ''mean_ran', 'std_ran', 'mean_strat', 'std_strat'])")
    exec(f"data_{name}.name = name")
    
del i, G_true, p, num_immune_nodes, mean_noImm, std_noImm, immune_nodes, unimmune_nodes, mean_ran, std_ran, dic, mean, std, name, c


# Analysis ###################################################################

# Pring mean for each error level for each centrality
for df in [data_degree, data_between, data_eigen, data_page]:
    df['robustness_Ros'] = df['mean_ran'] - df['mean_strat']
    print('\n', df.name)
    print(df.groupby('p')[['mean_ran', 'mean_strat', 'robustness_Ros']].mean())
    


#To Do for each centrality df: 
    #Calculate Rosenblatt-Robustness (difference to ran-immunization)
    #Include Niemeyer-Robustness





#Trying around

#mean_noImm, std_noImm = sim_infection(G_observed, [], list(G_observed.nodes), num_SIR, beta, gamma, seed)
graph = G_observed.copy()

immune_nodes = []
unimmune_nodes = list(G_observed.nodes)
numOfSimulations = num_SIR

