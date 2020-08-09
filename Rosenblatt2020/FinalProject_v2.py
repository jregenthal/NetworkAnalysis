#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: johannaregenthal

Date: August 2020

Goal: Experiment on random graphs
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

# Set working directory
os.chdir('/Users/johannaregenthal/Documents/GitHub/NetworkAnalysis/Rosenblatt2020')

# User-defined functions

def create_random_networks(original_graph, number):
    """
    Function to create several random graphs based on one original graph.
    Returns list of random graphs
    """
    # Get degree sequence from original graph
    degree_sequence = [d for n, d in original_graph.degree()]

    # Create true configuration networks without self-loops and parallel edges
    N = 0
    graphs = []
    while N < number:
        #Sample from original degree_sequence
        deg_sequence_sample = random.sample(degree_sequence, n)
        #Check for valid degree sequence
        if nx.is_graphical(deg_sequence_sample):
            G_random = nx.configuration_model(deg_sequence_sample, create_using = nx.Graph())
            graphs.append(G_random)
            N += 1
    return graphs


def create_networkcopy_with_missing_nodes(graph, list_missing_percentages):
    """
    Function to return a list of graphs which have missing nodes based on a given list.
        graph = true network
        list_missing_percentages = list that contains the percentages of which the observed networks should be created
    """
    lst = []
    for p in list_missing_percentages:
        G = graph.copy()
        G.remove_nodes_from(random.sample(list(G.nodes()), int(len(graph.nodes()) * p)))
        lst.append((G, p))
    return lst



def sim_infection(graph, immune_nodes, unimmune_nodes, numOfSimulations, beta, gamma):
    patient0s = [[patient0] for patient0 in np.random.choice(unimmune_nodes, size=numOfSimulations)] 
    final_size_list = []
    
    sp = SpreadingProcess(list(graph.edges()), beta, gamma, 0)
    
    # Start simulations
    for i in range(numOfSimulations): 
        #Inode_list = patient0s[i]
        sp.initialize(patient0s[i], immune_nodes, seed)
        sp.evolve(np.inf)
        
        finalSizeInclTargets = sp.get_Rnode_number_vector()[-1]
        adjustedFinalSize = finalSizeInclTargets -  len(immune_nodes)
        final_size_list.append(adjustedFinalSize) # final_size_list = Number of infected nodes per simulation
        sp.reset()
    
    #finalSizeDist = Counter(final_size_list) # Counting frequencies of number of infected nodes
    #finalSizeDist_noImm = dict(finalSizeDist) # Create dictionary of frequencies
    
    # Outbreak size = % of nodes that are recovered/immune and were infected
    pctTotalPopInf_mean = np.mean(final_size_list) / graph.number_of_nodes() # outbreak size in % 
    pctTotalPopInf_std = np.std(final_size_list) / graph.number_of_nodes() # std of outbreak size
    
    return [pctTotalPopInf_mean, pctTotalPopInf_std]
    #return [np.mean(final_size_list), np.std(final_size_list), graph.number_of_nodes()]


# Predefine parameters
seed = 42
P = np.arange(0, 0.5, 0.1) #Part of missing nodes (error level 0 - 80%)
q = 0.1 #Fraction of immunized nodes = 10%
n = 1000 #Number of nodes = 1000

# Outbreak parameters
beta = 0.95
gamma = 1

# centrality measure (functions)
# comparing to random 
#Degree
#Betweenness
#Eigenvector
#PageRank


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

graphs = create_random_networks(G, 1000)

# delete unnecessary variables
del data, nodes, edges, G


# Simulations ################################################################
listing = []
num_SIR = 2000
random.seed(42)
# For each true network
for G_true in graphs[:2]:
    # For each observed network with p% missing nodes
    for G_observed, p in create_networkcopy_with_missing_nodes(G_true, P):
        num_immune_nodes = int(len(G_observed.nodes) * q)
        print(G_observed.number_of_nodes(), p, num_immune_nodes)
        
        #Simulation with no immunization
        mean_noImm, std_noImm = sim_infection(G_observed, [], list(G_observed.nodes), num_SIR, beta, gamma)
        
        #Random immunization
        immune_nodes = random.sample(list(G_observed.nodes), num_immune_nodes) # Select q% of nodes for immunization
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_ran, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma)
        
        #Degree immunization
        dic_degree = nx.degree_centrality(G_observed) 
        immune_nodes = sorted(dic_degree, key=dic_degree.get, reverse=True)[:num_immune_nodes]
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_deg, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma)
        
        #Betweenness immunization
        dic_betweeness = nx.betweenness_centrality(G_observed) 
        immune_nodes = sorted(dic_betweeness, key=dic_degree.get, reverse=True)[:num_immune_nodes]
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_bet, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma)
        
        #Eigenvector immunization
        dic_eigen = nx.eigenvector_centrality(G_observed) 
        immune_nodes = sorted(dic_eigen, key=dic_degree.get, reverse=True)[:num_immune_nodes]
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_eig, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma)
        
        #Pagerank immunization
        dic_pagerank = nx.pagerank(G_observed) 
        immune_nodes = sorted(dic_pagerank, key=dic_degree.get, reverse=True)[:num_immune_nodes]
        unimmune_nodes = list(set(G_observed.nodes) - set(immune_nodes))
        mean_page, std = sim_infection(G_observed, immune_nodes, unimmune_nodes, num_SIR, beta, gamma)
     
        
        listing.append([G_observed.number_of_nodes(), p, mean_noImm, mean_ran, mean_deg, mean_bet, mean_eig, mean_page])

df = pd.DataFrame(listing, 
                  columns=['graphsize', 'p', 'mean_noImm', 'mean_ran', 'mean_deg', 'mean_bet', 'mean_eig', 'mean_page'])




#Next: robustness per 

# Calculate closeness centrality -> smallest average distance
dic_degree = nx.degree_centrality(G_true) 

# Get 10 nodes with highest 
for node in sorted(dic_degree, key=dic_degree.get, reverse=True)[:10]:
  print(node, dic_degree[node])