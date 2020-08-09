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
def create_networks_with_missing_nodes(graph, list_missing_percentages):
    """
    Function to return a list of graphs which have missing nodes based on a given list.
        graph = true network
        list_missing_percentages = list that contains the percentages of which the observed networks should be created
    """
    lst = []
    for p in list_missing_percentages:
        G = graph.copy()
        G.remove_nodes_from(random.sample(list(G.nodes()), int(len(graph.nodes()) * p)))
        lst.append(G)
    return lst

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

    

# Predefine parameters
P = np.arange(0, 0.9, 0.1) #Part of missing nodes (error level 0 - 80%)
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

graphs = create_random_networks(G, 500)

# Get degree sequence from Colorado Springs network
degree_sequence = [d for n, d in G.degree()]

# Create true configuration networks without self-loops and parallel edges
N = 0
graphs = []
while N < TrueNetworks:
    #Sample from original degree_sequence
    deg_sequence_sample = random.sample(degree_sequence, n)
    
    #Check for valid degree sequence
    if nx.is_graphical(deg_sequence_sample):
        G_random = nx.configuration_model(deg_sequence_sample, create_using = nx.Graph())
        graphs.append(G_random)
        N += 1

# delete unnecessary variables
del graphs


# Simulations ################################################################
# For each true network
for G in graphs[:2]:
    # For each observed network with p% missing nodes
    for G_observed in create_networks_with_missing_nodes(G, P):
        # Select q% of nodes for immunization
        imNodes = random.sample(list(G_observed.nodes), int(len(G_observed.nodes) * q))
        #print(imNodes)
        





# Creation of graph
G = nx.configuration_model(deg_sequence_sample, create_using = nx.Graph())
listEdges = list(G_random.edges())
listNodes = list(G_random.nodes())
unimmunizedNodes = listNodes[:]

# Initialize number of outbreak simulations
seed = 42
numOfSimulations = 2000
patient0s = [[patient0] for patient0 in np.random.choice(unimmunizedNodes, size=numOfSimulations)] 
Rnode_list = [] #no previous infections
final_size_list = []

sp = SpreadingProcess(listEdges, beta, gamma, 0)

# Start simulations
for i in range(numOfSimulations): 
    Inode_list = patient0s[i]
    sp.initialize(Inode_list, Rnode_list, seed)
    sp.evolve(np.inf)

    finalSizeInclTargets =  sp.get_Rnode_number_vector()[-1] 
    final_size_list.append(finalSizeInclTargets) # In this case there are no targets so all recovered nodes were infected, none were immunized
    sp.reset()

finalSizeDist = Counter(final_size_list)
finalSizeDist_noImm = dict(finalSizeDist)
pctTotalPopInf_noImm = np.mean(final_size_list) / G.number_of_nodes()









degDictG = {n: d for n, d in G.degree()}     

# Estimate R0: analytically

degs = list(degDictG.values())
deg_sqrds = [deg**2 for deg in degs]

deg_avg = np.mean(degs)# frequeently written as <k>
deg_sqr_avg = np.mean(deg_sqrds) # frequently written as <k^2>

R0_formula = (beta/(beta+gamma)) * (deg_sqr_avg/deg_avg -1) #the second term is the "excess degree"


# Estimate R0: empirically
seed = 42

R0_empirical_mean, R0_empirical_std = sp.estimate_R0(numOfSimulations, seed) 
sp.reset()

print("R0 formula: ", R0_formula)
print("R0 emprical: ", R0_empirical_mean)
print("R0 empirical std: ",R0_empirical_std)


