import msprime
import sys 
import numpy as np
import matplotlib.pyplot as pl
import os
from os import path
np.random.seed(12)

## Parse Inputs 

parameter_inputs = sys.argv[1]
folder_name = parameter_inputs.split("/")[2]
folder_name = folder_name.split(".")[0]
parameter_inputs = sys.argv[1]
parameters = []
for line in open(parameter_inputs, "r"):
    line = line.strip()
    parameters.append(line)

outpath = "../output/" + folder_name
    
## Functions 

# Simulate directly using MS Prime 
def msprime_simulate(sample_size, Ne, length, recombination_rate, mutation_rate, random_seed): 
    tree_sequence = msprime.simulate(
     sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate,
     mutation_rate=mutation_rate, random_seed=random_seed)
    return tree_sequence

# Sample 2 random haploids without replacement to get diploids
def get_diploids(G):
    Geno = np.ones(np.shape(G)[0])
    ids = np.arange(np.shape(G)[1])
    while len(ids) >= 2:
        random_index = np.random.choice(range(len(ids)), 2, replace=False)
        col1 = G[:,ids[random_index[0]]]
        col2 = G[:,ids[random_index[1]]]
        ids = np.delete(ids, random_index)
        geno = col1 + col2
        Geno = np.vstack((Geno, geno))
    return(Geno[1:,:])

# Calculate allele frequency in sample 
def calc_freq(G): 
    return np.sum(G, axis=1)/np.shape(G)[1]

# Shape data - first column is allele frequency followed by all the r^2 values for all that variants 
def shape_data(freq, cc):
    dat = np.ones((len(freq), np.shape(G)[0]+1))
    for i in range(len(freq)):
        dat[i,:] = np.append(np.array(freq[i]), cc[i])
    return dat

# Bin frequencies 
def bin_frequencies(dat, num_bins):
    bins = np.linspace(0,1, num_bins+1)
    out = []
    for j in range(len(bins)-1):
        r = np.array(1)
        for i in range(np.shape(dat)[0]):
            f = dat[i][0]
            if (f >= bins[j]) & (f < bins[j+1]):
                r = np.append(r,dat[i][1:])
        if r.size > 1:
            r = np.append(bins[j],r[1::])
            out.append(r)
    return out 


## Main program 

ts = msprime_simulate(int(parameters[0]), float(parameters[1]), float(parameters[2]), float(parameters[3]), float(parameters[4]), int(parameters[5]))
G = ts.genotype_matrix()
Genotype = get_diploids(G)
cc = np.corrcoef(np.transpose(Genotype))
freq = calc_freq(G)
dat = shape_data(freq, cc)
out = bin_frequencies(dat, int(parameters[6]))
print(len(out))

## Generate and save plots 
if not os.path.isdir(outpath):
    os.makedirs(outpath)

# SFS
num_variants = np.shape(G)[0]
fig = pl.figure()         
ax = fig.add_subplot()   
ax.hist(freq, bins=20)
pl.xlabel(num_variants)
fig.savefig(path.join(outpath,"SFS.png"))
    
# LD Plots    
for i in range(len(out)):
    print(i)
    x = out[i]
    fig = pl.figure()         
    ax = fig.add_subplot()   
    ax.hist(x[1::], bins=100)
    pl.xlabel('r')
    fig.savefig(path.join(outpath,"freq_{0}.png".format(x[0], 1)))


        
