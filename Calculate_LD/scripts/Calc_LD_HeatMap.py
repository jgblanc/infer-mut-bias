import msprime
import sys 
import numpy as np
import matplotlib.pyplot as pl
import os
from os import path
import csv
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

## Check Parameter inputs and assign name 
sample_size = int(parameters[0])
Ne = float(parameters[1])
length = float(parameters[2])
recombination_rate = float(parameters[3])
mutation_rate = float(parameters[4])
random_seed = float(parameters[5])
if parameters[6] != "None":
    start_time = int(parameters[6])
else:
    start_time = None
if parameters[7] != "None":
    end_time = int(parameters[7])
else:
    end_time = None
bin_size = int(parameters[8])
max_variants = int(parameters[9])
max_DAF = float(parameters[10])


## Functions 

# Simulate directly using MS Prime 
def msprime_simulate(sample_size, Ne, length, recombination_rate, mutation_rate, random_seed, start_time, end_time): 
    tree_sequence = msprime.simulate(
     sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate)
    ts_mutation = msprime.mutate(tree_sequence, rate=mutation_rate, random_seed=random_seed, start_time=start_time, end_time=end_time)
    return ts_mutation

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
def bin_frequencies(all_reps, num_bins):
    bins = np.linspace(0,1, num_bins+1)
    out = []
    for j in range(len(bins)-1):
        r = np.array(1)
        for k in range(len(all_reps)):
            dat = all_reps[k]
            for i in range(np.shape(dat)[0]):
                f = dat[i][0]
                if (f >= bins[j]) & (f < bins[j+1]):
                    r = np.append(r,dat[i][1:])
        if r.size > 1:
            r = np.append(bins[j],r[1::])
            out.append(r)
    return out 


# Bin r values with allele frequency bin 
def bin_r(out, num_bins):
    bins = np.linspace(-1,1, num_bins+1)
    freq_r_count = []
    for dat in out:
        freq = dat[0]
        for j in range(len(bins)-1):
            counter = 0
            f = np.array(1)
            for i in range(len(dat)):
                r = dat[i]
                if (r >= bins[j]) & (r < bins[j+1]):
                    counter = counter + 1
            freq_r_count.append([round(freq,2),round(bins[j],2),counter])  
    return freq_r_count



## Main program 
all_reps = []
num_hf_variants = 0
tot_variants = 0 
i = 0 
while num_hf_variants < max_variants:
    ts = msprime_simulate(sample_size, Ne, length, recombination_rate, mutation_rate, random_seed+i, start_time, end_time)
    G = ts.genotype_matrix()
    Genotype = get_diploids(G)
    print(np.shape(G)[0])
    tot_variants = np.shape(G)[0] + tot_variants
    cc = np.corrcoef(np.transpose(Genotype))
    freq = calc_freq(G)
    num_hf_variants = (freq > max_DAF).sum() + num_hf_variants 
    dat = shape_data(freq, cc)
    all_reps.append(dat)
    i = i+1
out = bin_frequencies(all_reps, bin_size)
freq_r_count = bin_r(out, bin_size)
a = np.array(freq_r_count)
print(freq_r_count)
print(np.sum(a, axis=0))
print(tot_variants)

## Generate and save plots 
if not os.path.isdir(outpath):
    os.makedirs(outpath)
cp = os.path.join(outpath, "HeatMapData.csv") 
with open(cp,"w") as f:
    wr = csv.writer(f)
    wr.writerows(freq_r_count)

    

