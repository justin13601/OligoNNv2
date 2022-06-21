"""
By Ariel -- June 2020
Contains functions for computing "weights" based on binding energies.
"""

from OligoNN.core import oligo_gen as gen, nupack_wrap as nupack
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ast
import random


def calc_prob(seq1, seq2):
    """
    Determine the minimum free energy, secondary structure, and probability of a secondary complex.
    
    Parameters -- two dna sequences
    Returns -- a tuple with the mfe, secondary stracture, and equilibirium probability
    """
    mfe = nupack.mfe((seq1, seq2), material="dna")
    seq_struct, energy = mfe[0]
    prob = nupack.prob((seq1, seq2), seq_struct, material="dna")
    return (float(energy), seq_struct, prob)


def spectrum(anchor, compare):
    """
    Determine a weight (0-1) between two strands.
    Based on where binding energy of compare to anchor falls within range of energies from anchor+anchor and anchor+complement.
    
    Parameters -- two dna sequences
    Returns -- a weight between 0 and 1 representative of binding energy scale
    """
    anchor_comp = gen.complement(anchor)

    # value of 1 is assigned to energy btw the anchor seq and its complement
    # value of 0 is assigned to energy btw the anchor seq and a copy of itself
    energy_100 = calc_prob(anchor, anchor_comp)[0]
    energy_0 = calc_prob(anchor, anchor)[0]
    energy_user = calc_prob(anchor, compare)[0]

    energy_percent = (energy_user - energy_0) / (energy_100 - energy_0)
    return energy_percent


def rev_spectrum(anchor, weight):
    """
    Determine the value of the binding energy corresponding to a specific weight.
    Returns binding energy a strand would need with the anchor to map to that weight.
    
    Parameters -- a dna sequence and a number from 0 to 1
    Returns -- a binding energy 
    """
    anchor_comp = gen.complement(anchor)
    struct_100, energy_100, prob_100 = calc_prob(anchor, anchor_comp)
    struct_0, energy_0, prob_0 = calc_prob(anchor, anchor)
    bind = energy_0 + weight * (energy_100 - energy_0)
    # print("Desired binding energy is: ", bind)
    return bind


def find_quanta(strand):
    """
    Removes bases from 3' end of complement strand and records the change in weight.
    Attempts to find 'units' of bases that correspond to a weight change of 0.1.
    
    Parameters -- strand, a dna sequence 
    Returns -- a dictionary with weights from 1 to 0 as keys and sequences as values
    """
    # generate complement
    complement = gen.complement(strand)
    percent = round(spectrum(strand, complement), 1)
    sequences = {}
    i = 1
    new_percent = 1
    sequences[new_percent] = complement
    # remove bases
    while new_percent > 0:
        new_percent = round(spectrum(strand, complement[:-i]), 1)
        sequences[new_percent] = complement[:-i]
        i += 1

    # checks if all possible weights have been found and fills in the missing ones
    perfect = set([round(i / 10, 1) for i in range(11)])
    seq_dict = set(sequences.keys())
    if perfect.difference(seq_dict) != set():
        #         print(sequences)
        #         print(perfect.difference(seq_dict))
        for w in perfect.difference(seq_dict):
            k = random.randint(1, 4)
            toe = sequences[round(w - 0.1, 1)] + "".join(random.choices(["A", "T", "C", "G"], k=k))
            while round(spectrum(strand, toe), 1) != w:
                #                 print(round(spectrum(strand, toe), 1))
                k = random.randint(1, 4)
                toe = sequences[round(w - 0.1, 1)] + "".join(random.choices(["A", "T", "C", "G"], k=k))
            sequences[w] = toe

    return sequences


def weights_database(filename, length, num):
    """
    Creates a pandas dataframe of sequences and their truncated complements with associated weights.
    Exports dataframe to csv.
    
    Parameters -- filename to save csv as, length of seqs, and number of seqs to produce
    Returns -- a dataframe with 2 cols, one a seq, the other a dict with weights: truncated complements
    """
    frame = [["Anchor", "Weight: Tf Strand"]]
    for i in range(num):
        row = []
        tfs = {}
        # generate an anchor sequence
        anchor = str(gen.oligo(length, fenergy=False, fmp=False).sequence)
        row.append(anchor)
        # remove bases from 3' end of complement one at a time and add weight and seq to dict
        complement = gen.complement(anchor)
        for x in range(length):
            strand = complement[:len(complement) - x]
            # weight is rounded to 2 decimal places
            weight = round(spectrum(anchor, strand), 2)
            if weight > 0:
                tfs[weight] = strand
        row.append(tfs)
        frame.append(row)
    # create dataframe and export as csv
    data = np.array(frame)
    df = pd.DataFrame(data=data[1:], columns=data[0])
    df.to_csv(filename + ".csv")
    return df


def find_weights(filename, weight):
    """
    Searches through csv file created by weights_database to find pairs of seqs with a specific weight.
    
    Parameters -- filename (include .csv) to read and a weight (num from 0 to 1 w two decimals)
    Returns -- prints pairs of seqs (an anchor and a trucated complement) with the weight from the database
    """
    df = pd.read_csv(filename)
    for row in df.itertuples(name=None):
        tfs = ast.literal_eval(row[3])
        if weight in tfs:
            print("Anchor: ", row[2])
            print("Tf: ", tfs[weight])
        else:
            print("No weight match found.")


def get_pairs_matrix(seq1, seq2):
    """
    Parameters -- two dna seqs
    Returns -- matrix of base pair probabilities between the two seqs
    """
    pairs_matrix = np.zeros((len(seq1), len(seq2)))
    print('getting pairs matrix...')
    pairs = nupack.pairs((seq1, seq2), material="dna")
    for i in pairs:
        if i[1] != len(seq1) + len(seq2) + 1 and i[1] > len(seq1) and i[0] <= len(seq1):
            pairs_matrix[(i[0]) - 1, (i[1]) - len(seq1) - 1] = i[2]
    print("done.")
    return pairs_matrix


def plot_pairs(pairs_matrix):
    """
    Plots the get_pairs_matrixs
    
    Parameters -- a matrix of bp probs
    Returns -- plots the matrix as a "heat map" of probs at positions
    """
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=4)

    fig, ax = plt.subplots()

    v = np.linspace(pairs_matrix.min(), pairs_matrix.max(), 5, endpoint=True)

    cax = ax.imshow(pairs_matrix, cmap='hot_r', interpolation='None')
    cbar = fig.colorbar(cax, ticks=v)
    cbar.locator = tick_locator
    ax.set_xlabel('Sequence 2 Position')
    ax.set_ylabel('Sequence 1 Position')
    plt.show()
