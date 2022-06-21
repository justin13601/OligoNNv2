#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 21:55:17 2017

@author: leo
"""

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 09:52:35 2017

@author: leo
"""
"""
Functions to produce and plot energy matrices, as well as calculate binding energy

Edits by Ariel
line 26-28: import my version of nupack, oligo_gen, and BioPython
calc_energy: fix vars to link to oligo_gen, allowed for seq and str objs as input
get_barcodes: updated how sequences are generated with oligo_gen
get_energy_matrix: allowed for a list of seq or str objs as input
"""

import numpy as np
import matplotlib.pyplot as plt
from OligoNN.core import oligo_gen as gen, nupack_wrap as nupack
from Bio.Seq import Seq, reverse_complement
import operator


def calc_energy(seq):
    """
    Calculates the binding energy of a sequence to its complement.
    
    Parameters: seq - a dna sequence
    Returns: the binding energy 
    """
    # takes in a Seq or str obj and returns its energy
    seq_struct = len(seq) * '(' + '+' + len(seq) * ')'
    seq_energy = nupack.energy((seq, gen.complement(seq)), seq_struct, material='dna')
    return seq_energy


def get_barcodes(len_dna, num_dna):
    """
    generate sequence using get_seq, calcualte Tms of candidates in pool, remove underperformers until criteria matched.
    I want a bigger initial pool that is independent of num_dna to get decent sampling
        
    returns: sorted_barcode_diffs = a list of dna sequences
    """

    num_barcodes, barcodes, barcode_diffs = 0, {}, {}

    # generate candidate sequences and binding energies
    print('generating initial sequence library...'),
    while num_barcodes <= num_dna:
        cand_seq = gen.oligo(len_dna).sequence
        cand_energy = calc_energy(cand_seq)
        # note - barcodes keys are seq objects
        barcodes[cand_seq] = cand_energy
        num_barcodes = len(barcodes)
    print('done.')

    print('optimizing sequences...'),

    # optmize sequence 100 times
    for _ in range(10):
        median = np.median(list(barcodes.values()))
        barcode_diffs = {seq: abs(median - energy) for (seq, energy) in barcodes.items()}
        del barcodes[max(barcode_diffs, key=barcode_diffs.get)]

        # add new sequence to dictionary
        cand_seq = gen.oligo(len_dna).sequence
        cand_energy = calc_energy(cand_seq)
        barcodes[cand_seq] = cand_energy

    # return the number of barcodes asked for
    sorted_barcode_diffs = sorted(barcode_diffs, key=operator.itemgetter(1))[
                           :num_dna]  # returns the top candidate sequences
    # barcodes_to_return = {seq : energy for (seq, energy) in barcodes.items() if seq in sorted_barcode_diffs} # fetch energy of top candidates

    print('done.')
    return sorted_barcode_diffs


def get_energy_matrix(seqs):
    """
    evaluate binding energy for each seq in seqs, generate a 2D energy matrix
    as proxy for sequence orthogonality
    
    seqs = list of DNA sequences
    returns:
        energy_matrix: a 2D np.array of binding energies (diagnoal = perfect complement), Units = kcal/mol
    """

    # seqs should be a list of Biopython seqs or strings
    # note - get_barcodes returns a list of BioPython seqs
    if type(seqs[0]) == Seq:
        seqs_rc = [seq.reverse_complement() for seq in seqs]
    elif type(seqs[0]) == str:
        seqs_rc = [reverse_complement(seq) for seq in seqs]
    else:
        raise TypeError("values in seqs must be of type seq or str.")

    energy_matrix = np.zeros((len(seqs), len(seqs_rc)))

    print('getting energy matrix...'),
    for idx, seq in enumerate(seqs):
        # use NUPACK's mfe to calculate all pairwise binding energies
        # returns a list of tuples ('structure', 'energy)
        dna_mfe = [float(nupack.mfe((seq, rc), material='dna')[0][1]) for rc in seqs_rc]
        energy_matrix[idx, :] = dna_mfe  # store results in 2D array for plotting
        # seq_energys[idx] = (seq, dna_mfe) # store results into another dictionary for sorting; idx s

    print('done.')
    # return seq_energys, energy_matrix
    return energy_matrix


def get_energy_matrix_WW(seqs, material):
    """
    evaluate binding energy for each seq in seqs, generate a 2D energy matrix
    as proxy for sequence orthogonality. This version only calculates Watson-
    Watson interactions.
    
    seqs = list of DNA sequences
    returns: energy_matrix: a 2D np.array of binding energies (diagnoal = perfect complement), Units = kcal/mol
    """
    # seqs_rc = [utils.revcom(seq) for seq in seqs]
    dim = len(seqs)
    energy_matrix = np.zeros((dim, dim))

    print('getting energy matrix...'),
    for idx, seq in enumerate(seqs):
        # use NUPACK's mfe to calculate all pairwise binding energies
        # returns a list of tuples ('structure', 'energy)
        dna_mfe = [float(nupack.mfe((seq, rc), material=material)[0][1]) for rc in seqs]
        energy_matrix[idx, :] = dna_mfe  # store results in 2D array for plotting
        # seq_energys[idx] = (seq, dna_mfe) # store results into another dictionary for sorting; idx s

    print('done.')
    # return seq_energys, energy_matrix
    return energy_matrix


def get_energy_matrix_material(seqs, material):
    """
    evaluate binding energy for each seq in seqs, generate a 2D energy matrix
    as proxy for sequence orthogonality
    
    seqs = list of DNA sequences
    returns:
        energy_matrix: a 2D np.array of binding energies (diagnoal = perfect complement), Units = kcal/mol
    """
    if type(seqs[0]) == Seq:
        seqs_rc = [seq.reverse_complement() for seq in seqs]
    elif type(seqs[0]) == str:
        seqs_rc = [reverse_complement(seq) for seq in seqs]
    else:
        raise TypeError("values in seqs must be of type seq or str.")

    energy_matrix = np.zeros((len(seqs), len(seqs_rc)))

    print('getting energy matrix...'),
    for idx, seq in enumerate(seqs):
        # use NUPACK's mfe to calculate all pairwise binding energies
        # returns a list of tuples ('structure', 'energy)
        dna_mfe = [float(nupack.mfe((seq, rc), material=material)[0][1]) for rc in seqs_rc]
        energy_matrix[idx, :] = dna_mfe  # store results in 2D array for plotting
        # seq_energys[idx] = (seq, dna_mfe) # store results into another dictionary for sorting; idx s

    print('done.')
    # return seq_energys, energy_matrix
    return energy_matrix


"""plotting and visualization tools
    1. plot_energy
"""


def plot_energy(energy_matrix, save_fig=False, fname=None):
    """
    function to visualize energy matrix using imshow
    """
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=4)

    fig, ax = plt.subplots()

    v = np.linspace(energy_matrix.min(), energy_matrix.max(), 5, endpoint=True)

    cax = ax.imshow(energy_matrix, cmap='hot', interpolation='None')
    cbar = fig.colorbar(cax, ticks=v)
    cbar.locator = tick_locator
    ax.set_xlabel('Sequence Complement')
    ax.set_ylabel('Sequence No.')
    plt.show()

    if save_fig and not fname == None:
        fig.savefig(fname, bbox_inches='tight')
    return
