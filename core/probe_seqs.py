import pandas as pd
from Bio import SeqIO
import csv


# top_genes = ['204439_at', '203153_at', '213797_at', '219519_s_at', '207329_at', '212768_s_at',
#              '205569_at', '214059_at', '202411_at', '200986_at', '209795_at', '204747_at', 'AFFX-HUMRGE/M10098_5_at']

def probes_df(top_genes, FASTA=False):
    """
    returns a dataframe with the probe names, list of all probe sequences from AFFY documentation, and optional FASTA gene sequences
    if FASTA = True, must have FASTA files downloaded with filenames {probe_name}.fasta
    """
    if FASTA:
        # put the FASTA gene sequences into a list
        seqs = []
        for gene in top_genes:
            record = SeqIO.read("{}.fasta".format(gene), "fasta")
            seqs.append(record.seq)

    # put the probe sequences into a list of lists
    AFFY_probes = []
    with open("OligoNN/data/HG-U133A_2.probe_tab") as probes:
        probe_reader = csv.reader(probes, delimiter="\t")
        # clean the data to only the probes in the top genes
        edited_probes = []
        for probe in probe_reader:
            if probe[0] in top_genes:
                edited_probes.append(probe)
        # make a list entry with a list of probe seqs for each gene
        for gene in top_genes:
            gene_probes = []
            for entry in edited_probes:
                if entry[0] == gene:
                    gene_probes.append(entry[4])
            AFFY_probes.append(gene_probes)

    data = {}
    data["Probes"] = top_genes
    if FASTA:
        data["Gene Sequences"] = seqs
    data["Probe Sequences"] = AFFY_probes

    df = pd.DataFrame(data)

    return df
