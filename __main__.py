from package.core import *
from Bio.Seq import Seq
import pandas as pd

# taken from 200128 RL001 Node Plan, these strands remain constant in all nodes
starting_genes = {"Tether": Seq("GCTACTCACTCAGATAGGTGACTGA"),
                  "Cage Sense": Seq("TCAGTCACCTATCTGTTTCAAATTAATACGACTCACTATA"),
                  "Cage Antisense": Seq("TATAGTGAGTCGTATTAATTTG"), "T7 Promoter": Seq("TAATACGACTCACTATAG")}

json_file = 'OligoNN/models/json/pytorch_model_12-24-12-3.json'
output_dir = 'OligoNN/'

# example run to store sequences in a csv
NN = gen_nodes(json_file, starting_genes)
for layer in range(len(NN)):
    pd.DataFrame.from_dict(NN[layer]).to_csv(
        output_dir + 'weight_to_seq_outputs/molecular_model_12-24-12-3_layer_{}.csv'.format(layer))
