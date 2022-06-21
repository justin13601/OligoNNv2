import numpy as np
from Bio.Seq import Seq
import json
from OligoNN.core import energy_matrix as energy, oligo_gen as oligo, oligo_database as db, \
    quant_energy as quant

from OligoNN.core import probe_seqs as inputs


def read_json(filename):
    """
    Reads a json file containing parameters for a specific neural network
    
    Parameters: filename - a json file formatted as a dictionary containing keys 'NN_shape', 'Gene_List', 'Weights', 'One-Hot_Encoded_Output', and 'Biases'
    Returns: variables for each value in the json dictionary
    """
    with open(filename) as json_file:
        file = json.load(json_file)
        shape = file["NN_Shape"]
        top_genes = file["Gene_List"]
        weights = file["Weights"]
        output_key = file["One-Hot_Encoded_Output"]
        biases = file["Biases"]
    return shape, top_genes, weights, output_key, biases


def gen_tf(anchor, tether, threshold):
    """
    Generates two transcription factors.
    
    Parameters: anchor - a dna sequence, tether - a dna sequence the same length as anchor,
    threshold - an integer that represents the energy difference between one tf and the cage bound to the tether, as well as the difference between both tfs and the cage bound to the tether
    Returns: two tf sequences, the static and the dynamic
    """
    disp_energy = energy.calc_energy(tether[10:])
    length = len(tether)
    min_energy = energy.calc_energy(tether)
    if disp_energy - threshold < min_energy:
        raise TypeError("The threshold energy is too large.")
    elif disp_energy + threshold > 0:
        raise TypeError("The threshold energy is too large.")
    for i in range(1, length):
        # energy of 1st alpha tf must be weaker than the cage energy by a certain threshold
        if energy.calc_energy(tether[i:]) > disp_energy + threshold:
            a_t1 = oligo.complement(tether[i:])
            b_t1 = oligo.complement(anchor[:len(a_t1)])
            t1 = a_t1 + b_t1
            for j in range(1, i):
                # energy of both alpha tfs must be stronger than the cage energy by a certain threshold
                if energy.calc_energy(tether[i - j:]) < disp_energy - threshold:
                    a_t2 = oligo.complement(tether[i - j:i])
                    b_t2 = oligo.complement(anchor[len(a_t1):len(a_t1) + len(a_t2)])
                    t2 = b_t2 + a_t2
                    return t1, t2


def shape_from_weights(top_genes, weights):
    """
    If the shape is not given in the json file, interprets it from the 3D weight matrix.
    
    Parameters: top_genes - a list of the probe names, weights - a 3D matrix or nested list of all the weights in a model
    Returns: a list with the number of nodes in each layer of the NN, including the input
    """
    shape = [len(top_genes)]
    for weight_matrix in weights:
        num_rows, num_cols = weight_matrix.shape
        shape.append(num_rows)
    return shape


def gen_nodes(modelfile, starting_genes):
    """
    Generates all of the dna sequences for a specific neural net model. All strands are BioPython Seq objects, 5'-3'.
    
    Parameters: modelfile - a json file with all the NN parameters, starting_genes - a dictionary of sequences to incorporate into the model
    Returns: a list where each entry is a dictionary for one layer of the neural net, in a dictionary the keys are the names of strands and the values are a list of dna sequences for each node in that layer,
    one node can be retrieved by indexing each dictioary value at the same position
    """
    # read json file with final model variables
    shape, top_genes, weights, output_key, biases = read_json(modelfile)

    # initialize database
    database = db.Database()

    # create list to store all layers
    NN = []

    # get input probe sequences
    input_seqs_df = inputs.probes_df(top_genes)
    # each layer is a dictionary with keys as names of strands and values as a list of seqs
    l_0 = {}
    probe_seqs = []
    for probe in input_seqs_df["Probe Sequences"]:
        index = 0
        size = database.size
        while database.size < size + 1:
            try:
                database.database_insert(Seq(probe[index]))
                index += 1
            # except block handles case that NONE of the probe sequences were accepted into the database
            # ***TEMPORARY FIX***
            except IndexError:
                index -= 1
                break
        probe_seqs.append(Seq(probe[index]))
    l_0["Probe Sequence"] = probe_seqs
    print("Layer 0: ", l_0)
    NN.append(l_0)

    # add the tether and promotor to the database
    database.database_insert(starting_genes["Tether"])
    database.database_insert(starting_genes["T7 Promoter"])

    # generate all the sequences for every node in each layer
    for layer in range(1, len(shape)):
        # add the cage and tether sequences to the layer dictionary
        l_i = {}
        l_i["Cage Sense"] = [starting_genes["Cage Sense"]] * shape[layer]
        l_i["Cage Antisense"] = [starting_genes["Cage Antisense"]] * shape[layer]
        l_i["Tether"] = [starting_genes["Tether"]] * shape[layer]

        print("getting anchor strands")
        tether_length = len(starting_genes["Tether"])
        size = database.size
        # generate anchor strands until all of them have been accepted into the database
        while database.size < size + shape[layer]:
            anchor = oligo.oligo(tether_length)
            database.database_insert(anchor)
        anchor_seqs = [Seq(x) for x in database.contents['Strand'][size:]]
        print("DONE")

        print("getting transcription factors")
        threshold_energy = 9  # variable that can be changed, pos integer, see gen_tf for description
        static_tf_seqs = []
        tf_seqs = []
        for anchor in anchor_seqs:
            static_tf, tf = gen_tf(anchor, starting_genes["Tether"], threshold_energy)
            static_tf_seqs.append(static_tf)
            tf_seqs.append(tf)
        print("DONE")

        print("getting outputs")
        output_length = 25  # length of dna transcript from one node
        size = database.size
        while database.size < size + shape[layer]:
            output = oligo.oligo(output_length).sequence
            database.database_insert(output)
        transcript_seqs = [Seq(x) for x in database.contents['Strand'][size:]]
        print("DONE")

        # assemble longer strands in the node
        l_i["Static TF + Transcript Sense"] = [static_tf_seqs[i] + starting_genes["T7 Promoter"] + transcript_seqs[i]
                                               for i in range(shape[layer])]
        l_i["Transcript Antisense + Anchor"] = [
            oligo.complement(transcript_seqs[i]) + oligo.complement(starting_genes["T7 Promoter"]) + anchor_seqs[i] for
            i in range(shape[layer])]

        # intermediates are the strands that determine weights in toehold-mediated displacement
        print("getting intermediate")
        toe_length = 20  # standard length for all toehold sequences
        # get the 2D matrix for this layer and round the values to one decimal place
        weight_matrix = np.array(weights[layer - 1])
        weight_matrix = np.round(weight_matrix, 1)
        intermediate_seqs = []
        tf_appendage_seqs = []
        for i in range(shape[layer - 1]):
            if layer == 1:
                output = NN[0]["Probe Sequence"][i]
            else:
                output = NN[layer - 1]["Static TF + Transcript Sense"][i][-output_length:]
            inters = []
            top_toe = output[:toe_length]
            b_dom = output[toe_length:]
            tf_appendage_seqs.append(b_dom)
            # get all the possible sequences for toehold weights between 0 and 1
            weight_dict = quant.find_quanta(top_toe)
            for j in range(shape[layer]):
                w = weight_matrix[j, i]
                tf = tf_seqs[j]
                a_star_tf = tf[:len(tf) // 2]
                if w < 0:
                    # negative weights
                    inters.append(a_star_tf + oligo.complement(b_dom) + weight_dict[w * -1])
                else:
                    # positive weights
                    inters.append(oligo.complement(a_star_tf) + oligo.complement(b_dom) + weight_dict[w])

            intermediate_seqs.append(inters)
        # each list in the nested list is for one node in the layer, get nodes row-wise
        l_i["Intermediate"] = np.array(intermediate_seqs).T.tolist()
        print("DONE")

        # TF and TF Inhibitor are products of toehold-mediated displacement for pos and neg weights, respectively
        full_tf_seqs_2D = []
        attack_seqs_2D = []
        for tf in tf_seqs:
            full_tf_seqs = []
            attack_seqs = []
            for appendage in tf_appendage_seqs:
                full_tf_seq = appendage + tf
                attack_seq = appendage + oligo.complement(tf[:len(tf) // 2])
                full_tf_seqs.append(full_tf_seq)
                attack_seqs.append(attack_seq)
            full_tf_seqs_2D.append(full_tf_seqs)
            attack_seqs_2D.append(attack_seqs)
        l_i["TF"] = full_tf_seqs_2D
        l_i["TF Inhibitor"] = attack_seqs_2D

        print("Layer {}: ".format(layer), l_i)
        # add the completed layer to the NN list
        NN.append(l_i)

    return NN
