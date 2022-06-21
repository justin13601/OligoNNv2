"""
Generates a database of oligos, sorted by length

Edits by Ariel:
- allowed Biopython Seq objects or strings to be inserted into the database
"""

from OligoNN.core import *
import pandas as pd
from Bio.Seq import Seq
from .oligo_gen_string import Oligo

# define
complement_filter_size = 5  # 5 or longer complement nt sequences not allowed (ie. strands can't stick together)
sort_condition = 1  # 1 is decreasing length order, 0 is increasing length order


class Database:
    def __init__(self):
        self.size = 0
        self.contents = {'Strand': [], 'Numbers': [], 'Length': [], 'Strand Energy': [], 'Binding Energy': []}
        self.database = pd.DataFrame(self.contents)

    def __str__(self):
        return str(self.database)

    def update_strand_energy(self):
        self.contents['Strand Energy'] = []
        for i in range(0, len(self.contents['Strand'])):
            self.contents['Strand Energy'].append(calc_energy(self.contents['Strand'][i]))
        return True

    def update_energy_matrix(self):
        self.contents['Binding Energy'] = []
        matrix = get_energy_matrix(self.contents['Strand'])
        for i in range(0, len(matrix)):
            result = []
            for j in range(0, len(matrix[0])):
                result.append((self.contents['Strand'][i], self.contents['Strand'][j], matrix[i][j]))
            self.contents['Binding Energy'].append(result)

    def database_insert(self, strand_input):
        """
        Calls insert_helper to insert oligonucleotides into database
        :param strand_input: (list of Oligo/oligo objects OR singular Oligo/oligo object)
        :return: None
        """
        if isinstance(strand_input, list):  # Input Type: List of objects
            for i in range(len(strand_input)):
                if not self.insert_helper(strand_input[i]):
                    print("Strand not inserted.")
                    return
                self.size += 1
        else:  # Input Type: One object
            if not self.insert_helper(strand_input):
                print("Strand not inserted.")
                return
            self.size += 1
        self.database = pd.DataFrame(self.contents)
        self.reorder_by_length()

    def insert_helper(self, oligo_strand):
        """
        Inserts oligonucleotides by adding each strand into dictionary that is converted into dataframe using pandas
        :param oligo_strand: (Oligo/oligo object)
        :return: bool - True is inserted, False is not inserted (ie. same strands or sticking strands present)
        """
        if isinstance(oligo_strand, oligo):  # Ariel's oligo objects (biopython)
            if self.check_for_same_strand(oligo_strand) or self.check_for_sticking(oligo_strand) or check_hamming(
                    self.contents["Strand"], oligo_strand.sequence):
                return False
            self.contents['Strand'].append(str(oligo_strand.sequence))
            self.contents['Numbers'].append(oligo_strand.numbers)
            self.contents['Length'].append(oligo_strand.length)
            self.update_energy_matrix()
            self.update_strand_energy()
            return True
        elif isinstance(oligo_strand, Oligo):  # Justin's Oligo objects (numbers)
            if self.check_for_same_strand(oligo_strand) or self.check_for_sticking(oligo_strand) or check_hamming(
                    self.contents["Strand"], oligo_strand.sequence):
                return False
            self.contents['Strand'].append(oligo_strand.sequence)
            self.contents['Numbers'].append(oligo_strand.numbers)
            self.contents['Length'].append(oligo_strand.length)
            self.update_energy_matrix()
            self.update_strand_energy()
            return True
        elif isinstance(oligo_strand, Seq):
            if self.check_for_same_strand(oligo_strand) or self.check_for_sticking(oligo_strand) or check_hamming(
                    self.contents["Strand"], oligo_strand):
                return False
            self.contents['Strand'].append(str(oligo_strand))
            self.contents['Numbers'].append(str_to_num(oligo_strand))
            self.contents['Length'].append(len(oligo_strand))
            self.update_energy_matrix()
            self.update_strand_energy()
            return True

    def check_for_same_strand(self, oligo_strand):
        """
        Checks for same existing strand in database
        :param oligo_strand: (Oligo/oligo object)
        :return: bool - True is same strand detected, False is no same strand detected
        """
        if isinstance(oligo_strand, Oligo):  # Justin's Oligo objects (numbers)
            oligo_sequence = oligo_strand.sequence
        elif isinstance(oligo_strand, oligo):  # Ariel's oligo objects (biopython)
            oligo_sequence = str(oligo_strand.sequence)
        elif isinstance(oligo_strand, Seq):
            oligo_sequence = str(oligo_strand)
        for each_sequence in self.contents['Strand']:
            if each_sequence == oligo_sequence:
                return True
        return False

    def check_for_sticking(self, oligo_strand):
        """
        Checks for strand sticking within existing database as defined by complement_filter_size
        :param oligo_strand: (Oligo/oligo object)
        :return: bool - True is sticking strand detected, False is no sticking strand detected
        """
        if isinstance(oligo_strand, Oligo):  # Justin's Oligo objects (numbers)
            oligo_sequence = oligo_strand.sequence
        elif isinstance(oligo_strand, oligo):  # Ariel's oligo objects (biopython)
            oligo_sequence = str(oligo_strand.sequence)
        elif isinstance(oligo_strand, Seq):
            oligo_sequence = str(oligo_strand)
        try:
            for x in range(0, len(oligo_sequence)):
                inserted_strand_sequences = []
                for y in range(x, x + complement_filter_size):
                    inserted_strand_sequences.append(oligo_sequence[y])

                inserted_strand_sequences.reverse()
                reverse_complement = ''.join(Oligo.complement(inserted_strand_sequences))

                for each in self.database['Strand']:
                    if reverse_complement in each:
                        return True
            return False
        except IndexError:
            return False

    def reorder_by_length(self):
        """
        Reorders database based on length; sort_condition defines ascending or descending order
        :return: None
        """
        if sort_condition:
            self.database.sort_values('Length', ascending=False, inplace=True)
        else:
            self.database.sort_values('Length', ascending=True, inplace=True)

    def export_to_csv(self):
        """
        Exports database to .csv file
        :return: None
        """
        self.database.to_csv('Oligo_database.csv', index=False)

# if __name__ == "__main__":
#     database = Database()
# 
#     for i in range(6, 16):  # Justin's Oligo
#         strand = Oligo(10, True, True, True, False)
#         database.database_insert(strand)
# 
#     strand2 = oligo(10)  # Ariel's oligo
#     database.database_insert(strand2)
#     
#     database.update_energy_matrix()
# 
#     print(database)
#     print("\nThis database has " + str(database.size) + " oligonucleotides.")
#     database.export_to_csv()
