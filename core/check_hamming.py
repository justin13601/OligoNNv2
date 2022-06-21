"""
Hamming distance function, implemented as a filter in oligo_database.
"""

min_ham = 2  # all strands must have a hamming distance of 2


def check_hamming(strands, oligo):
    """
    Checks if the new strand has a minimum hamming distance with all other strands.
    
    Parameters: strands - a list of dna sequences in the database, oligo - a dna sequence to check
    Returns: A boolean value, False if oligo passes the filter
    """
    oligo_seq = str(oligo)
    for strand in strands:
        # only checks strands of the same length
        if len(strand) == len(oligo_seq):
            hamming = sum(s1 != s2 for s1, s2 in zip(oligo_seq, strand))
            if hamming < min_ham:
                return True
    return False
