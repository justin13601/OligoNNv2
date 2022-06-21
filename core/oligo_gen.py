"""
Generates oligo sequences as BioPython sequence objects
- added .numbers attribute for consistency as well as str_to_nums funct adapted from Justin
"""

from Bio.Seq import Seq, reverse_complement
from Bio.SeqUtils import GC 
from Bio.SeqUtils import MeltingTemp as mt
import random
from OligoNN.core import nupack_wrap as nupack

length = 10 #nt of oligos generated
gc_percent = 40 #fgc percent
gc_accuracy = 5 #fgc +/- x percent
rep_len = 3 #min repeat length
pin_len = 6 #min palindrome length
mp_method = "Wallace" #melting point function to use: Wallace, GC, or NN
mp_gen = 200 #how many sequences to test to get melting point conditions
mp_accuracy = 1 #fmp +/- x temp
"""dictionary whose key is the length of a dna sequence, and its values is a tuple
containing the median and the standard deviation of binding energy for sequence of that length
values obtained from generating 10,000 strands and calculating their median and stdev energies.
"""
dna_median_energy = {
        6  : (-7.45135141908, 0.45939899710),
        7  : (-8.88135141908, 0.71063569771),
        8  : (-10.2213514191, 0.49569376251),
        9  : (-11.6713514191, 0.73127004127),
        10 : (-13.0713514191, 0.54781790184),
        11 : (-14.4763514191, 0.80621698512),
        12 : (-15.8813514191, 1.07738262915),
        13 : (-17.3113514191, 0.81228370511),
        14 : (-18.6013514191, 1.08627160858),
        15 : (-20.1213514191, 0.83970005591),
        16 : (-21.5213514191, 1.13527658375),
        17 : (-22.9013514191, 1.38393291279),
        18 : (-24.2813514191, 1.13470857364),
        19 : (-25.6663514191, 1.37779371983),
        20 : (-27.0813514191, 1.17127609452),
        21 : (-28.4313514191, 1.45233172395),
        22 : (-29.9063514191, 1.66054894053),
        23 : (-31.2413514191, 1.48392673778),
        24 : (-32.7513514191, 1.72627418807),
        25 : (-34.0613514191, 1.48566292647),
        26 : (-35.4363514191, 1.82423875016),
        27 : (-36.9113514191, 2.02626664336),
        28 : (-38.3363514191, 1.70276119418),
        29 : (-39.6913514191, 2.04216766772),
        30 : (-41.2263514191, 1.79164414366),
        31 : (-42.3913514190, 2.04311030431),
        32 : (-44.0313514190, 2.30575292441)
        }

class oligo:
    def __init__(self, length, fgc=True, frepeat=True, fhairpin=True, fenergy=True, fmp=False, fprint=False):
        self.length = length
        self.fgc = fgc
        self.frepeat = frepeat
        self.fhairpin = fhairpin
        self.fenergy = fenergy
        self.fmp = fmp
        self.sequence = oligo.generate(self)
        self.numbers = str_to_num(self.sequence)
        self.gccontent = GC(self.sequence)
        self.meltingpt = melt(self.sequence, mp_method)
        # generate new sequences until all filters return True
        # note if a filter is turned off it will always return True (see in methods)
        while oligo.fgc(self, gc_percent, gc_accuracy)==False or oligo.frepeat(self, rep_len)==False \
              or oligo.fhairpin(self, pin_len)==False or oligo.fenergy(self)==False or oligo.fmp(self)==False:
            
            self.sequence = oligo.generate(self)
            self.gccontent = GC(self.sequence)
            self.meltingpt = melt(self.sequence, mp_method)
        if fprint:
            print(self.sequence)
    
    def generate(self):
        # generate a BioPython sequence of desired length by randomly choosing bases
        # weights correspond to desired gccontent to increase prob of passing fgc
        new_seq = Seq("".join(random.choices(["A", "T", "C", "G"], weights=[100-gc_percent, 100-gc_percent, gc_percent, gc_percent], k=self.length)))
        return new_seq
    
    def fgc(self, gc_percent, gc_accuracy):
        if not self.fgc:
            return True
        # if gccontent is above or below the threshold, does not pass filter
        if self.gccontent < (gc_percent - gc_accuracy) or self.gccontent > (gc_percent + gc_accuracy):
            return False
        return True     
    
    def frepeat(self, rep_len):
        if not self.frepeat:
            return True 
        # if repeats of any base are greater than the min length, does not pass filter
        # only has to check if a repeat exists with length +1 of min (if "AAAA" exists then so does "AAA")
        for base in [x*(rep_len+1) for x in "ATCG"]:
            if base in self.sequence:
                return False
        return True
    
    def fhairpin(self, pin_len):
        if not self.fhairpin:
            return True
        # iterate through each seq w lenth of half min palindrome length
        # check if reverse complement is in the remainder of the sequence
        pin = pin_len//2
        for x in range((self.length-pin_len)+1):
            frag = self.sequence[x:x+pin]
            if frag.reverse_complement() in self.sequence[x+pin:]:
                return False
        return True
    
    def fenergy(self):
        if not self.fenergy:
            return True
        # following has been adapted from TF generator
        # check if the difference btw the energy of the sequence and the median energy is less than or equal to the stdev
        dna_struct = self.length*"("+"+"+ self.length*")"
        seq=str(self.sequence)
        comp=complement(seq)
        seq_energy = nupack.energy((seq, comp), dna_struct, material ='dna')
        energy_median, energy_stdev = dna_median_energy[self.length][0],dna_median_energy[self.length][1]
        if abs(seq_energy - energy_median) > energy_stdev:
            return False
        return True
    
    def fmp(self):
        if not self.fmp:
            return True
        # get average melting point to use as filter conditions
        # check if melting point is within acceptable range
        sum_mp=0
        for _ in range(mp_gen):
            sum_mp += melt(oligo(length).sequence, mp_method)
        avg_mp= sum_mp/mp_gen
        if self.meltingpt > (avg_mp + mp_accuracy) or self.meltingpt < (avg_mp - mp_accuracy):
            return False
        return True

def translate(oligo):
    # assumes oligo is sense strand (just replaces T with U)
    return oligo.sequence.transcribe()

def reverse_translate(rna):
    # turns rna into dna (replace U with T) returns sense strand
    return rna.back_transcribe()

def complement(oligo):
    if type(oligo) == oligo:
        return oligo.sequence.reverse_complement()
    elif type(oligo) == Seq:
        return oligo.reverse_complement()
    elif type(oligo) == str:
        return reverse_complement(oligo)

def sequence(oligo):
    # sequence function can be accessed by self.sequence attribute
    return oligo.sequence

def count_bases(oligo):
    # returns a dictionary with nt as key and frequence (%) as values
    d = {}
    for base in "ATGC":
        if type(oligo) == oligo:
            d[base] = oligo.sequence.count(base)/oligo.length
        elif type(oligo) == Seq or type(oligo) == str:
            d[base] = oligo.count(base)/len(oligo)
    return d

def str_to_num(oligo):
    # Converts seq obj to list of numbers
    oligo=list(oligo)
    str_list = []
    str_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, 'U': 5}
    num_list = list(str_to_num.keys())
    nt_list = list(str_to_num.values())
    for base in oligo:
        str_list.append(nt_list[num_list.index(base)])
    return str_list

def melt(oligo, method = "Wallace"):
    if method == "Wallace":
        return(mt.Tm_Wallace(oligo))
    elif method == "GC":
        return(mt.Tm_GC(oligo))
    elif method == "NN":
        return(mt.Tm_NN(oligo))
    else:
        raise ValueError("Please specify a valid method: Wallace, GC, or NN.")

def fhairpin(oligo, pin_len):
    # oligo can be Seq obj or str
    # returns True if no hairpins
    # iterate through each seq w lenth of half min palindrome length
    # check if reverse complement is in the remainder of the sequence
    pin = pin_len//2
    for x in range((len(oligo)-pin_len)+1):
        frag = oligo[x:x+pin]
        if complement(frag) in oligo[x+pin:]:
            return False
    return True




