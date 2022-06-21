import numpy as np
import time

# define
repeat_filter = 4  # Greater than 3 repeats
first_pass_filter = 0.4  # 40%
gc_filter_tolerance = 0.1  # 10%
palindromic_filter = 6  # Sequences greater or equal to 6 n.t.


# Nucleotides
# A = 1     T = 2       G = 3       C = 4       U = 5


class Oligo():
    def __init__(self, length, fgc, frepeat, fhairpin, fprint):  # Filter Toggles
        self.numbers = self.generate(length, fgc, frepeat, fhairpin)
        self.length = length
        self.gccontent = self.find_gc(self.numbers, length)

        if fprint:  # Print Property
            print("Generated Oligonucleotide:")
            print(self.sequence(self.numbers))

    def __str__(self):
        return self.sequence(self.numbers)

    @staticmethod
    def num_to_str(numbers):  # Converts array of integers to array of letters
        nt_list = []
        num_to_str = {1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'U'}
        values_nt = list(num_to_str.keys())
        keys_num = list(num_to_str.values())
        for num in numbers:
            nt_list.append(keys_num[values_nt.index(num)])
        return nt_list

    @staticmethod
    def str_to_num(nts):  # Converts array of letters to array of integers
        str_list = []
        str_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, 'U': 5}
        num_list = list(str_to_num.keys())
        nt_list = list(str_to_num.values())
        for nt in nts:
            str_list.append(nt_list[num_list.index(nt)])
        return str_list

    @staticmethod
    def find_gc(numbers, length):  # Finds percentage of GC
        gc = numbers.count(3) + numbers.count(4)
        return round(gc / length, 2)

    def generate(self, length, fgc, frepeat, fhairpin):
        numbers = list(np.random.randint(low=1, high=5, size=6))  # Generates random 6 nt strand
        while self.fgc(numbers, len(numbers)) is False or self.fhairpin(numbers) is False or self.frepeat(
                numbers) is False:  # Tests if sequence fits initial conditions
            numbers = list(np.random.randint(low=1, high=5, size=6))
        for x in range(0, length - 6):  # Inserts until length is reached
            numbers = self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
        return numbers

    def generate_helper(self, numbers, length, fgc, frepeat, fhairpin):
        if fgc:
            if frepeat:
                if fhairpin:
                    # if len(numbers) == 6:
                    # print("All 3 Filters Applied")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
                else:
                    # if len(numbers) == 6:
                    # print("Filters 1 & 2")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
            else:
                if fhairpin:
                    # if len(numbers) == 6:
                    # print("Filters 1 & 3")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
                else:
                    # if len(numbers) == 6:
                    # print("Filter 1 Only")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
        else:
            if frepeat:
                if fhairpin:
                    # if len(numbers) == 6:
                    # print("Filters 2 and 3")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
                else:
                    # if len(numbers) == 6:
                    # print("Filter 2 Only")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
            else:
                if fhairpin:
                    # if len(numbers) == 6:
                    # print("Filter 3 Only")
                    numbers.append(np.random.randint(1, 5))
                    if self.frepeat(numbers) and self.fhairpin(numbers):
                        # print("Filter(s) Passed")
                        # print(self.num_to_str(numbers))
                        return numbers
                    else:
                        # print("Filter(s) Failed")
                        # print(self.num_to_str(numbers))
                        numbers.pop()
                        self.generate_helper(numbers, length, fgc, frepeat, fhairpin)
                else:
                    # if len(numbers) == 6:
                    # print("No Filters Applied")
                    numbers = list(np.random.randint(low=1, high=5, size=length - 6))
                    return numbers
        return numbers

    def fgc(self, numbers, length):  # GC content filter; True if passes, False if not
        if first_pass_filter + gc_filter_tolerance > self.find_gc(numbers,
                                                                  length) > first_pass_filter - gc_filter_tolerance:
            return True
        return False

    def frepeat(self, numbers):  # Repeat filter; True if passes, False if not
        try:
            final_nt_index = len(numbers)
            test_for_repeats = []
            if final_nt_index == 6:  # Checks initial 6 nt strand for repeats
                return not self.frepeat_helper(numbers)
            else:
                for j in range(final_nt_index - repeat_filter, final_nt_index):
                    test_for_repeats.append(numbers[j])
            if len(test_for_repeats) > 0:
                return not (all(elem == test_for_repeats[0] for elem in test_for_repeats))
        except IndexError:
            return True

    @staticmethod
    def frepeat_helper(numbers):  # Checks initial 6 nt strand for repeats; True if repeats exists, False if not
        try:
            for i in range(0, len(numbers)):
                repeats = [numbers[i]]
                for j in range(1, repeat_filter):
                    repeats.append(numbers[i + j])
                if all(elem == repeats[0] for elem in repeats):
                    return True
        except IndexError:
            return False

    def fhairpin(self, numbers):  # Palindrome Filter; True if passes, False if not
        try:
            final_nt_index = len(numbers)
            if final_nt_index == 6:  # Checks initial 6 nt strand for palindromes
                return not self.fhairpin_helper(numbers)
            else:
                return not self.fhairpin_helper_2(numbers, 0)
        except IndexError:
            return True

    def fhairpin_helper(self, numbers):  # Checks initial 6 nt strand for palindromes; True if hairpin exists, False if not
        try:
            for x in range(len(numbers), 0, -1):
                j = x - palindromic_filter
                while j >= 0:
                    palindrome = []
                    for k in range(j, x):
                        palindrome.append(numbers[k])
                    j -= 1
                    if self.checkPalindrome(palindrome):
                        return True
            return False
        except IndexError:
            return False

    def fhairpin_helper_2(self, numbers, index):  # Recursively checks for hairpins; True if exists, False if not
        try:
            reversed_complement = self.complement(reversed(numbers[-palindromic_filter:]))
            sequence_1 = []
            for x in range(0, palindromic_filter):
                sequence_1.append(numbers[index + x])
            if sequence_1 == reversed_complement:
                return True
            else:
                index += 1
                return self.fhairpin_helper_2(numbers, index)
        except IndexError:
            return False

    @staticmethod
    def sequence(numbers):  # Converts strand to letters
        nt_list = []
        num_to_str = {1: 'A', 2: 'T', 3: 'G', 4: 'C', 5: 'U'}
        values_nt = list(num_to_str.keys())
        keys_num = list(num_to_str.values())
        for num in numbers:
            nt_list.append(keys_num[values_nt.index(num)])
        return ''.join(nt_list)

    @staticmethod
    def translate(numbers):  # Transcribes to RNA when oligonucleotide is coding strand
        result = []
        rna_pairs = {1: 1, 4: 4, 3: 3, 2: 5}
        for base in numbers:
            result.append(rna_pairs[base])
        return result

    def checkPalindrome(self, numbers):  # Checks if sequence is Palindrome
        complementary_numbers = self.complement(numbers)
        reversed_numbers = numbers[::-1]
        if reversed_numbers == complementary_numbers:
            return True
        else:
            return False

    @staticmethod
    def complement(numbers):  # Finds complement strand
        result = []
        pairs = {1: 2, 4: 3, 3: 4, 2: 1, "A": "T", "T": "A", "C": "G", "G": "C"}
        for base in numbers:
            result.append(pairs[base])
        return result


# Testing
#if __name__ == "__main__":
    #start_time = time.time()

    #test_strand = Oligo(30, True, True, True, True)

    # print("Generating 50k sequences...")
    # for i in range(0, 50000):
    #    test_strand_2 = Oligo(15, True, True, True, False)

    # test_strand_3 = Oligo(15, True, True, True, False)
    # print("\nSequenced Oligonucleotide & Translated Oligonucleotide:")
    # print(test_strand_3)
    # print(sequence(translate(test_strand_3.numbers)))

    #print("\n--- %s seconds ---" % (time.time() - start_time))
    # 32 seconds on avg over 10 runs
    # Time is based on computer performance and printing?
