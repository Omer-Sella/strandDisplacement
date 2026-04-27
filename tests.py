from utils import *
def test_Similarity():
    
    seq1 = 'AAAAAAAAAACAGACAAAAAAAAAAAAACTTTAAACAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    seq2 = 'CGAAGCGCGAACAGGCGCCCGCCATAGAAGTGGGCCGTATCCCGGCTCCAACAACGGTCATGGTCCCGAACGATAACGCTGAACGTGGTCATCGGTGCAG'
    seq3 = 'CTGCACCGATGACCACGTTCAGCGTTATCGTTCGGGACCATGACCGTTGTTGGAGCCGGGATACGGCCCACTTCTATGGCGGGCGCCTGTTCGCGCTTCG'
    #print(gibbsFreeEnergy(seq1, seq2))
    #print(gibbsFreeEnergy(seq1, seq3))
    #print(gibbsFreeEnergy(seq2, seq3))
    print(binds_too_strongly(seq1, seq2))
    print(binds_too_strongly(seq1, seq3))
    print(binds_too_strongly(seq2, seq3))

test_Similarity()
