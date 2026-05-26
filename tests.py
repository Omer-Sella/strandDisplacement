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

def test_deltaGCalculation():
    from nupack import mfe
    dnaLibrary                 = ['CGAAGCGCGAACAGGCGCCCGCCATAGAAGTGGGCCGTATCCCGGCTCCAACAACGGTCATGGTCCCGAACGATAACGCTGAACGTGGTCATCGGTGCAG', 'CTGCACCGATGACCACGTTCAGCGTTATCGTTCGGGACCATGACCGTTGTTGGAGCCGGGATACGGCCCACTTCTATGGCGGGCGCCTGTTCGCGCTTCG']
    dnaLibraryANDnextCandidate = ['CGAAGCGCGAACAGGCGCCCGCCATAGAAGTGGGCCGTATCCCGGCTCCAACAACGGTCATGGTCCCGAACGATAACGCTGAACGTGGTCATCGGTGCAG', 'CTGCACCGATGACCACGTTCAGCGTTATCGTTCGGGACCATGACCGTTGTTGGAGCCGGGATACGGCCCACTTCTATGGCGGGCGCCTGTTCGCGCTTCG', 'AAAAAAAAAACAGACAAAAAAAAAAAAACTTTAAACAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']
    dnaLibraryAsComplex = Complex([Strand(string=s, name=f'strand_{s}') for s in dnaLibraryANDnextCandidateDNA])
    model = Model(material='dna', celsius=37, sodium=0.05)
    stabilityScore = mfe(dnaLibraryAsComplex, model=model)[0].energy
    print('OK')

def test_scrambling_and_unscrambling():
    png_path = 'imperialBlue_with_plte_2_colours.png'


    # 2. Set the limit on the number of bits per chunk. This will affect how many bases are in every DNA strand.
    K_BITS = 200  # Placeholder constant; you can change this later. The number of bases per strand == K_BITS / 2

    dna_chunks = png_to_dna_chunks(png_path, K_BITS)
    
    for c in dna_chunks:
        localRandom = np.random.RandomState(c['seed'])
        scramble1 = "".join(localRandom.choice(['A'], size = len(c['dna'])))
        print(f"Scramble generated from seed:")
        print(scramble1)
        print(f"Scramble logged in the chunk:")
        print(c['scramble'])
        assert(c['scramble'] ==  scramble1)
        assert(c['dna'] == addScrambleToDnaSequence(c['unscrambled_dna'], scramble1))

test_Similarity()
