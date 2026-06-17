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

def validate_free_gibbs_energy_in_chunks(fasta_path):
    """
    Load DNA chunks from a FASTA file and calculate the minimum free Gibbs energy
    between any two different DNA strands.
    
    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file containing DNA chunks (saved with save_dna_chunks_to_fasta).
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'min_energy': Lowest free Gibbs energy found
        - 'strand1_idx': Index of first strand in minimum energy pair
        - 'strand2_idx': Index of second strand in minimum energy pair
        - 'strand1_seq': Sequence of first strand
        - 'strand2_seq': Sequence of second strand
        - 'total_comparisons': Total number of strand pairs compared
    """
    # Load DNA chunks from FASTA
    dna_chunks = load_dna_chunks_from_fasta(fasta_path)
    print(f"Loaded {len(dna_chunks)} DNA chunks from {fasta_path}")
    
    # Extract DNA sequences
    dna_sequences = [chunk['dna'] for chunk in dna_chunks]
    print(f"Extracted {len(dna_sequences)} DNA sequences")
    
    min_energy = float('inf')
    min_pair = None
    comparison_count = 0
    
    # Compare all pairs of different DNA strands
    for i in range(len(dna_sequences)):
        for j in range(i + 1, len(dna_sequences)):
            strand1 = dna_sequences[i]
            strand2 = dna_sequences[j]
            
            # Create complex with the two strands
            complexToBeChecked = Complex([
                Strand(string=strand1, name=f'strand_{i}'),
                Strand(string=strand2, name=f'strand_{j}')
            ])
            
            # Calculate free Gibbs energy
            energy = mfe(complexToBeChecked, model=MODEL_FOR_NUPACK)[0].energy
            comparison_count += 1
            
            # Track minimum energy
            if energy < min_energy:
                min_energy = energy
                min_pair = (i, j, strand1, strand2)
    
    result = {
        'min_energy': min_energy,
        'strand1_idx': min_pair[0] if min_pair else None,
        'strand2_idx': min_pair[1] if min_pair else None,
        'strand1_seq': min_pair[2] if min_pair else None,
        'strand2_seq': min_pair[3] if min_pair else None,
        'total_comparisons': comparison_count
    }
    
    print(f"\nFree Gibbs Energy Analysis Results:")
    print(f"Total strand pairs compared: {comparison_count}")
    print(f"Minimum free Gibbs energy: {min_energy:.3f} kcal/mol")
    print(f"Between strands {min_pair[0]} and {min_pair[1]}")
    print(f"Strand {min_pair[0]}: {min_pair[2][:50]}..." if len(min_pair[2]) > 50 else f"Strand {min_pair[0]}: {min_pair[2]}")
    print(f"Strand {min_pair[1]}: {min_pair[3][:50]}..." if len(min_pair[3]) > 50 else f"Strand {min_pair[1]}: {min_pair[3]}")
    
    return result

test_Similarity()
