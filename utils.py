import numpy as np
# Omer@Ilaria - is this a good enough metric ?
from Bio import Align, pairwise2
from Bio.Seq import Seq as seq
import numpy as np
import pandas as pd
import zlib
import copy

# Declaring / instantiating an alignment function
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
from pathlib import Path
from PIL import Image
from IPython.display import display
#AAAAAAAAAAAAGGGGGGGGGGGGGGGG
#            CCCCCCCCCCCCCCAAAAAAAAAAAAAAAAA
# See potential penalties to be defined here: https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html#chapter-pairwise
aligner.match_score = 1 # Each match counts as 1
aligner.mismatch_score = -1
#aligner.gap_score = 0 # No penalty or gain for gaps
aligner.open_gap_score = -1
aligner.extend_gap_score = -1

def basesToBytes(inputBases):
  if type(inputBases) != str:
    raise TypeError("Input must be a string")
  outputBitStream = ''
  for i in range(len(inputBases)):
    if inputBases[i] == "A":
        twoBits = '00'
    elif inputBases[i] == "C":
        twoBits = '01'
    elif inputBases[i] == "T":
        twoBits = '11'
    elif inputBases[i] == "G":
        twoBits = '10'
    else:
        raise ValueError(f"Invalid base encountered {inputBases[i]}")

    outputBitStream = outputBitStream + twoBits
  return outputBitStream

def byteToBases(inputByte, numOfBits = 8):
    if type(inputByte) != str:
        bitStream = (bin(int(inputByte.hex(), 16))[2:]).zfill(numOfBits)
    else:
        bitStream = inputByte
    fourNuceleotides = ""
    for i in range(numOfBits // 2):
        twoBits = bitStream[2 * i : 2 * i + 2]
        #print(twoBits)
        if twoBits == "00":
            nuceleotide = 'A'
        elif twoBits == "01":
            nuceleotide = 'C'
        elif twoBits == "11":
            nuceleotide = 'T'
        elif twoBits == "10":
            nuceleotide = 'G'
        else:
            print(f"Bad value {twoBits}")
            raise ValueError(f"Invalid string encountered {twoBits}")
        fourNuceleotides = fourNuceleotides + nuceleotide
    return fourNuceleotides

def DNAAddition(lhs, rhs):
    lhsBinary = basesToBytes(lhs)
    rhsBinary = basesToBytes(rhs)
    twoBits = str((int(lhsBinary[0]) + int(rhsBinary[0])) %2) + str((int(lhsBinary[1]) + int(rhsBinary[1])) %2)
    # Convert 2-bit value back to one DNA base.
    result = byteToBases(twoBits, len(twoBits))
    return result

def addScrambleToDnaSequence(dnaSequence, scramble):
  
  if len(dnaSequence) != len(scramble):
    raise
  else:
    return ''.join([DNAAddition(dnaSequence[i], scramble[i]) for i in range(len(dnaSequence))])


def checkHomopolymer(dnaString):
    if dnaString == '' or len(dnaString) < 1:
      return 0
    else:
      c = dnaString[0]
      maxRun = 1
      currentRun = 1
      for i in range(1,len(dnaString)):
        if dnaString[i] == c:
          currentRun = currentRun + 1
          if currentRun > maxRun:
            maxRun = currentRun
        else:
          c = dnaString[i]
          currentRun = 1
      return maxRun

#Let's check that DNAAddition of a scramble is an order 2 operator (so adding the same scramble twice is the identity):
#dna1 = "AAATTTTTTT"
#print(basesToBytes(dna1))
#scramble = "TTTAAAAAAA"
#scrambledDNA = addScrambleToDnaSequence(dna1, scramble)
#print(f"{dna1} + {scramble} == {scrambledDNA}")
#print(f"{scrambledDNA} + {scramble} == {addScrambleToDnaSequence(scrambledDNA, scramble)}")


def checkGC(dnaString):
  return np.sum([1 for n in dnaString if (n=='C' or n=='G')]) / len(dnaString)

def commitSingleBinaryToDnaLibrary(binaryToBeCommitted, dnaLibrary):
    K = 50
    MAX_HOMOPOLYMER = 8
    GC_LOW = 0.2
    GC_HIGH = 0.8
    # We're going to use seeds between 0 and maximumSeed-1 (you can replace numberOfBits with something different, but the point is store this either at the beginning of the sequence, or offline, so it should be as small as possible)
    maximumSeed = 256
    #for s in textLibrary[1:]:
    similarityScore = K + 1
    unscrambledDNA = byteToBases(binaryToBeCommitted, len(binaryToBeCommitted))
    nextCandidateDNA = copy.copy(unscrambledDNA)
    localRandom = np.random.RandomState(0)
    scramble1 = "".join(localRandom.choice(['A'], size = len(unscrambledDNA)))
    #scrambledDNA = seq(addScrambleToDnaSequence(unscrambledDNA, scramble1))
    #print(f"Attempting to commit {unscrambledDNA}")
    i = 0
    bestSimilarityScore = np.inf
    bestSeed = 0
    commitable = False
    
    while (i< maximumSeed) and not commitable:
        # Get the similarity scores between the candidate and every DNA already commited to the library
        if len(dnaLibrary) == 0:
            # There are no dna strands in the library, so there is no similarity problem, only GC and homopolymer
            similarityScore = 0;
        else:
            scores = [aligner.score(seq(nextCandidateDNA), targetSeq) for targetSeq in dnaLibrary]
            #For debug purposes we can look at the alignments found, but this takes a very long time.
            #local_alignments = aligner.align(seq(nextCandidateDNA), dnaLibrary[0])
            #[print(l) for l in local_alignments]
            #the highest similarity score is the worst score
            similarityScore = max(scores)
        commitable = (similarityScore < K) and (checkGC(nextCandidateDNA) > GC_LOW) and (checkGC(nextCandidateDNA) < GC_HIGH) and (checkHomopolymer(nextCandidateDNA) < MAX_HOMOPOLYMER)
        #print(f"Attempting to commit using seed == {i}")
        if commitable: 
            scrambledDNA = seq(nextCandidateDNA)
        # The while loop will break after this
        else:
            # Prepare next attempt of DNA sequence
            i = i + 1
            # There is no real need to instantiat the rng here, just set it to a seed, but I figured this is easier to understand
            localRandom = np.random.RandomState(i)
            # Now we have a choice, either permute or one-time-pad the DNA sequence. Right now I am only supporting permutation
            scramble1 = "".join(localRandom.choice(['A' ,'C' ,'T' ,'G'], size = len(nextCandidateDNA)))
            #nextCandidateDNABeforeScrambling = nextCandidateDNA
            nextCandidateDNA = addScrambleToDnaSequence(unscrambledDNA, scramble1)
            # Data validation step: make sure that the sequence can be unscrambled:
            localRandom = np.random.RandomState(i)
            scramble2 = "".join(localRandom.choice(['A' ,'C' ,'T' ,'G'], size = len(nextCandidateDNA)))
            assert(scramble1 == scramble2)
            checkUnscrambledDNA = addScrambleToDnaSequence(nextCandidateDNA, scramble2)
            assert(checkUnscrambledDNA == unscrambledDNA)
        if i == maximumSeed:
            print(f"Failed to find a seed to meet all constraints for the sequence {binaryToBeCommitted}")
            raise Exception("Failed to encode the binary library into a DNA library with sufficient dissimilarity. Consider changing the number of possible seeds, or allowing for more similarity.")

    return scrambledDNA, unscrambledDNA, i, scramble1, similarityScore

def binaryLibraryToDnaLibrary(binaryLibrary):
    dnaLibrary = [] #We start with an empty library seq0, seq0.reverse_complement()]
    # Now we're finally going to populate the dna librtary
    # This is where we set the similarity score
    
    # By convention, if there was no need for a permutation / scrambling then we assign seed == 0
    seeds = []#[0]
    scrambles = []
    # permute == True then we try making the DNA strand different by permutation.
    # If permute == False, then we try making the DNA different by bitwise XOR with a random pattern.
    permute = False
    scores = []
    sourceDNAUnscrambled = []

    # Now we're going to go over the binaery elements one by one, try the simple map, if it works - great ! if not, try to scramble/ shuffle / alter the DNA data so that it does work.
    # What we're going to end up with, is:
    # dnaLibrary - a library of dna strands, sufficiently dissimilar.
    # seeds - a list of seeds that we used to shuffle the original DNA, in order to get dissimilar dna strands.
    # scores - a list of scores that the sligner found.
    
    for j, binaryCandidate in zip(range(len(binaryLibrary)),  binaryLibrary): 
        # First we try using the simple mapping in byteToBases:
        scrambledDna, unscrambledDna, seed, scramble, score = commitSingleBinaryToDnaLibrary(binaryCandidate, dnaLibrary = dnaLibrary)
        dnaLibrary.append(scrambledDna)
        dnaLibrary.append(scrambledDna.reverse_complement())
        sourceDNAUnscrambled.append(unscrambledDna)
        seeds.append(seed)
        scrambles.append(scramble)
        scores.append(score)
        
        
    return dnaLibrary, seeds, scrambles, scores, sourceDNAUnscrambled

def dnaLibraryToBinaryUnscrambler(dnaLibrary, seeds):
    if not (len(dnaLibrary) == len(seeds)):
        raise ValueError("Number of dna strands must match the number of seeds")
    binaryLibrary = []
    for strand, seed in zip (dnaLibrary, seeds):
        if seed == 0:
            binaryLibrary.append(basesToBytes(strand))
        else:
            localRandom = np.random.RandomState(seed)
            scramble2 = "".join(localRandom.choice(['A' ,'C' ,'T' ,'G'], size = len(strand)))
            unscrambledDNA = addScrambleToDnaSequence(strand, scramble2)
            binaryLibrary.append(basesToBytes(unscrambledDNA))
    return binaryLibrary
    
def bytes_to_bitstring(data: bytes) -> str:
    return ''.join(f'{b:08b}' for b in data)


def split_bits(bits: str, k_bits: int):
    if k_bits <= 0:
        raise ValueError('k_bits must be a positive integer')
    return [bits[i:i + k_bits] for i in range(0, len(bits), k_bits)] if bits else []

def pad_bits(bits: str, k_bits: int):
    """
    Pad a bitstring with trailing zeros to exactly k_bits.
    Returns (padded_bits, valid_bits_len).
    """
    if len(bits) > k_bits:
        raise ValueError(f"Cannot pad: bitstring length {len(bits)} exceeds k_bits={k_bits}.")
    return bits + ('0' * (k_bits - len(bits))), len(bits)


def parse_png_chunks(png_bytes: bytes):
    """
    Parse a PNG byte stream into logical PNG chunks.

    Each returned item contains raw bytes for:
    - length (4 bytes)
    - type (4 bytes)
    - data (variable length)
    - crc/checksum (4 bytes)

    The function first validates the PNG signature, then walks chunk-by-chunk
    until it reaches IEND (or the input ends), raising an error on truncation.
    """
    png_signature = b'\x89PNG\r\n\x1a\n'
    if len(png_bytes) < 8 or png_bytes[:8] != png_signature:
        raise ValueError('Input is not a valid PNG file (invalid signature).')

    offset = 8
    parsed = []
    while offset + 12 <= len(png_bytes):
        length = int.from_bytes(png_bytes[offset:offset + 4], 'big')
        chunk_type = png_bytes[offset + 4:offset + 8]
        data_start = offset + 8
        data_end = data_start + length
        crc_start = data_end
        crc_end = crc_start + 4

        if crc_end > len(png_bytes):
            raise ValueError('Malformed PNG: truncated chunk encountered.')

        parsed.append({
            'type': chunk_type.decode('ascii', errors='replace'),
            'length_bytes': png_bytes[offset:offset + 4],
            'type_bytes': chunk_type,
            'data_bytes': png_bytes[data_start:data_end],
            'crc_bytes': png_bytes[crc_start:crc_end],
        })

        offset = crc_end
        if chunk_type == b'IEND':
            break

    return parsed


def build_png_bit_chunks(png_bytes: bytes, k_bits: int):
    """
    Convert parsed PNG chunks into ordered bit chunks of size <= k_bits.

    For each PNG chunk, this function emits subchunks for:
    - header (length + type)
    - payload data (or 'palette' for PLTE)
    - checksum (CRC)

    Special handling:
    - PLTE payload is tagged as field='palette' so it stays separate.
    - CRC is always tagged as field='checksum' so it remains dedicated.

    Returns a flat list of dictionaries, each carrying chunk metadata and a
    bitstring slice ready for DNA mapping.
    """
    parsed_chunks = parse_png_chunks(png_bytes)
    all_chunks = []

    for i, ch in enumerate(parsed_chunks):
        chunk_type = ch['type']

        # Normal fields (length + type) are chunked by K bits.
        header_bits = bytes_to_bitstring(ch['length_bytes'] + ch['type_bytes'])
        for j, sub in enumerate(split_bits(header_bits, k_bits)):
            padded_bits, valid_bits = pad_bits(sub, k_bits)
            all_chunks.append({
                'chunk_index': i,
                'png_chunk_type': chunk_type,
                'field': 'header',
                'subchunk_index': j,
                'bits': padded_bits,
                'valid_bits': valid_bits,
            })

        data_bits = bytes_to_bitstring(ch['data_bytes'])
        if chunk_type == 'PLTE':
            # Palette is a dedicated chunk group.
            for j, sub in enumerate(split_bits(data_bits, k_bits)):
                padded_bits, valid_bits = pad_bits(sub, k_bits)
                all_chunks.append({
                    'chunk_index': i,
                    'png_chunk_type': chunk_type,
                    'field': 'palette',
                    'subchunk_index': j,
                    'bits': padded_bits,
                    'valid_bits': valid_bits,
                })
        else:
            for j, sub in enumerate(split_bits(data_bits, k_bits)):
                padded_bits, valid_bits = pad_bits(sub, k_bits)
                all_chunks.append({
                    'chunk_index': i,
                    'png_chunk_type': chunk_type,
                    'field': 'data',
                    'subchunk_index': j,
                    'bits': padded_bits,
                    'valid_bits': valid_bits,
                })

        # Checksum is always a dedicated chunk group.
        crc_bits = bytes_to_bitstring(ch['crc_bytes'])
        for j, sub in enumerate(split_bits(crc_bits, k_bits)):
            padded_bits, valid_bits = pad_bits(sub, k_bits)
            all_chunks.append({
                'chunk_index': i,
                'png_chunk_type': chunk_type,
                'field': 'checksum',
                'subchunk_index': j,
                'bits': padded_bits,
                'valid_bits': valid_bits,
            })

    return all_chunks


def _encode_bits_to_dna_with_seed(bitstring: str, dna_library_context):
    """
    Encode one bitstring and return (dna, seed).

    If dna_library_context is provided, commitSingleBinaryToDnaLibrary is used
    so encoding considers existing DNA strands.
    """
    if dna_library_context is None:
        dna_library_context = []
    scrambled_dna, _, seed, _, _ = commitSingleBinaryToDnaLibrary(
        bitstring,
        dnaLibrary=dna_library_context
    )
    return str(scrambled_dna), int(seed)


def png_to_dna_chunks(png_path: str, k_bits: int):
    png_bytes = Path(png_path).read_bytes()
    bit_chunks = build_png_bit_chunks(png_bytes, k_bits)

    dna_chunks = []
    dna_library_context = []
    for c in bit_chunks:
        dna_strand, _, seed, _, _ = commitSingleBinaryToDnaLibrary(
            c['bits'],
            dnaLibrary=dna_library_context
        )
        dna_chunks.append({
            **c,
            'dna': str(dna_strand),
            'seed': int(seed),
        })
        # Keep both strand and reverse-complement in context, matching library rules.
        dna_library_context.append(seq(str(dna_strand)))
        dna_library_context.append(seq(str(dna_strand)).reverse_complement())

    return dna_chunks



def bitstring_to_bytes(bits: str) -> bytes:
    if len(bits) % 8 != 0:
        raise ValueError(f'Bitstring length must be a multiple of 8, got {len(bits)}.')
    return bytes(int(bits[i:i + 8], 2) for i in range(0, len(bits), 8))


def reconstruct_png_bytes_from_dna_chunks(dna_chunks):
    png_signature = b'\x89PNG\r\n\x1a\n'

    # Group by original PNG chunk index while preserving field separation.
    grouped = {}
    for item in dna_chunks:
        idx = item['chunk_index']
        grouped.setdefault(idx, {'header': [], 'data': [], 'palette': [], 'checksum': []})

        seed = int(item.get('seed', 0))
        recovered_bits = dnaLibraryToBinaryUnscrambler([item['dna']], [seed])[0]
        expected_bits_len = len(item.get('bits', recovered_bits))
        if len(recovered_bits) != expected_bits_len:
            raise ValueError(
                f"Recovered bits length mismatch at chunk {idx}, field {item['field']} "
                f"(expected {expected_bits_len}, got {len(recovered_bits)})."
            )

        valid_bits_len = item.get('valid_bits', len(recovered_bits))
        if valid_bits_len > len(recovered_bits):
            raise ValueError(
                f"valid_bits ({valid_bits_len}) exceeds decoded bit length ({len(recovered_bits)}) "
                f"at chunk {idx}, field {item['field']}."
            )
        grouped[idx][item['field']].append((item['subchunk_index'], recovered_bits[:valid_bits_len]))

    assembled = bytearray(png_signature)

    for idx in sorted(grouped.keys()):
        g = grouped[idx]

        header_bits = ''.join(bits for _, bits in sorted(g['header'], key=lambda x: x[0]))
        checksum_bits = ''.join(bits for _, bits in sorted(g['checksum'], key=lambda x: x[0]))

        # PLTE chunk data is kept in 'palette'; other chunk payloads are in 'data'.
        payload_field = 'palette' if g['palette'] else 'data'
        payload_bits = ''.join(bits for _, bits in sorted(g[payload_field], key=lambda x: x[0]))

        chunk_bytes = bitstring_to_bytes(header_bits) + bitstring_to_bytes(payload_bits) + bitstring_to_bytes(checksum_bits)
        assembled.extend(chunk_bytes)

    return bytes(assembled)


def reconstruct_and_display_png(dna_chunks, output_png_path='reconstructed.png'):
    png_bytes = reconstruct_png_bytes_from_dna_chunks(dna_chunks)
    out_path = Path(output_png_path)
    out_path.write_bytes(png_bytes)

    print(f'Reconstructed PNG written to: {out_path.resolve()}')
    img = Image.open(out_path)
    display(img)
    return out_path


def extract_plte_table(png_path: str, display_table: bool = True):
    """
    Extract PLTE payload from a PNG and return it as a palette table.

    Returns a pandas DataFrame with columns:
    - index: palette entry index
    - R, G, B: color channel values (0-255)
    """
    png_bytes = Path(png_path).read_bytes()
    parsed_chunks = parse_png_chunks(png_bytes)

    plte_chunk = next((ch for ch in parsed_chunks if ch['type'] == 'PLTE'), None)
    if plte_chunk is None:
        raise ValueError(f'No PLTE chunk found in PNG: {png_path}')

    plte_data = plte_chunk['data_bytes']
    if len(plte_data) % 3 != 0:
        raise ValueError(f'Malformed PLTE length: {len(plte_data)} is not divisible by 3.')

    rows = []
    for idx in range(len(plte_data) // 3):
        base = idx * 3
        rows.append({
            'index': idx,
            'R': plte_data[base],
            'G': plte_data[base + 1],
            'B': plte_data[base + 2],
        })

    table = pd.DataFrame(rows, columns=['index', 'R', 'G', 'B'])
    if display_table:
        display(table)
    return table


def update_plte_bits_in_dna_chunks(dna_chunks, new_plte_bits=None, bit_modifier=None, in_place=False):
    """
    Update PLTE payload bits in dna_chunks, fix CRC bits, and re-encode DNA.

    Parameters
    ----------
    dna_chunks : list[dict]
        Output of png_to_dna_chunks().
    new_plte_bits : str | None
        Full replacement bitstring for the PLTE payload.
    bit_modifier : callable | None
        Function that receives current PLTE payload bits and returns updated bits.
    in_place : bool
        If True, mutate the input list. Otherwise return a modified copy.

    Notes
    -----
    Exactly one of new_plte_bits / bit_modifier must be provided.
    The updated PLTE bitstring must preserve length and be byte-aligned.
    """
    if (new_plte_bits is None) == (bit_modifier is None):
        raise ValueError("Provide exactly one of new_plte_bits or bit_modifier.")

    target = dna_chunks if in_place else [dict(item) for item in dna_chunks]

    # Locate PLTE chunk index.
    plte_indices = sorted({c['chunk_index'] for c in target if c.get('png_chunk_type') == 'PLTE'})
    if not plte_indices:
        raise ValueError("No PLTE chunk found in dna_chunks.")
    if len(plte_indices) > 1:
        raise ValueError("Multiple PLTE chunks found; expected a single PLTE chunk.")
    plte_idx = plte_indices[0]

    plte_palette_items = [c for c in target if c['chunk_index'] == plte_idx and c['field'] == 'palette']
    checksum_items = [c for c in target if c['chunk_index'] == plte_idx and c['field'] == 'checksum']
    if not plte_palette_items or not checksum_items:
        raise ValueError("PLTE chunk is missing palette or checksum fields.")

    plte_palette_items.sort(key=lambda x: x['subchunk_index'])
    checksum_items.sort(key=lambda x: x['subchunk_index'])
    changed_item_ids = {id(item) for item in plte_palette_items + checksum_items}
    unchanged_context = [seq(c['dna']) for c in target if id(c) not in changed_item_ids]

    original_plte_bits = ''.join(c['bits'][:c.get('valid_bits', len(c['bits']))] for c in plte_palette_items)
    if new_plte_bits is None:
        updated_plte_bits = bit_modifier(original_plte_bits)
    else:
        updated_plte_bits = new_plte_bits

    if not isinstance(updated_plte_bits, str):
        raise TypeError("Updated PLTE bits must be a string.")
    if any(b not in '01' for b in updated_plte_bits):
        raise ValueError("Updated PLTE bits must contain only '0' and '1'.")
    if len(updated_plte_bits) != len(original_plte_bits):
        raise ValueError(
            f"Updated PLTE bits length mismatch: expected {len(original_plte_bits)}, got {len(updated_plte_bits)}."
        )
    if len(updated_plte_bits) % 8 != 0:
        raise ValueError("Updated PLTE bits must be byte-aligned (multiple of 8).")

    # Write modified PLTE bits back using the original subchunk sizes.
    cursor = 0
    for item in plte_palette_items:
        n_total = len(item['bits'])
        n_valid = item.get('valid_bits', n_total)
        sub_bits_valid = updated_plte_bits[cursor:cursor + n_valid]
        if len(sub_bits_valid) != n_valid:
            raise ValueError("Failed to partition updated PLTE bits into original subchunks.")
        item['bits'] = sub_bits_valid + ('0' * (n_total - n_valid))
        item['dna'], item['seed'] = _encode_bits_to_dna_with_seed(item['bits'], dna_library_context=unchanged_context)
        cursor += n_valid
    if cursor != len(updated_plte_bits):
        raise ValueError("Updated PLTE bits were not fully consumed.")

    # Recompute CRC for PLTE: CRC(type_bytes + data_bytes)
    plte_type_bytes = b'PLTE'
    plte_data_bytes = bitstring_to_bytes(updated_plte_bits)
    crc_value = zlib.crc32(plte_type_bytes + plte_data_bytes) & 0xFFFFFFFF
    crc_bits = f'{crc_value:032b}'

    checksum_total_len = sum(c.get('valid_bits', len(c['bits'])) for c in checksum_items)
    if checksum_total_len != 32:
        raise ValueError(f"Unexpected PLTE checksum bit length: {checksum_total_len} (expected 32).")

    cursor = 0
    for item in checksum_items:
        n_total = len(item['bits'])
        n_valid = item.get('valid_bits', n_total)
        sub_bits_valid = crc_bits[cursor:cursor + n_valid]
        item['bits'] = sub_bits_valid + ('0' * (n_total - n_valid))
        item['dna'], item['seed'] = _encode_bits_to_dna_with_seed(item['bits'], dna_library_context=unchanged_context)
        cursor += n_valid

    return target


def update_plte_entry_in_dna_chunks(dna_chunks, palette_index: int, new_rgb, in_place=False):
    """
    Update one PLTE palette entry by index, then fix CRC and DNA encoding.

    Parameters
    ----------
    dna_chunks : list[dict]
        Output of png_to_dna_chunks().
    palette_index : int
        Palette slot to update (0-based).
    new_rgb : tuple[int, int, int] | list[int]
        Replacement RGB values in [0, 255].
    in_place : bool
        If True, mutate input list. Otherwise return a modified copy.
    """
    if not isinstance(palette_index, int) or palette_index < 0:
        raise ValueError("palette_index must be a non-negative integer.")
    if not isinstance(new_rgb, (tuple, list)) or len(new_rgb) != 3:
        raise ValueError("new_rgb must be a tuple/list of 3 integers: (R, G, B).")
    if any((not isinstance(v, int)) or v < 0 or v > 255 for v in new_rgb):
        raise ValueError("Each RGB component must be an integer in [0, 255].")

    # Determine palette size first for bounds checking.
    plte_indices = sorted({c['chunk_index'] for c in dna_chunks if c.get('png_chunk_type') == 'PLTE'})
    if not plte_indices:
        raise ValueError("No PLTE chunk found in dna_chunks.")
    if len(plte_indices) > 1:
        raise ValueError("Multiple PLTE chunks found; expected a single PLTE chunk.")
    plte_idx = plte_indices[0]

    plte_palette_items = [c for c in dna_chunks if c['chunk_index'] == plte_idx and c['field'] == 'palette']
    if not plte_palette_items:
        raise ValueError("PLTE chunk is missing palette field.")
    plte_palette_items.sort(key=lambda x: x['subchunk_index'])
    plte_bits = ''.join(c['bits'][:c.get('valid_bits', len(c['bits']))] for c in plte_palette_items)
    if len(plte_bits) % 24 != 0:
        raise ValueError(f"PLTE payload bit length is invalid: {len(plte_bits)}.")

    palette_len = len(plte_bits) // 24
    if palette_index >= palette_len:
        raise IndexError(f"palette_index {palette_index} out of range for palette size {palette_len}.")

    r, g, b = new_rgb
    replacement_bits = f"{r:08b}{g:08b}{b:08b}"
    start = palette_index * 24
    end = start + 24

    def _modify(bits):
        return bits[:start] + replacement_bits + bits[end:]

    return update_plte_bits_in_dna_chunks(
        dna_chunks,
        bit_modifier=_modify,
        in_place=in_place
    )
