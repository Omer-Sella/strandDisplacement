import sys
PATH_TO_PROJECT = "/rds/general/user/osella/home/strandDisplacement/"
from utils import *
from pathlib import Path
sys.path.insert(0, PATH_TO_PROJECT)
import argparse
import json

def encodeImage(png_path,save_to_path, K_BITS = 200):
    # 1. Provide a Path to png file with a PLTE field
    #png_path = './strandDisplacement/imperialBlue_with_plte_2_colours.png'
    #png_path = 'imperialBlue_with_plte_2_colours.png'
    #png_path = 'imperialBlue_with_plte.png'

    # 2. Set the limit on the number of bits per chunk. This will affect how many bases are in every DNA strand.
    K_BITS = 200  # Placeholder constant; you can change this later. The number of bases per strand == K_BITS / 2
    
    # 3. Encode the binary chunks into DNA chunks 
    dna_chunks = png_to_dna_chunks(png_path, K_BITS)
    
    # 4. Save to a file
    with open(save_to_path, "w") as f:
        json.dump(dna_chunks, f, indent=4)
    
    return dna_chunks


#parser = argparse.ArgumentParser(description = "Encoding and decoding png images to and from DNA strands for Ilaria's strand displacement project.")

dna_chunks = encodeImage(png_path = '/rds/general/user/osella/home/strandDisplacement/imperialBlue_with_plte_2_colours.png', save_to_path = '/rds/general/user/osella/home/strandDisplacement/chunks_imperialBlue_with_plte_2_colours.json')
