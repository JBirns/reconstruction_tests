from sage.all import *
from cospectrality import *
from copy import copy
from key_cospectrum import *
import time 

def binary_strings_to_sage_matrix(binary_strings):
    matrix_rows = []
    
    # Process each binary string
    for line in binary_strings.split("\n"):
        if line.strip():  # Skip empty lines
            row = [int(bit) for bit in line.strip()]
            matrix_rows.append(row)
    
    # Return actual Sage matrix
    return Matrix(matrix_rows)

def extract_graph_blocks(input_string):
    # Split the input string into lines
    lines = input_string.strip().split('\n')
    
    # Remove the first and last lines (which contain sentences)
    lines = lines[1:-1]
    
    # Initialize a list to hold the blocks of 0's and 1's
    blocks = []
    current_block = []
    
    # Loop through each line and group them into blocks separated by empty lines
    for line in lines:
        if line.strip():  # If the line is not empty
            current_block.append(line.strip())
        else:
            # If there's an empty line, it signifies the end of the current block
            if current_block:
                blocks.append('\n'.join(current_block))
                current_block = []
    
    # If the last block is not added yet (no empty line at the end)
    if current_block:
        blocks.append('\n'.join(current_block))
    
    return blocks


def tree_structure(G, k):
    tree = key_cospectrum(G, k)[1]
    pen_bottom = [[] for i in range(k+1)]

    # layer = k
    for layer in range(k, 0, -1):
        for tup in Permutations([n for n in range(G.nrows())],Integer(layer-1)):
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            tstr = "tree"
            for i in tup2:
                tstr = tstr + "[1]" + str([i]) + "[0]"
            
            # tstr = tstr + "[1]"
            foot = eval(tstr)
            if foot not in pen_bottom[layer]:
                pen_bottom[layer].append(foot)
        
            colour = []
            uniques = []
            for x in foot[1]:
                if x not in uniques:
                    uniques.append(x)
            for val in uniques:
                colour.append(foot[1].count(val))
            colour.sort()

            if layer < k:
                exec(tstr + "[0] = " + "colour")
            if layer == k:            
                exec(tstr + " = " + "colour")


    #order tree branches canonically
    for layer in range(k-1, 0, -1):
        for tup in Permutations([n for n in range(G.nrows())],Integer(layer-1)):#####
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            tstr = "tree"
            for i in tup2:
                tstr = tstr + "[1]" + str([i]) + "[0]"
            tstr = tstr + "[1]"
            exec(tstr + ".sort()")
    return(tree)


def tree_compressor(G, k, compression_level, colour_type = 'position'):
    tree = key_cospectrum(G, k)[1]
    pen_bottom = [[] for i in range(k+1)]

    # layer = k
    for layer in range(k, k-compression_level, -1):
        for tup in Permutations([n for n in range(G.nrows())],Integer(layer-1)):#####
            tup2 = copy(tup)
            tup2 = re_index(tup2)
            tstr = "tree"
            for i in tup2:
                tstr = tstr + "[1]" + str([i]) + "[0]"
            
            # tstr = tstr + "[1]"
            foot = eval(tstr)
            if foot not in pen_bottom[layer]:
                pen_bottom[layer].append(foot)
            
            if colour_type == 'position':
                # #colour by position
                i = 0
                while pen_bottom[layer][i] != foot:
                    i = i + 1
                colour = -(i + 1) # to distinguish from existing colours
            if colour_type == 'shape':#colour by shape
                colour = []
                uniques = []
                for x in foot[1]:
                    if x not in uniques:
                        uniques.append(x)
                for val in uniques:
                    colour.append(foot[1].count(val))
                colour.sort()
            else:
                raise TypeError("colour only by shape or position")
            
            exec(tstr + "= " + "colour")

    t = 0
    while pen_bottom[0] == []:
        del pen_bottom[0]
        t = t+1

    layers = []
    for i in range(len(pen_bottom)):
        layers.append("Layer " + str(i + t) + ":")
        layers.append(pen_bottom[i])
    return(tree, "<-- tree | layers compressed -->" ,layers)