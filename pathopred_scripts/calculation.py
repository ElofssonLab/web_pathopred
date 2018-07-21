import os
import os.path
import numpy as np
import math
import pickle
import urllib2
import sys

from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import AlignInfo
from Bio.SubsMat import MatrixInfo
from goatools import obo_parser
from urllib2 import URLError, HTTPError
from keras.models import Model


import _load_data

def get_go_terms(identifier, go):

    #to get the go terms for this identifier, access uniprot
    #we want e.g. http://www.uniprot.org/uniprot/?query=accession:Q13635&format=tab&columns=id,go-id

    identifier = identifier.strip()
    GOurl = "http://www.uniprot.org/uniprot/?query=accession:" + identifier + "&format=tab&columns=go-id"

    try:
        response = urllib2.urlopen(GOurl)
    except HTTPError as e:
        print('Uniprot: The server could not fulfill the request.')
        print('Uniprot: Error code: ', e.code)
        GOresult = None
        sys.exit(1)
    except URLError as e:
        print('Uniprot: We failed to reach a server.')
        print('Uniprot: Reason: ', e.reason)
        GOresult = None
        sys.exit(1)
    else:
        GOresult = response.read()
        response.close()

    #print(GOresult)

    #we want all go ids and ancestor ids
    all_go_ids = []

    #if not GOresult:
    #    return None

    nLineSplit = GOresult.split('\n')
    #then  we need grab second line
    nLineSplit = nLineSplit[1]
    #then we split on ; and strip the ids..
    go_ids = nLineSplit.split(';')
    go_ids = [goid.strip() for goid in go_ids]
    #print(go_ids)

    for goid in go_ids:
        #some ids are obsolete; check if they exist in the current go term database
        if goid in go:
            go_term = go[goid]
            all_go_ids.append(go_term.id)
            parents = go_term.get_all_parents()
            for term in parents:
                all_go_ids.append(go[term].id)

    return all_go_ids


def get_lr_sum(seq_name, go, go_term_freq_pathogenic, go_term_freq_neutral, refseq_mapping):

    #now we need to get the GO terms for this identifier. first map to uniprot if this mapping exists
    #if we are using a mapping in the first place
    if refseq_mapping is not None:
        if seq_name.startswith('NP') or seq_name.startswith('ENSG'):
            if seq_name in refseq_mapping:
                identifier = refseq_mapping[seq_name]
            else:
                return 'error'
        else:
            identifier = seq_name
    else:
        identifier = seq_name

    #then get the terms
    all_go_ids = get_go_terms(identifier, go)
    #if all_go_ids is None or len(all_go_ids) == 0:
    #    return 'error'

    all_go_ids = np.array(all_go_ids)
    all_go_ids = np.unique(all_go_ids)

    lr_sum = 0
    for term in all_go_ids:
        freq_path = 1
        if term in go_term_freq_pathogenic: 
            freq_path += go_term_freq_pathogenic[term]
        freq_neut = 1
        if term in go_term_freq_neutral:
            freq_neut += go_term_freq_neutral[term]

        lr_sum += math.log(freq_path / freq_neut, 2)

    return lr_sum



def get_ldata_matrix(a3m_name, mut_pos):

    #then we also want the features themselves, for both ranges
    self_info, part_entr, seq_feat = _load_data.process_a3m(a3m_name)

    #prep for adding below
    self_info = self_info[mut_pos-11:mut_pos+10, :]
    part_entr = part_entr[mut_pos-11:mut_pos+10, :]
    seq_feat = seq_feat[mut_pos-11:mut_pos+10, :]

    #if we have a size less than the window, pad out the array (it'll be the same for each)
    if self_info.shape[0] < 21:
        to_pad = 21 - self_info.shape[0]
        self_info = np.pad(self_info, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        part_entr = np.pad(part_entr, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        seq_feat = np.pad(seq_feat, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))

    return self_info, part_entr, seq_feat


def get_ss_matrix(a3m_name, mut_pos, sspred):

    self_info, part_entr, seq_feat = _load_data.process_a3m(a3m_name)

    #prepare and generate prediction
    self_info = np.expand_dims(self_info, axis=0)
    part_entr = np.expand_dims(part_entr, axis=0)
    seq_feat = np.expand_dims(seq_feat, axis=0)

    indata = [seq_feat, self_info, part_entr]
    predictions = sspred.predict(indata, batch_size=1)

    ss3 = np.array(predictions[0])
    ss6 = np.array(predictions[1])
    rsa = np.array(predictions[2])
    dihed = np.array(predictions[3])

    #then grab those around the position
    ss3 = ss3[0, mut_pos-11:mut_pos+10]
    ss6 = ss6[0, mut_pos-11:mut_pos+10]
    rsa = rsa[0, mut_pos-11:mut_pos+10]
    dihed = dihed[0, mut_pos-11:mut_pos+10]
    #if we have a size less than the window, pad out the array (it'll be the same for each)
    if ss3.shape[0] < 21:
        to_pad = 21 - ss3.shape[0]
        ss3 = np.pad(ss3, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        ss6 = np.pad(ss6, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        rsa = np.pad(rsa, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        dihed = np.pad(dihed, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))

    ss_matrix = np.concatenate((ss3,ss6,rsa,dihed), axis=1)

    return ss_matrix


def get_ss_matrix_sequence(a3m_name, mut_pos, sspred):

    _, _, seq_feat = _load_data.process_a3m(a3m_name)

    #prepare and generate prediction
    seq_feat = np.expand_dims(seq_feat, axis=0)

    #indata = [seq_feat]
    predictions = sspred.predict(seq_feat, batch_size=1)

    ss3 = np.array(predictions[0])
    ss6 = np.array(predictions[1])
    rsa = np.array(predictions[2])
    dihed = np.array(predictions[3])

    #then grab those around the position
    ss3 = ss3[0, mut_pos-11:mut_pos+10]
    ss6 = ss6[0, mut_pos-11:mut_pos+10]
    rsa = rsa[0, mut_pos-11:mut_pos+10]
    dihed = dihed[0, mut_pos-11:mut_pos+10]
    #if we have a size less than the window, pad out the array (it'll be the same for each)
    if ss3.shape[0] < 21:
        to_pad = 21 - ss3.shape[0]
        ss3 = np.pad(ss3, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        ss6 = np.pad(ss6, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        rsa = np.pad(rsa, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))
        dihed = np.pad(dihed, ((0, to_pad), (0,0)), 'constant', constant_values=(0,0))

    ss_matrix = np.concatenate((ss3,ss6,rsa,dihed), axis=1)

    return ss_matrix


def get_freq_matrix(a3m_name, mut_pos, mapping):

    alphabet = Gapped(IUPAC.extended_protein)
    alignment = AlignIO.read(open(a3m_name), "fasta", alphabet=alphabet)

    #collect frequencies as well as entr around position
    freq_matrix = np.zeros((21,24))
    i = 0
    for current_pos in range(mut_pos-10, mut_pos+11):

        freq_vector = np.zeros((1,24))
        all_symbols = []

        for record in alignment:
            if current_pos < len(record.seq):
                #print('Trying to access ' + str(current_pos))
                #print('Record length is ' + str(len(record.seq)))
                aa = record.seq[current_pos-1]
                all_symbols.append(aa)
            else:
                all_symbols.append('x')

        #make the list into a string
        sym_string = ''.join(all_symbols)
        all_counts = []
        for letter in mapping.keys():
            letter_count = sym_string.count(letter)
            all_counts.append(int(letter_count))
            #keep all symbols as frequencies
            letter_frequency = float(letter_count) / float(len(all_symbols))
            freq_vector[0,mapping[letter]] = letter_frequency

        #also get the entropy
        pos_entropy = 0
        for count in all_counts:
            #skip zeros
            if count is not 0:
                prop = float(count)/float(len(all_symbols))
                part_ent = prop * math.log(prop, 2)
                pos_entropy -= part_ent
        freq_vector[0,-1] = pos_entropy

        freq_matrix[i,:] = freq_vector

        i += 1

    return freq_matrix


def create_entropy_matrix(a3m_name):

    entr_vector = get_all_entropies_normalized(a3m_name)
    #then read the a3m and get the sequence..
    with open(a3m_name, 'r') as f:
        f.readline()
        sequence = f.readline()
        sequence = sequence.strip()

    positions = range(0, len(entr_vector))
    entr_vector = np.array(entr_vector)
    sequence = np.array(list(sequence))
    positions = np.array(positions)
    data_array = np.array([positions, entr_vector, sequence])
    data_array = np.transpose(data_array)
    #np.savetxt('data.csv', data_array, delimiter=',', fmt="%s,%s,%s", header="pos,entropy,residue", comments='')

    return data_array


def get_all_entropies_normalized(a3m_name):

    _, part_entr, _ = _load_data.process_a3m(a3m_name)
    entr_vector = -np.sum(part_entr, axis=1)
    return entr_vector


def get_all_entropies_non_normalized(a3m_name):


    mapping = {'-': 0, 'A': 1, 'B': 3, 'C': 2, 'D': 5, 'E': 4, 'F': 7,
                         'G': 6,'H': 9, 'I': 8, 'K': 10, 'L': 12, 'M': 11,
                         'N': 13, 'P': 15,'Q': 14, 'R': 17, 'S': 16, 'T': 18,
                         'V': 20, 'W': 19, 'Y': 21,
                         'O': 10, 'U': 2, 'Z': 4, 'X': 22}

    alphabet = Gapped(IUPAC.extended_protein)
    alignment = AlignIO.read(open(a3m_name), "fasta", alphabet=alphabet)

    entropy_vector = np.zeros((1,alignment.get_alignment_length()))
    all_symbols = []
    
    for current_pos in range(0, alignment.get_alignment_length()):
        for record in alignment:
	    aa = record.seq[current_pos]
	    all_symbols.append(aa)

        sym_string = ''.join(all_symbols)
        all_counts = []
        for letter in mapping.keys():
            letter_count = sym_string.count(letter)
            all_counts.append(int(letter_count))


        pos_entropy = 0
        for count in all_counts:
            if count is not 0:
                prop = float(count)/float(len(all_symbols))
                part_ent = prop * math.log(prop,2)
                pos_entropy -= part_ent

	entropy_vector[0, current_pos] = pos_entropy

    return entropy_vector


def generate_indata(seq_name, mutation, sspred, align_file_low, align_file_high, go, 
    go_term_freq_pathogenic, go_term_freq_neutral):


    mapping = {'-': 0, 'A': 1, 'B': 3, 'C': 2, 'D': 5, 'E': 4, 'F': 7,
                         'G': 6,'H': 9, 'I': 8, 'K': 10, 'L': 12, 'M': 11,
                         'N': 13, 'P': 15,'Q': 14, 'R': 17, 'S': 16, 'T': 18,
                         'V': 20, 'W': 19, 'Y': 21,
                         'O': 10, 'U': 2, 'Z': 4, 'X': 22}

    mut_string = mutation.strip()
    ref_aa = mut_string[0]
    alt_aa = mut_string[-1]
    mut_pos = int(mut_string[1:-1])

    ref_aa_idx = mapping[ref_aa]
    alt_aa_idx = mapping[alt_aa]

    ref_encode = np.zeros((1, 22))
    ref_encode[0, ref_aa_idx-1 ] = 1
    alt_encode = np.zeros((1,22))
    alt_encode[0, alt_aa_idx-1 ] = 1

    refseq_mapping = None
    lr_sum = get_lr_sum(seq_name, go, go_term_freq_pathogenic, go_term_freq_neutral, refseq_mapping)
    if lr_sum == 'error':
        return 'error'

    ref_alt_vector = np.concatenate((ref_encode, alt_encode), axis=1)

    #add the LR sum to the vector too
    ref_alt_vector = np.append(ref_alt_vector, lr_sum)
    ref_alt_vector = np.expand_dims(ref_alt_vector, axis=0)
    #then we also want the features themselves, for both ranges

    a3m_name = align_file_high
    #check if we have it or not
    if not os.path.isfile(a3m_name):
        print('Alignment file does not exist')
        sys.exit(1)
    self_info_1, part_entr_1, seq_feat_1 = get_ldata_matrix(a3m_name, mut_pos)
    a3m_name = align_file_low
    if not os.path.isfile(a3m_name):
        print('Alignment file does not exist')
        sys.exit(1)
    self_info_2, part_entr_2, seq_feat_2 = get_ldata_matrix(a3m_name, mut_pos)

    pentr_matrix = np.concatenate((part_entr_1, part_entr_2), axis=1)
    sinfo_matrix = np.concatenate((self_info_1, self_info_2), axis=1)
    
    a3m_name = align_file_low
    ss_matrix = get_ss_matrix(a3m_name, mut_pos, sspred)

    a3m_name = align_file_high
    freq_matrix_1 = get_freq_matrix(a3m_name, mut_pos, mapping)
    a3m_name = align_file_low
    freq_matrix_2 = get_freq_matrix(a3m_name, mut_pos, mapping)
    freq_matrix = np.concatenate((freq_matrix_1, freq_matrix_2), axis=1)

    ref_alt_vector = np.expand_dims(ref_alt_vector, axis=0)
    freq_matrix = np.expand_dims(freq_matrix, axis=0)
    sinfo_matrix = np.expand_dims(sinfo_matrix, axis=0)
    pentr_matrix = np.expand_dims(pentr_matrix, axis=0)
    seq_feat_2 = np.expand_dims(seq_feat_2, axis=0)
    ss_matrix = np.expand_dims(ss_matrix, axis=0)
    indata = [ref_alt_vector, freq_matrix, sinfo_matrix, pentr_matrix, seq_feat_2, ss_matrix]

    return indata

def calculate_sequence_version(seq_name, mutation, sample, i, sspred, alignFolder_a3m, h5file, go, 
    go_term_freq_pathogenic, go_term_freq_neutral, refseq_mapping):


    mapping = {'-': 0, 'A': 1, 'B': 3, 'C': 2, 'D': 5, 'E': 4, 'F': 7,
                         'G': 6,'H': 9, 'I': 8, 'K': 10, 'L': 12, 'M': 11,
                         'N': 13, 'P': 15,'Q': 14, 'R': 17, 'S': 16, 'T': 18,
                         'V': 20, 'W': 19, 'Y': 21,
                         'O': 10, 'U': 2, 'Z': 4, 'X': 22}

    mut_string = mutation.strip()
    ref_aa = mut_string[0]
    alt_aa = mut_string[-1]
    mut_pos = int(mut_string[1:-1])

    ref_aa_idx = mapping[ref_aa]
    alt_aa_idx = mapping[alt_aa]

    ref_encode = np.zeros((1, 22))
    ref_encode[0, ref_aa_idx-1 ] = 1
    alt_encode = np.zeros((1,22))
    alt_encode[0, alt_aa_idx-1 ] = 1

    lr_sum = get_lr_sum(seq_name, go, go_term_freq_pathogenic, go_term_freq_neutral, refseq_mapping)
    if lr_sum == 'error':
        return 'error'

    ref_alt_vector = np.concatenate((ref_encode, alt_encode), axis=1)

    #get a BLOSUM matrix
    alphabet = IUPAC.protein
    blosum_matrix = MatrixInfo.blosum62

    #add the LR sum to the vector too
    ref_alt_vector = np.append(ref_alt_vector, lr_sum)
    #let's also append the BLOSUM score for the subsitution
    pair = (ref_aa, alt_aa)
    if pair not in blosum_matrix:
        pair = (alt_aa, ref_aa)
    blosum_sub = blosum_matrix[pair]
    ref_alt_vector = np.append(ref_alt_vector, blosum_sub)
    ref_alt_vector = np.expand_dims(ref_alt_vector, axis=0)


    #this really should only be the sequence version; but do this for a trial run
    alignfolder_2 = alignFolder_a3m + '90-0/'
    a3m_name = alignfolder_2 + seq_name + "_90-0.a3m"
    ss_matrix = get_ss_matrix_sequence(a3m_name, mut_pos, sspred)
    _, _, seq_feat = _load_data.process_a3m(a3m_name)

    
    #we should get the BLOSUM scores for each AA in the window
    #we also need to read the sequence
    #manual specify
    prot_folder = "/scratch_ssd/Data/proteins_p2/"
    with open(prot_folder + seq_name + '.fasta','r') as file:
        file.readline()
        full_sequence = file.read()
        full_sequence = full_sequence.strip()
        full_sequence = ''.join(full_sequence.split())

    blosum_window = []
    for current_pos in range(mut_pos-10, mut_pos+11):
        if current_pos < len(full_sequence):
            position_aa = full_sequence[current_pos]
        else:
            position_aa = 'X'
        #sometimes this aa is not in standard 20 -- can we replace all these with X?
        if position_aa not in alphabet.letters:
            position_aa = 'X'
        blosum_vector = []
        for letter in alphabet.letters:
	    pair = (position_aa, letter)
	    if not pair in blosum_matrix:
	        pair = (letter, position_aa)
            blosum_vector.append( blosum_matrix[pair] )
        blosum_window.append(blosum_vector)
    
    blosum_window = np.array(blosum_window)

    h5file.create_array(sample, 'refalt', ref_alt_vector, "Ref alt LR")
    h5file.create_array(sample, 'blosum', blosum_window, "freqs and entr")
    h5file.create_array(sample, 'seq', seq_feat, "AA seq")
    h5file.create_array(sample, 'ssmatrix', ss_matrix, "predicted ss")

    return 'ok'
