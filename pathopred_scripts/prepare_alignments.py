from Bio.Blast import NCBIXML
from Bio.Align.Applications import MafftCommandline
import os
import os.path
import subprocess
from subprocess import check_output
import numpy as np
import glob
import sys

#split the results up into blocks of similarity ranges and e-val thresholds, store in fasta files
def mask_blast(out_folder, seq_name, xml_file, fasta_file, E_VALUE_THRESH, identity_range):


    blast_save_name = out_folder + seq_name + '_' + str(identity_range[0]) + '-' + str(identity_range[1]) + '._blast.fasta'

    print("Masking " + seq_name + " blast with E-val " + str(E_VALUE_THRESH) + " and idp range " + str(identity_range[0]) + ' - ' + str(identity_range[1]))
    xml_handle = open(xml_file)
    blast_record = NCBIXML.read(xml_handle)
    xml_handle.close()
    
    blast_save_file = open(blast_save_name, 'w')

    #remember to include the original sequence
    with open(fasta_file,'r') as file:
        file.readline()
        full_sequence = file.read()
        full_sequence = full_sequence.strip()
        full_sequence = ''.join(full_sequence.split())
        blast_save_file.write(">%s\n%s\n" % ('ORIGINAL_' + seq_name,full_sequence))

    tot_in_range = 0
    #get _all_ the identifiers, then batch search them.
    identifier_list = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= E_VALUE_THRESH:
                hsp_idp = float(hsp.identities) / float(hsp.align_length) * 100
                #print(hsp_idp)
                if hsp_idp <= identity_range[0] and hsp_idp > identity_range[1]:
                    tot_in_range += 1
                    #grab the GI identifier. split on | in hit id, grab second column for gi number
                    id_split = alignment.hit_id.split('|')
                    gi_id = id_split[1]
                    #print(gi_id)
                    identifier_list.append(gi_id)
                        

    if tot_in_range == 0:
        print("None in range for " + seqName)

    #then let's make it into a comma separated string
    id_string = ','.join(identifier_list)

    try:
        out = check_output(['blastdbcmd', '-entry', id_string, '-outfmt', '%f'])
        blast_save_file.write(out)
    except KeyboardInterrupt:
        print('Exiting..')
        sys.exit(0)
    except:
        print('Could not find identifier in database')

    blast_save_file.close()

    return blast_save_name

def align_sequence(out_folder, blast_file):

    #make sure downloaded file is valid
    with open(blast_file,'r') as file:
        file.readline()
        checkline = file.readline()
        if 'html' in checkline:
            print('invalid file')
            sys.exit(1)

    out_file = blast_file.replace('_blast.fasta','clustal')
    script_folder = os.path.dirname(os.path.realpath(__file__))

    mafft_exec = script_folder + "/mafft-linux64/mafft.bat"
    print(mafft_exec)
    #thread option is not available on older versions of biopython
    mafft_cline = MafftCommandline(cmd=mafft_exec, input=blast_file, auto=True, clustalout=True)
    assert os.path.isfile(mafft_exec), "MAFFT executable missing"
    print("Aligning with..")
    print(mafft_cline)
    try:
        stdout, stderr = mafft_cline()
        with open(out_file, 'w') as outmafft:
            outmafft.write(stdout)
    except:
        print(blast_base + " need anysymbol")
        with open(out_file, 'w') as outmafft:
            print('Running with anysymbol..')
            subprocess.call([mafft_exec, '--auto', '--thread', '-1', '--clustalout', '--anysymbol', blast_file], stdout = outmafft)

    return out_file



def convert_alignment(out_folder, align_file):

    script_folder = os.path.dirname(os.path.realpath(__file__))
    out_file = align_file.replace('clustal','a3m')
    pipe = subprocess.call(['perl', script_folder + '/reformat.pl', 'clu', 'a3m', align_file, out_file, '-r', '-l', '30000'])

    
if __name__ == "__main__":

    num_args = len(sys.argv)

    if num_args < 3:
        print('Invalid arguments\nUsage: %s <fasta> <blastxml> <outfolder>' % sys.argv[0])
        sys.exit(1)

    #fasta first arg, xml is given as second arg, outfolder as third arg
    fasta_file = sys.argv[1]
    xml_file = sys.argv[2]
    out_folder = sys.argv[3]

    identifier = os.path.basename(fasta_file)
    identifier = os.path.splitext(identifier)[0]
    #given a single xml, sort and prepare alignment

    #First run the 100-90 range, then the 90-0 range
    b_file_1 = mask_blast(out_folder, identifier, xml_file, fasta_file, 0.001, [100,90])
    b_file_2 = mask_blast(out_folder, identifier, xml_file, fasta_file, 0.001, [90,0])

    #then align everything
    align_file_1 = align_sequence(out_folder, b_file_1)
    align_file_2 = align_sequence(out_folder, b_file_2)

    #also convert to a3m
    convert_alignment(out_folder, align_file_1)
    convert_alignment(out_folder, align_file_2)




                


    
