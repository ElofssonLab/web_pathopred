import os
import sys
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from keras.models import load_model
from goatools import obo_parser

import calculation


if __name__ == "__main__":


    num_args = len(sys.argv)

    #usage is predictor.py fasta align1 align2 variants
    if num_args < 4:
        print('Invalid arguments\nUsage: %s <fasta> <align_low> <align_high> <variant file> <outdir>' % sys.argv[0])
        sys.exit(1)


    fasta_file = sys.argv[1]
    align_low = sys.argv[2]
    align_high = sys.argv[3]
    variant_file = sys.argv[4]
    outdir = sys.argv[5]

    #Get identifier from fasta query
    with open(fasta_file, 'r') as f:
        identifier = f.readline()
        identifier = identifier.strip()
        #query_XXXX
        identifier = identifier[7:]
    print(identifier)
    #Read variant file and get a prediction for each variant for this identifier
    all_variants = []
    current_identifier = ''
    
    with open(variant_file, 'r') as file:

        for line in file:

	    line = line.strip()
	    if line.startswith('>'):
    	    	current_identifier = line[1:]
    	    else:
            	variant = line
        	if current_identifier == identifier:
            	    all_variants.append(line)
    print('All variants')
    print(all_variants)

    #Get the entropies along the whole sequence
    entropy_vector = calculation.get_all_entropies_normalized(align_low)
    fig, ax = plt.subplots(figsize=(20,10), dpi=80)
    variant_positions = [int(var_line[1:-1]) for var_line in all_variants]
    print(variant_positions)
    line, = ax.plot(entropy_vector, color='cadetblue', marker='o', markevery=variant_positions, markersize=12,
    		    mec='black', mfc='None')
    		    
    for variant in all_variants:
        var_position = int(variant[1:-1])
       	ax.annotate(variant, xy=(var_position, entropy_vector[var_position]), xycoords='data',
       		    xytext=(0,-100), textcoords='offset points', arrowprops=dict(facecolor='black', shrink=0.05),
       		    horizontalalignment='right', verticalalignment='bottom')
    ax.set_xlabel('Residue')
    ax.set_ylabel('Entropy')

    plt.savefig(outdir + '/' + identifier + '_entropy.png')

    #Generate a csv with entropy and residue data
    data_array = calculation.create_entropy_matrix(align_low)
    np.savetxt(outdir + '/entropy_data.csv', data_array, delimiter=',', fmt="%s,%s,%s",
    		  header="pos,entropy,residue", comments='')
    #Now run a prediction on each variant

    #Load the GO obo file
    script_folder = os.path.dirname(os.path.realpath(__file__))
    go_obo = script_folder + "/go-basic.obo"
    go = obo_parser.GODag(go_obo)

    #Load frequency pickles
    go_term_freq_neutral = pickle.load( open( script_folder + '/go_term_freq_neutral_p2_refseq.p', 'rb') )
    go_term_freq_pathogenic = pickle.load( open( script_folder + '/go_term_freq_pathogenic_p2_refseq.p', 'rb') )

    #Load the sspred model
    sspred = load_model(script_folder + '/ss_pred_resnet_elu_nolr_dropout01_l26_large_v3_saved_model.h5')

    #Load prediction model
    prediction_model = load_model(script_folder + '/p2_la4_dr7.hdf5')
    #Load severity prediction model
    sev_prediction_model = load_model(script_folder + '/ps_transfer.h5')
    print('Running predictions for variants')

    #We need to make sure we get identifier in uniprot, refseq, ensembl, what we have / can get mappings for.
    #Currently they will enter as "query" from server -- need to add some solution for this
    for variant in all_variants:
	print('On variant: ' + variant)
        indata = calculation.generate_indata(identifier, variant, sspred, align_low, align_high, go, go_term_freq_pathogenic,
            go_term_freq_neutral)

        prediction = prediction_model.predict(indata, batch_size=1)
        print('Prediction:')
        print(prediction)
	severity_str = "N/A"
	if prediction[0,0] > 0.5:
    	    prediction_str = 'Pathogenic'
    	    #Also run the severity model in this case
    	    severity_prediction = sev_prediction_model.predict(indata, batch_size=1)
    	    if severity_prediction < 0.5:
        	    sev_str = 'Non-severe'
       	    else:
           	    sev_str = 'Severe'
    	else:
       	    prediction_str = 'Neutral'
        #Write the prediction to file
        with open(outdir + '/output_predictions', 'a') as out_file:
            out_file.write("%s\t%s\t%s\t%f\t%s\t%f\n" % (identifier, variant, prediction_str, prediction[0,0], sev_str,
            	severity_prediction[0,0]))


