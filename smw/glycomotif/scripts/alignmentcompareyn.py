import csv
import sys


def compare(first, second):
    
    sharedKeys = set(first.keys()).intersection(second.keys())
    for key in sharedKeys:
        if first[key] != second[key]:
            
            print('Key: {}, Value 1: {}, Value 2: {}'.format(key, first[key], second[key]))
            diff = [i for i, (x, y) in enumerate(zip(first[key], second[key])) if x != y]
            print(diff)
      
import gzip
def myopen(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename)
    return open(filename)

def match_tsv_to_dictionary(input_file,with_indices):
    alignments = {}
    with myopen(input_file) as csv_f:
        for row in csv.DictReader(csv_f, delimiter='\t'):
            Motif = row['Motif']
            Motif_Structure = row['Motif'] + '_' + row['Structure']
            if with_indices == False:
                alignments[Motif_Structure] = [row['Core_Inclusive'], row['Substructure_Inclusive'],row['Whole_Inclusive'], row['Non_Red_Inclusive'],row['Core_Strict'], row['Substructure_Strict'],row['Whole_Strict'], row['Non_Red_Strict']]
            if with_indices == True:
                alignments[Motif_Structure] = [row['Core_Inclusive'][0], row['Substructure_Inclusive'][0],row['Whole_Inclusive'][0], row['Non_Red_Inclusive'][0],row['Core_Strict'][0], row['Substructure_Strict'][0],row['Whole_Strict'][0], row['Non_Red_Strict'][0]]
        
    return(alignments)

    
def get_glycan_acc(input_file):
    glycan_acc_set = set()
    with myopen(input_file) as csv_f:
        for row in csv.DictReader(csv_f, delimiter='\t'):
            structure = row['Structure']
            if structure not in glycan_acc_set:
                glycan_acc_set.add(structure)
    return(glycan_acc_set)
            
new_file = sys.argv[1]
old_file = sys.argv[2]

# glycan_acc_set_old = get_glycan_acc(old_file)
# glycan_acc_set_new= get_glycan_acc(new_file)
# glyc_acc_diff = glycan_acc_set_new.difference(glycan_acc_set_old)

new_alignmnets = match_tsv_to_dictionary(new_file,with_indices=True)
old_alignments = match_tsv_to_dictionary(old_file,with_indices=False)


compare(new_alignmnets,old_alignments)



