
import csv
import sys

def get_indices(value):
    indices = set()
    if value == 'Y':
        for ind in value.split(':')[1:]:
            indices.update(ind.split(','))
    return(indices)  

keys = filter(None,"""
loose_core
loose_substructure
loose_whole
loose_nred
strict_core
strict_substructure
strict_nred
""".split())
    
def row_to_dictionary(row):
    alignments_indices={}
    for i,k in enumerate(keys):
        alignments_indices[k] = get_indices(row[i])
    return alignments_indices

comparisons = filter(None,"""
loose_whole <= loose_core
loose_whole <= loose_nred
loose_whole <= loose_substructure
loose_core <= loose_substructure
loose_nred <= loose_substructure
strict_whole <= strict_core
strict_whole <= strict_nred
strict_whole <= strict_substructure
strict_core <= strict_substructure
strict_nred <= strict_substructure
strict_core <= loose_core
strict_nred <= loose_nred
strict_substructure <= loose_substructure
""".split())

def set_checker(input_list):
    match_indices =  row_to_dictionary(input_list)

    for comp in comparisons:
        splcomp = comp.split()
        if len(splcomp) != 3:
            continue
        if splcomp[1] == "<=":
            if not match_indices[splcomp[0]] <= match_indices[splcomp[2]]:
                print >>sys.stderr, "Error: %s:%s NOT <= %s:%s"%(splcomp[0],",".join(sorted(match_indices[splcomp[0]])),
                                                                 splcomp[2],",".join(sorted(match_indices[splcomp[2]])))
                return False
        else:
            return False
    return True
 
if __name__ == "__main__":
    input_file = sys.argv[1]
    with open(input_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        next(tsv_file, None)
        for line in tsv_file:
            #print(line[2:])
            set_checker(line[2:])

