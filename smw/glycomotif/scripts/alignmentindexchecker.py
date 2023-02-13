
import csv
import sys


def get_indices(row):
    if row == 'N':
        indices = set()
    else:
        indices = set(row[2:].split(','))   
    return(indices)  

    
def row_to_dictionary(row):

    
    alignments_indices={}

    alignments_indices.update( {'loose_core' : get_indices(row[0])} )

    alignments_indices.update( {'loose_substructure' : get_indices(row[1])} )

    alignments_indices.update( {'loose_whole' : get_indices(row[2])} )

    alignments_indices.update( {'loose_nred' : get_indices(row[3])} )

    alignments_indices.update( {'strict_core' : get_indices(row[4])} )

    alignments_indices.update( {'strict_substructure' : get_indices(row[5])} )

    alignments_indices.update( {'strict_whole' : get_indices(row[6])} )

    alignments_indices.update( {'strict_nred' : get_indices(row[7])} )

    return(alignments_indices)


def set_checker(input_list):
    match_indices =  row_to_dictionary(input_list)
  
    error = False
    if not (match_indices["loose_core"].issubset(match_indices["loose_substructure"])):
    #print("error:", row)
        print(1)
        error = True
    if not (match_indices["loose_nred"].issubset(match_indices["loose_substructure"])):
    #print("error:", row)
        print(2)
        error = True
        
    if not (match_indices["strict_core"].issubset(match_indices["loose_core"])):
    #print("error:", row)
       print(3)
       error = True

    if not (match_indices["strict_core"].issubset(match_indices["strict_substructure"])):
    #print("error:", row)
        print(4) 
        error = True  

    if not (match_indices["strict_nred"].issubset(match_indices["strict_substructure"])):
    #print("error:", row)
        print(5)
        error = True

    if not (match_indices["strict_nred"].issubset(match_indices["loose_nred"])):
    #print("error:", row)
        print(6) 
        error = True 
    if not (match_indices["strict_substructure"].issubset(match_indices["loose_substructure"])):
    #print("error:", row)
        print(7)
        error = True
        
    return(error)
       
 
if __name__ == "__main__":
    input_file = sys.argv[1]
    with open(input_file) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        next(tsv_file, None)
        for line in tsv_file:
            #print(line[2:])
            set_checker(line[2:])

