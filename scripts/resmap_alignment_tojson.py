import csv
import gzip
import json
from itertools import islice
import sys
import glob
import io
import os

import time


start_time = time.time()

#alignment_file = '/home/ddragovic/git_pull_032023/PyGly/smw/glycomotif/data/motif_alignment.tsv.gz'

#alignment_file = '/home/ddragovic/PyGly/smw/glycomotif/data/motif_alignment.tsv.gz'

alignment_file = sys.argv[1]

def myopen(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename)
    return open(filename)



name_params = [0,1,2,3,4,5,6,7,8,9]

for name_param in name_params:  
    whole_alignfile_dict = {} 
    #print(name_param)
    with myopen(alignment_file) as f:
        next(f) 
        reader = csv.DictReader(f, delimiter='\t', fieldnames=['Motif', 'Structure' ,'Core_Inclusive' ,'Substructure_Inclusive' ,'Whole_Inclusive' ,'Non_Red_Inclusive' ,'Core_Strict' ,'Substructure_Strict' , 'Whole_Strict' ,'Non_Red_Strict'])
        ignore_key = ['Motif','Structure']  
        #print(name_param)
        for row in reader:
            #if row["Structure"][1] == str(name_param):
            if row["Structure"][1] == str(name_param):
            
                #print(row["Structure"][1])
                #print(row["Structure"])
                if row["Structure"] not in whole_alignfile_dict:
                    whole_alignfile_dict[row['Structure']] = {'motif_alignments': {},}
                    
                    
                    for key in row:
                        if key not in ignore_key:
                            if row[key] != 'N':
                                id_list = row[key].replace("Y:", "").replace(":", ",").split(",")
                                id_list = [elem for elem in id_list if elem != '']
                                #alignment_dict[key] = id_list
                                motif_alignment_info = row["Motif"]+"_" + key
                                whole_alignfile_dict[row['Structure']]['motif_alignments'][motif_alignment_info] = id_list
                                # print(whole_alignfile_dict)
                                
                                
                                # print("motif_alignment_info",motif_alignment_info)
                                # print("id_list",id_list)
                                # print("key",key, row["Motif"])
                    #print(alignment_dict)
                    
                     
                    
                else: 
                    motif_dict = whole_alignfile_dict[row['Structure']]['motif_alignments']
                    alignment_dict ={}
                    for key in row:
                        if key not in ignore_key:
                            if len(row[key]) > 1:
                                id_list = row[key].replace("Y:", "").replace(":", ",").split(",")
                                id_list = [elem for elem in id_list if elem != '']
                                #alignment_dict[key] = id_list
                                motif_alignment_info = row["Motif"]+"_" + key
                                whole_alignfile_dict[row['Structure']]['motif_alignments'][motif_alignment_info] = id_list
                    
                    # if row['Structure'] == 'G00026MO':
                    #     print(row['Structure'])
                    #     print(whole_alignfile_dict[row['Structure']]['motif_alignments'])
    
                           
                    
        #data = whole_alignfile_dict
        
            #directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/prac6json'

        #directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/jsonpracadir'

        #directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/resmapjsons0405'
        
        #directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/resmapjsons0405'
        
        #directory_path = '/home/ddragovic/PyGly/scripts/mappedjson'
        
        #directory_path = '/home/ddragovic/PyGly/scripts/prac1'
        
        directory_path = sys.argv[2]
        
        file_paths = glob.glob(directory_path + "/*.json")
        
       
       
     

        for path in file_paths:
            with open(path, "r") as file:
                resmap_json = json.load(file)
                acc = os.path.splitext(os.path.basename(path))[0]
                #print(acc)
                #print(resmap_json)
                #resmap_json['annotations'] = {}
                resmap_json['annotations']['motif_alignments'] ={}
                #print(resmap_json)
                for entry in whole_alignfile_dict:
                    #print(entry)
                    if entry == acc:
                        
                        #print(whole_alignfile_dict[entry])
                        resmap_json['annotations']['motif_alignments']= whole_alignfile_dict[entry]['motif_alignments']
                        #print(whole_alignfile_dict[entry])
                        print(entry)
                        cwd = os.getcwd()
                        print(cwd)
                        parent_dir = "/home/ddragovic/PyGly/scripts" ## change no need for this..
                        new_dir = sys.argv[3]
                        new_dir_path = os.path.join(parent_dir, new_dir)
                        if not os.path.exists(new_dir_path):
                            print(path)
                            os.mkdir(new_dir_path)  
                        jsonfilename = acc + ".json"
                        jsonfilepath = os.path.join(new_dir_path, jsonfilename)
                        
                        with open(jsonfilepath, "w") as json_file:
                            json.dump(resmap_json, json_file,indent=4)
                
            
        break
        
end_time = time.time()

total_time = end_time - start_time
print("Total time taken: " + str(total_time) + " seconds")                
