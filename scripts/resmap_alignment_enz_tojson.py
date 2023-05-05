import gzip 
import csv
import json
import requests
import sys, glob, hashlib, os, os.path, traceback
import findpygly
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan, GlyTouCanNoPrefetch
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from pygly.GlycanBuilderSVGParser import GlycanBuilderSVG

wp = WURCS20Format()
gp = GlycoCTFormat()
sp = GlycanBuilderSVG()

glyeq = GlycanImageEqual()

gtc = GlyTouCan()

enzyme_alignment_file = sys.argv[1]
f=gzip.open(enzyme_alignment_file,'rb')

#f=gzip.open('/home/ddragovic/PyGly/smw/glycomotif/data/enzymealignments.nlinked.tsv.gz','rb')

json_dict = {} 

line_count = 0  
for row in csv.DictReader(f, delimiter='\t'):
    line_count = line_count +1
    print(line_count)
    # if line_count >= 30:
    #     break
    if row["Structure"] in json_dict:
        #print("hey")
        if row["StructureResidue"] != "-":
            if row["StructureResidue"] not in json_dict[row["Structure"]]["residue_map"]:
                json_dict[row["Structure"]]["residue_map"].append(row["StructureResidue"])
                if row["HumanEnzyme"] != "-":
                    human_enzymes = row["HumanEnzyme"].split(",")
                    for enzyme in human_enzymes:
                        enzyme = "Enzyme_" + enzyme
                        if enzyme in json_dict[row["Structure"]]["annotation"]:
                            #print(row["Motif"])
                            #print(json_dict[row["Motif"]]["annotations"])
                            json_dict[row["Structure"]]["annotation"][enzyme].append(row["StructureResidue"])
                            #print(json_dict[row["Motif"]]["annotations"])   
                        else: 
                           json_dict[row["Structure"]]["annotation"][enzyme] = [row["StructureResidue"]]
                if row["MouseEnzyme"] != "-":
                    mouse_enzymes = row["MouseEnzyme"].split(",")
                    for enzyme in mouse_enzymes:
                        enzyme = "Enzyme_" + enzyme
                        if enzyme in json_dict[row["Structure"]]["annotation"]:
                            #print(row["Motif"])
                            #print(json_dict[row["Motif"]]["annotations"])
                            json_dict[row["Structure"]]["annotation"][enzyme].append(row["StructureResidue"])
                            #print(json_dict[row["Motif"]]["annotations"])   
                        else: 
                           json_dict[row["Structure"]]["annotation"][enzyme] = [row["StructureResidue"]]
                            
                    #print(json_dict[row["Motif"]]["annotations"])    
                    #print(".")     
    else:
        if row["StructureResidue"] != "-" :
            #print(row["Motif"])
            if row["HumanEnzyme"]!= "-" or row["MouseEnzyme"]!= "-" :
                json_dict[row["Structure"]] = {}
                row_dict = {}
                row_dict['accession'] = row["Structure"]
                row_dict["residue_map"] = [row["StructureResidue"]]
                annotations_dict = {}
                if row["HumanEnzyme"] != "-":
                    human_enzymes = row["HumanEnzyme"].split(",")
                    for enzyme in human_enzymes:
                        enzyme = "Enzyme_" + enzyme
                        
                        #print(human_enzymes)  
                        #print(annotations_dict)
                        if enzyme not in annotations_dict:   
                            #print(enzyme)                             
                            annotations_dict[enzyme] = [row["StructureResidue"]]
                            #print(annotations_dict)
                            #print(annotations_dict)
                        else:
                            annotations_dict[enzyme].append(row["StructureResidue"])
                            #print(annotations_dict)
                            #print(annotations_dict)
                if row["MouseEnzyme"] != "-":
                    
                    mouse_enzymes = row["MouseEnzyme"].split(",")
                    #print(mouse_enzymes)
                    for enzyme in mouse_enzymes:
                        enzyme = "Enzyme_" + enzyme
                        #print(mouse_enzymes)  
                        #print(annotations_dict)
                        if enzyme not in annotations_dict:   
                            #print(enzyme)                             
                            annotations_dict[enzyme] = [row["StructureResidue"]]
                            #print(annotations_dict)
                            #print(annotations_dict)
                        else:
                            annotations_dict[enzyme].append(row["StructureResidue"])
                            #print(annotations_dict)
                            #print(annotations_dict)
                            
                #print(annotations_dict)            
                row_dict["annotation"] = annotations_dict
                json_dict[row["Structure"]] = row_dict
                #print(row_dict["annotations"])
                
        # line_count += 1
                #row_dict["annotations"] = {row["HumanEnzyme"]:row["MotifResidue"],row["MouseEnzyme"]:row["MotifResidue"]}
                #print(row_dict)
                #json_dict[row["Motif"]] = row_dict
   
   
first_value = next(iter(json_dict.values()))

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input
    
    
#directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/try'

#directory_path = '/home/ddragovic/git_pull_032023/PyGly/scripts/resmapalignjson04051'

#directory_path = '/home/ddragovic/PyGly/scripts/alignedjson'

#directory_path = '/home/ddragovic/PyGly/scripts/mappedalignedjson'

#directory_path = '/home/ddragovic/PyGly/scripts/prac2'

directory_path = sys.argv[2]

file_paths = glob.glob(directory_path + "/*.json")

for path in file_paths:
    with open(path, "r") as file:
        data = json.load(file)
        acc = os.path.splitext(os.path.basename(path))[0]
        try:
            #print(json_dict[acc]["annotation"])
            print(acc)
            sys.exit()
            enzyme_annotation = json_dict[acc]["annotation"]
            print(enzyme_annotation)
            data["annotations"]["EnzymeAlignments"] = enzyme_annotation
            #data_str = json.dumps(data, ensure_ascii=False)
            
        except KeyError:
            continue
            #data["annotations"]["enzymes"]="No Enzyme information found"
            #data_str = json.dumps(data, ensure_ascii=False)
            
        # if not os.path.exists('enz_data'):
        #     os.mkdir('enz_data')
        
        
        print(data["accession"])
        
       
        
        if not os.path.exists(sys.argv[3]):
            os.mkdir(sys.argv[3])
            
        filename = data["accession"] + ".json"
        filepath = os.path.join(sys.argv[3], filename)



        with open(filepath, "w") as json_file:
            json.dump(data, json_file,json_file,indent=4)
            
            


