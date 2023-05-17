#!/bin/env python2

import json
import glob 
import io
import os
import requests
import sys, glob, hashlib, os, os.path, traceback
import findpygly
from pygly.Glycan import *
from pygly.GlycanResource import *
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual
from pygly.GlycanFormatter import *
from pygly.GlycanBuilderSVGParser import GlycanBuilderSVG

ctable = ResidueCompositionTable()
iupacSym = IUPACSym()

ip = IUPACLinearFormat()
#directory_path = "/data2/projects/GlyGen/APIFramework/src/Application/Glymage/work/wurcs"

directory_path = sys.argv[1]

#svg_directory_path = "/data2/projects/GlyGen/APIFramework/src/Application/Glymage/work/snfg/extended"

wp = WURCS20Format()
sp = GlycanBuilderSVG()
glyeq = GlycanImageEqual()



# Use glob to get a list of file paths in the directory
file_paths = glob.glob(directory_path + "/*.txt")



# Loop through the file paths and read the contents of each file
for path in file_paths:
    with open(path, "r") as file:
        WURCS_Seq = file.read()
        acc = os.path.splitext(os.path.basename(path))[0]
        structure_dict = {}
        print(acc)
        structure_dict['accession'] = acc
        idmapsdict={}
        svg_directory_path = sys.argv[2]
        
# Set the file name
        svg_file_name = acc + ".svg"

# Navigate to the directory
        os.chdir(svg_directory_path)

# Open the file
        if not os.path.exists(svg_file_name):
            idmapsdict='Error'
            print(acc, "no path found")
            continue
        if os.path.exists(svg_file_name):  
            with open(svg_file_name, "r") as svg_file:
            
               
    # Read the contents of the file
                svg_file_contents = svg_file.read()
          
            seqstr = svg_file_contents
            
                
                
            seqhash = hashlib.md5(seqstr).hexdigest().lower()
            try:
                gly = sp.toGlycan(seqstr)
                print("222222222")
                fmt = 'SVG'
            except GlycanParseError:
                idmapsdict='Error'
                print("Error")
                traceback.print_exc(file=sys.stdout)
                continue
            except KeyError: #added
                idmapsdict='Error'
                print("Error")
                continue
            
            
            try: 
                cannongly = wp.toGlycan(WURCS_Seq)
                idmapsdict='Error'
            except WURCS20ParseError:
                idmapsdict='Error'
                continue
            
            for m in cannongly.all_nodes():
                print(m.composition(ctable))
           
            print(cannongly.iupac_composition())
            nodeiterable = [ cannongly.root() ]
            
            for m in cannongly.all_nodes():
                try:
                    sym = iupacSym.toStr(m)
                    print(sym)
                except KeyError:
                    sym = None
                    print(sym)
            
            
            for m in nodeiterable:
                try:
                    sym = iupacSym.toStr(m)
                    print(sym)
                except KeyError:
                    sym = None
                    print(sym)
            
            
            
            
            if not cannongly:
                idmapsdict='Error'
                print("Error")
                # idmapids_dict = "Error"
                # data['annotations']['residue_mapping'] = idmapids_dict
                # res_map_dict["residueids"] = idmapids_dict
                # json_dict[json_acc]['residue_map'] = res_map_dict
                continue
            idmap = []
            
            if glyeq.eq(gly,cannongly,idmap=idmap):
            
                bad = False
                glyids = [ m.external_descriptor_id() for m in gly.all_nodes() ]
                idmapglyids = [ t[0].external_descriptor_id() for t in idmap ]
                cannonglyids = [ m.id() for m in cannongly.all_nodes() ]
                cannongly_iupacids = [iupacSym.toStr(m) for m in cannongly.all_nodes() ]
                print(cannonglyids)
                print(cannongly_iupacids)
                iupac_dict = {}
                ####################################
                for i in range (0,len(cannonglyids)):
                    print(cannonglyids[i])
                    print(cannongly_iupacids[i])
                    can_id = str(cannonglyids[i])
                    iupac_id = cannongly_iupacids[i]  
                    
                    if iupac_id not in iupac_dict: 
                        iupac_dict[iupac_id] = [can_id]
                    else: 
                        iupac_dict[iupac_id].append(can_id)
                print(iupac_dict)
    
                #     if can_id not in iupac_dict: 
                #         iupac_dict[can_id] = iupac_id
                #     else: 
                #         iupac_dict[can_id].append(iupac_id)
                # print(iupac_dict)
            
                
                
                
                cannonidmapglyids = [ t[1].id() for t in idmap ]
                idmapids = [ (t[0].external_descriptor_id(),t[1].id()) for t in idmap ]
                if len(glyids) != len(idmap) or len(cannonglyids) != len(idmap):
    
                    print("Sequence %s:%s has an incorrect number of nodes in idmap: %d, %d, %d"%(acc,seqhash,len(glyids),len(idmap),len(cannonglyids)))
                    bad = True
                if len(glyids) != len(set(glyids)) or len(cannonglyids) != len(set(cannonglyids)):
                    print("Sequence %s:%s has non-unique monosaccharide ids: %s,%s"%(acc,seqhash,glyids,cannonglyids))
                    bad = True
                if len(idmapglyids) != len(set(idmapglyids)) or len(cannonidmapglyids) != len(set(cannonidmapglyids)):
                    print("Sequence %s:%s idmap has non-unique monosaccharide ids: %s"%(acc,seqhash,idmapids))
                    bad = True
                if len(idmapids) != len(set(idmapids)):
                    print("Sequence %s:%s idmap has non-unique monosaccharide id pairs: %s"%(acc,seqhash,idmapids))
                    bad = True
                for mi,mj in idmap:
                    if not glyeq.monoeq(mi,mj):
                        print("Sequence %s:%s:%s has a mismatched monosaccharide"%(acc,seqhash,fmt))
                        print(mi)
                        print(mj)
                        bad = True
                if not bad:
                    idmapsdict=idmapids
                    
            
                    idmapids_dict = {}
                    ####################################
                    for i in idmapids:
                        can_id = str(i[1])
                        print(type(can_id)) 
                        svg_id = i[0]  
                        if i[1] not in idmapids_dict: 
                            idmapids_dict[i[1]] = [i[0]]
                        else: 
                            idmapids_dict[i[1]].append(i[0])
                            
                  
                    
                                       
                    for l in cannongly.all_links():
                        for parent_svgid in idmapids_dict[l.parent().id()]:
                            for child_svgid in idmapids_dict[l.child().id()]:
                                svgidbase,parent_svgid1 = parent_svgid.rsplit(':',1)
                                svgidbase = svgidbase.split('-',1)[1]
                                child_svgid1 = child_svgid.rsplit(':',1)[1]
                                #print(parent_svgid1,child_svgid1)
                               #print(svgidbase)
                                link_id = str(l.parent().id()) + "-" + str(l.child().id())
                                #print(link_id)
                                svg_link_id = "l-1:" + str(parent_svgid1) + ","+ str(child_svgid1)
                                #print(svg_link_id)
                                if link_id not in idmapids_dict:
                                    idmapids_dict[link_id] = [svg_link_id]
                                else: 
                                    idmapids_dict[link_id].append(svg_link_id)
                    print(idmapids_dict)
                            
                    idmapsdict = idmapids_dict
                    
                    
                            
                            
                    
                        
                    
                    
                
                else:
                    idmapsdict='Error'
                    print("Error")
                

            else:
                idmapsdict='Error'
                print("Error")
               
        else: 
            idmapsdict='Error'
            print("Error")
            
            
    structure_dict['residue_map'] = {}
    structure_dict['annotations'] = {}   ## added       
    structure_dict['residue_map']['residueids'] = {}   
    structure_dict['residue_map']['seqhash'] = seqhash    
    structure_dict['residue_map']['residueids'] = idmapsdict 
    structure_dict['annotations']['IUPAC'] = iupac_dict 
    
    if structure_dict['residue_map']['residueids'] !="Error":
        print(structure_dict)
        cwd = os.getcwd()
        print(cwd)
        parent_dir = "/home/ddragovic/PyGly/scripts"
        new_dir = sys.argv[3]
        new_dir_path = os.path.join(parent_dir, new_dir)
        if not os.path.exists(new_dir_path):
            print(path)
            os.mkdir(new_dir_path)  
        jsonfilename = acc + ".json"
        jsonfilepath = os.path.join(new_dir_path, jsonfilename)
        with open(jsonfilepath, "w") as json_file:
            json.dump(structure_dict, json_file,indent=4)
                
