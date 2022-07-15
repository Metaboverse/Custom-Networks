""" Generate a custom network for Metaboverse for M. leprae
"""

import os 
import sys
import re
import json
import numpy as np
import pandas as pd

__path__ = os.getcwd()

reactions_url   = os.path.join(__path__, "source", "pntd.0007871.s003.reactions.txt")
pathways_url    = os.path.join(__path__, "source",  "pntd.0007871.s003.pathways.txt")
metabolites_url = os.path.join(__path__, "source",  "pntd.0007871.s003.metabolites.txt")
output_file = os.path.join("..", "..", "curated_organisms", "m_leprae", "m_leprae.json")

reactions   = pd.read_csv(reactions_url, sep="\t", header=None, encoding='windows-1252')
pathways    = pd.read_csv(pathways_url, sep="\t", header=None, encoding='windows-1252')
metabolites = pd.read_csv(metabolites_url, sep="\t", header=None, encoding='windows-1252')

model = {}
model["species"] = {}
model["reactions"] = {}
model["synonyms"] = {}
model["pathways"] = {}
model["components_database"] = {}

metabolite_dictionary = {}
for index, row in metabolites.iterrows():
    metabolite_dictionary[row[0]] = row[1]
    
model["synonyms"] = metabolite_dictionary

pathway_dictionary = {}
for index, row in pathways.iterrows():
    if row[0] != np.nan and row[1] != np.nan and str(row[1]) != "nan":
        pathway_dictionary[row[0]] = str(row[1]).strip("\"").split(", ")
    
model["pathways"] = pathway_dictionary

for index, row in reactions.iterrows():
    reaction_name = row[0]
    reaction = row[1]
    modifier = row[2].replace("(", "").replace(")", "")
    description = row[3]
    
    inputs = reaction.split("=")[0]
    if " + " in inputs:
        inputs = inputs.split(" + ")
    else: 
        inputs = [inputs]
    inputs = [i.replace(" ", "")  for i in inputs]
    
    outputs = reaction.split("=")[1]
    if " + " in outputs:
        outputs = outputs.split(" + ")
    else: 
        outputs = [outputs]
    outputs = [o.replace(" ", "")  for o in outputs]
    
    if "AND" in modifier.upper() or "OR" in modifier.upper():
        modifier = re.split(' AND | OR ', modifier)
    else:
        modifier = [modifier]
        
    for i in inputs:
        if i in metabolite_dictionary:
            name = metabolite_dictionary[i]
        else:
            name = i
        model["species"][i] = {"id": i, "is": i, "hasPart": [], "name": name, "type": "metabolite"}
        model["components_database"][i] = model["species"][i]
        
    for o in outputs:
        if o in metabolite_dictionary:
            name = metabolite_dictionary[o]
        else:
            name = o
        model["species"][o] = {"id": o, "is": o, "hasPart": [], "name": name, "type": "metabolite"}
        model["components_database"][o] = model["species"][o]
    
    modifier_list = []
    for m in modifier:
        m_stripped = m.replace(" ", "")
        m_name = m_stripped
        if "#Exchange:" in m_stripped:
            m_stripped = m_stripped.replace("#Exchange:", "")
        modifier_list.append([m_stripped, "catalyst"])
        model["species"][m_stripped] = {"id": m_stripped, "is": m_stripped, "hasPart": [], "name": m_name, "type": "modifier"}
        model["components_database"][m_stripped] = model["species"][m_stripped]
    
    model["reactions"][reaction_name] = {
        "name": reaction_name, 
        "reactants": inputs,
        "products": outputs,
        "modifiers": modifier_list,
        "description": description
    }

with open(output_file, 'w') as output:
    json.dump(model, output, indent=4)