#!/usr/bin/env python3

import os
import sys
import csv
import pprint
import requests
import pandas as pd

###########################
#This drops all of the rows that have an NA in the modifications thus no mods to XML
datafile2 = pd.read_csv('AllProteinGroups.csv')
datafile3 = datafile2.loc[datafile2['Modification Info List'].dropna().index]
 
#This function will be applied row by row and clean the modification column    
def process_modifications(dataframe_accession, dataframe_column_split):
    df = pd.DataFrame(columns = ['accession', 'aa#', 'modification', 'aa'])
    for item in dataframe_column_split:
        if 'DECOY' in dataframe_accession:
            continue
        else:
            modification_info = []
            modification_info.append(dataframe_accession)
        if 'Fe' in item:
            data = item.split(',info')[0]
            aa = data.split('[')[0]
            modification = data.split('[',1)[1].split(' on ')[0]
            site = data.split('[',1)[1].split(' on ')[1]
            modification_info.append(aa)
            modification_info.append(modification)
            modification_info.append(site)
            df_length = len(df)
            df.loc[df_length] = modification_info
        else:
            data = item.split(',info')[0]
            aa = data.split('[')[0]
            mod = data.split('[')[1]
            modification = mod.split(' on ')[0]
            site = mod.split(' on ')[1]
            aa = aa.strip('#aa')
            modification_info.append(aa)
            modification_info.append(modification)
            modification_info.append(site)
            df_length = len(df)
            df.loc[df_length] = modification_info
    return df
   
   
#Applying the function to the dataframe
final_df = pd.concat(datafile3.apply(lambda x: process_modifications(x['Protein Accession'], x['Modification Info List'].split(';')),axis=1).values.tolist())


#Dictionary of AA and their full name
residue_dictionary = {'A':'Alanine','R':'Arginine','N':'Asparagine',
                      'D':'Aspartic acid','C':'Cysteine','E':'Glutamic Acid',
                      'Q':'Glutamine','G': 'Glycine','H':'Histidine','I': 'Isoleucine',
                      'L':'Leucine','K':'Lysine','M':'Methionine','F':'Phenylalanine',
                      'P':'Proline','S':'Serine','T':'Threonine','W':'Tryptophan',
                      'Y':'Tyrosine','V': 'Valine'}

#Add full aa name to dataframe
final_df['aa_name'] = final_df.apply(lambda x: residue_dictionary[x['aa']],axis=1)

def aa_sub_split(modification_entry):
    if '->' in modification_entry:
        l = modification_entry.split('->')[1]
        return l
    else:
        return modification_entry

final_df['modification'] = final_df.apply(lambda x: aa_sub_split(x['modification']),axis=1)

#Request from uniprot and write the XML files. Need to edit this because it does pull down 0 xmls
#Can also add the xml edits from below in the function
def uniprot_request(accession_list):
    for accession in accession_list:
        if 'DECOY' in accession:
            continue
        else:
            xml_file = accession + '.xml'
            f = open(xml_file, 'w+')
            url = 'https://www.uniprot.org/uniprot/%s.xml'  % accession
            result = requests.get(url)
            f.write(result.text)
            f.close()
   



if os.path.isdir("appended_xml_files") == False:
    os.mkdir("appended_xml_files")
    os.chdir("appended_xml_files")
else:
    os.chdir("appended_xml_files")

unique_accessions = final_df.accession.unique()

uniprot_request(unique_accessions)

for value in unique_accessions:
    df = final_df[final_df['accession'] == value]
    f = open((value + '.xml'), 'r')
    first_line = f.readline()
    if first_line == '':
        f.close()
        continue
    else:
        print(value + '.xml')
        f.close()




###########################