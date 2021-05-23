import os
import sys
import csv
import pprint
import requests
import pandas as pd


#This drops all of the rows that have an NA in the modifications thus no mods to XML
targetaccession = []
if (len(sys.argv) > 1):
    sysarg = pd.DataFrame(columns= ['Accessions'])
    for x in range(1,len(sys.argv)):
        targetaccession.append(sys.argv[x].strip())
    sysarg['Accessions'] = targetaccession
    datafile1 = pd.read_csv('AllProteinGroups.csv')
    datafile2 = datafile1[datafile1['Protein Accession'].isin(targetaccession)]
    datafile3 = datafile2.loc[datafile2['Modification Info List'].dropna().index]
    datafile4 = datafile2.loc[~datafile2['Modification Info List'].index.isin(datafile2['Modification Info List'].dropna().index)]
else:
    datafile2 = pd.read_csv('AllProteinGroups.csv')
    datafile3 = datafile2.loc[datafile2['Modification Info List'].dropna().index]
    datafile4 = datafile2.loc[~datafile2['Modification Info List'].index.isin(datafile2['Modification Info List'].dropna().index)]

    
#This function will be applied row by row and clean the modification column    
def process_modifications(dataframe_accession, dataframe_column_split):
    df = pd.DataFrame(columns = ['accession', 'aa#', 'modification', 'aa'])
    for item in dataframe_column_split:
        if 'DECOY' in dataframe_accession:
            continue
        elif '|' in dataframe_accession:
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
    
def mod_classify(mod_entry):
    if len(mod_entry)>1:
        typemod = "modified residue"
    else:
        typemod = "sequence variant"
        
    return typemod
          

final_df['modification'] = final_df.apply(lambda x: aa_sub_split(x['modification']),axis=1)

final_df['modification type'] = final_df.apply(lambda x: mod_classify(x['modification']),axis=1)




modsfile = open('Mods.txt').readlines()

# Data Management
holder = []
for row in modsfile[2:]:
    holder.append(row.split('\n'))

modstxtdata = []
for lines in holder:
    modstxtdata.append(lines[0])
    

for x in range(0,len(modstxtdata)):
    modstxtdata[x] = modstxtdata[x].split('   ')

tempcount = 0
tempcountdat = []
for x in range(0,len(modstxtdata)):
    tempfind = '//'
    tempcount = tempcount + 1
    if tempfind in modstxtdata[x]:
        tempcountdat.append(tempcount)

####
tracker = 0
targetmoddata = []
for x in range(0,len(tempcountdat)):
    keepdatatemp = []
    button = 0
    for k in range(tracker,tempcountdat[x]):
        fieldone = 'ID'
        fieldtwo = 'TG'
        fieldthree = 'DR'
        if fieldone in modstxtdata[k]:
            keepdataone = modstxtdata[k][1]

        if fieldtwo in modstxtdata[k]:
            keepdatatwo = modstxtdata[k][1]

        if fieldthree in modstxtdata[k]:
            keepdatathree = modstxtdata[k][1]
            button = 1
        if (button == 0) and (k == tempcountdat[x]-1): 
            keepdatathree = 'N/A'
    
        
    tracker = tempcountdat[x]
    
    keepdatatemp.append(keepdataone)
    keepdatatemp.append(keepdatatwo)
    keepdatatemp.append(keepdatathree)
    targetmoddata.append(keepdatatemp)


# Reference to PTM List
ptmlist = open('ptmlist.txt').readlines()

holdertwo = []
for row in ptmlist[62:9747]:
    holdertwo.append(row.split('\n'))

ptmlistdata = []
for lines in holdertwo:
    ptmlistdata.append(lines[0])

for x in range(0,len(ptmlistdata)):
    ptmlistdata[x] = ptmlistdata[x].split('   ')

tempcounttwo = 0
tempcountdattwo = []
for x in range(0,len(ptmlistdata)):
    tempfindtwo = '//'
    tempcounttwo = tempcounttwo + 1
    if tempfindtwo in ptmlistdata[x]:
        tempcountdattwo.append(tempcounttwo)


trackertwo = 0
targetptmdata = []
for x in range(0,len(tempcountdattwo)):
    keepdatatemptwo = []
    buttontwo = 0
    totaldat = []
    testbutton = 0
    unimodtracker = 'Unimod;'
    
    for k in range(trackertwo,tempcountdattwo[x]-1):

        if unimodtracker in ptmlistdata[k][1]:
            testbutton = 1
            
    
    if testbutton == 1:
        for k in range(trackertwo,tempcountdattwo[x]-1):
            ptmfieldone = 'ID'
            ptmfieldtwo = 'KW'
            ptmfieldthree = 'TG'
            ptmfieldfour = 'Unimod'

            
            if ptmfieldone in ptmlistdata[k]:
                ptmkeepdataone = ptmlistdata[k][1]

            if ptmfieldtwo in ptmlistdata[k]:
                ptmkeepdatatwo = ptmlistdata[k][1]

            if ptmfieldthree in ptmlistdata[k]:
                ptmkeepdatathree = ptmlistdata[k][1]
            
            if ptmfieldfour in ptmlistdata[k][1]:
                ptmkeepdatafour = ptmlistdata[k][1]
            else:
                ptmkeepdatafour = 'N/A'
        
    
        keepdatatemptwo.append(ptmkeepdataone)
        keepdatatemptwo.append(ptmkeepdatatwo)
        keepdatatemptwo.append(ptmkeepdatathree)
        keepdatatemptwo.append(ptmkeepdatafour)
        targetptmdata.append(keepdatatemptwo)
    trackertwo = tempcountdattwo[x]


def modtxtinfo(modtype,target):
    valid = 0
    if len(modtype)>1:
        for x in range(0,len(targetmoddata)):
            if (valid == 0):
                namemod = targetmoddata[x][0]
                locationmod = targetmoddata[x][1]
                unimod = targetmoddata[x][2]
              
                if unimod is None:
                    unimod = 'N/A'
                    valid = 1
                elif 'N6,N6-dimethyllysine' in modtype:
                    unimod = 'Unimod; 36.'
                    valid = 1
                elif 'N6-methyllysine' in modtype:
                    unimod = 'Unimod; 34.'
                    valid = 1
                elif 'Phosphoserine' in modtype:
                    unimod = 'Unimod; 21.'
                    valid = 1
                elif 'N6-acetyllysine' in modtype:
                    unimod = 'Unimod; 1.'
                    valid = 1
                elif 'Omega-N-methylarginine' in modtype:
                    unimod = 'Unimod; 34.'
                    valid = 1 
                elif modtype in namemod:
                    if target in locationmod:
                        valid = 1
    else:
        unimod = 'N/A'
    return unimod

            
final_df['Unimod'] = final_df.apply(lambda x: modtxtinfo(x['modification'],x['aa']),axis=1)
            
            
def ptmlinker(unimod, modification, aa_name):
    valid = 0
    if (len(modification)>1):
        for g in range(0,len(targetptmdata)):
            if (valid == 0):
               
                if unimod in targetptmdata[g][3]:
                    if modification in targetptmdata[g][1]:
                        if aa_name in targetptmdata[g][2]:
                            description = targetptmdata[g][0]
                            valid = 1
                            return description
            
            
final_df['Residue Description'] = final_df.apply(lambda x: ptmlinker(x['Unimod'],x['modification'],x['aa_name']),axis=1)

unique_modified_accessions = final_df['accession'].unique()

counter = 0
xmldata = []

for x in range (0,len(unique_modified_accessions)):
    counter += 1
    target = unique_modified_accessions[x]
    mod_data = []
    for w in range(0,len(final_df['accession'])):
        if (final_df['accession'].iloc[w] == target):
            if (final_df['modification type'].iloc[w] == 'sequence variant'):
                mod_data.append('<feature type="sequence variant" description="XXXX" id="XXXX" evidence="XX">')
        
                originalline = "<original>" + final_df['aa'].iloc[w] + "</original>"
                mod_data.append(originalline)
                
                variationline = "<variation>" + final_df['modification'].iloc[w] + "</variation>"
                mod_data.append(variationline)
                
                mod_data.append('<location>')

                positionline = '<position position="' + final_df['aa#'].iloc[w] + '"/>'
                mod_data.append(positionline)

                mod_data.append('</location>')
                mod_data.append('</feature>')
                
            if (final_df['modification type'].iloc[w] == 'modified residue') and (not final_df['Residue Description'].iloc[w] == None):
                featuretype = '<feature type="modified residue" description="' + final_df['modification'].iloc[w] + '" evidence="XX">'
                mod_data.append(featuretype)

                mod_data.append('<location>')

                modpositionline = '<position position="' + final_df['aa#'].iloc[w] + '"/>'
                mod_data.append(modpositionline)

                mod_data.append('</location>')
                mod_data.append('</feature>')
                
    xmldata.append(mod_data)


uniqueaccessions = pd.DataFrame(columns= ['Accessions'])
uniqueaccessions['Accessions'] = unique_modified_accessions
uniqueaccessions['XML Data'] = xmldata

combinedfile = open('combined_xml_file.xml', 'w')

# Requesting from uniprot
def unique_uniprot_request(accession, xmldata):
    url = 'https://www.uniprot.org/uniprot/%s.xml'  % accession
    result = requests.get(url)
    uniprotxmldata = result.text.splitlines()
     
    fileindex = -1
    counter = 0

    success = '[200]>'
    failure = '[404]>'

    if failure in str(result):
        print("")
         
    elif success in str(result):
        if (len(uniprotxmldata)>5):
            
            for lines in uniprotxmldata[0:]:
                keyword = '<evidence key'
                fileindex = fileindex + 1
                    
                if (counter==0):
                    if keyword in lines:
                        keyindex = fileindex
                        counter = 1         

            finalfile = uniprotxmldata[0:keyindex] + xmldata + uniprotxmldata[keyindex:]
            outputfilename = accession + '.xml'
        
            output = open(outputfilename, "w")

            for x in range(0,len(finalfile)):
                output.write(finalfile[x] + "\n")
                combinedfile.write(finalfile[x] + "\n")
                
            combinedfile.write('\n')
            output.write('\n')
            output.close()
         
         
if os.path.isdir("appended_xml_files") == False:
    os.mkdir("appended_xml_files")
    os.chdir("appended_xml_files")
else:
    os.chdir("appended_xml_files")
    
    
def ordinary_uniprot_request(accession):
    if not 'DECOY' in accession:
        if not '|' in accession:
            xml_file = accession + '.xml'
            f = open(xml_file, 'w')
            url = 'https://www.uniprot.org/uniprot/%s.xml'  % accession
            result = requests.get(url)
    
            success = '[200]>'
            failure = '[404]>'
    
            if failure in result:
                print('')
            elif success in result:
                f.write(result.text)
                combinedfile.write(result.text)
                combinedfile.write('\n')
                f.close()
    
    

uniqueaccessions.apply(lambda x: unique_uniprot_request(x['Accessions'],x['XML Data']),axis=1)
datafile4.apply(lambda x: ordinary_uniprot_request(x['Protein Accession']),axis=1)

combinedfile.close()

print(datafile4['Modification Info List'])








