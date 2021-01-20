#!/usr/bin/env python3

import os
import sys
import csv
import pprint
import requests 



if os.path.isdir("resource_files") == True:
    os.chdir("resource_files")
else:
    print("Error: Resource file folder could not be found.")
    sys.exit()

datafile = open('AllProteinGroups.csv').readlines()

# Modification Dictionary
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


def ptmmodinfo(modtype,target):
    for x in range(0,len(targetmoddata)):
        namemod = targetmoddata[x][0]
        locationmod = targetmoddata[x][1]
        unimod = targetmoddata[x][2]

        if modtype in namemod:
            if target in locationmod:
                return unimod


modinfo = []
accession = []
submoddata = []

for row in datafile[1:]:
    rawdataone = row.split(',info:')
    junkrawdata = row.split(',')

    accession.append(junkrawdata[0])

    placeholder = []

    for row in rawdataone[0:]:
        rawdatatwo = row.split(';')

        for row in rawdatatwo[0:]:
            unformatdata = row.split('"')
            sub = '#aa'
            subtwo = 'on'
           
            for text in unformatdata:
                
                if (sub in text) and (subtwo in text):
                    placeholder.append(text)

    for x in range(len(placeholder[0:])):
        placeholder[x] = placeholder[x] + ']'  

    modinfo.append(placeholder)


numaccessions = len(accession)
nummods = []
subsdata = []
datatwo = []
moddata = []
numsubs = []
numothermods = []


for j in range(0,numaccessions):
    nummod = len(modinfo[j])
    nummods.append(nummod)

for k in range(0,numaccessions):
    numsub = 0
    datafour = []
    datathree = []
    numothermod = 0 
    datafive =[]
    datasix = []
    tempdat = []

    for w in range(0,nummods[k]):
        referencedata = modinfo[k][w]
        search = '->'
        keydata = [0,0]
        keydatatwo = [0,0,0]

        if search in referencedata:
            numsub = numsub + 1
            splitdataone = referencedata.split('[')


            for lines in splitdataone[0:]:
                searchtwo = '->'
                aminosearch = '#aa'

                if aminosearch in lines:
                    datasix.append(lines)

                
                if searchtwo in lines:
                    datatwo.append(lines)

            
            refdat = datatwo[len(datatwo)-1]
            find = '->'
           
            if find in refdat:
                splitrefdat = refdat.split('->')
                keydata[0] = splitrefdat[0]
                keydata[1] = splitrefdat[1][0]
                datathree.append(keydata)
                
                
            datafour = datathree
        
        else:
            numothermod = numothermod + 1
            splitdatatwo = referencedata.split('[')
            keepdata = splitdatatwo[1]

            keydatatwo[2] = splitdatatwo[0].split('#aa')[1] 
    
            targetdata = keepdata.split(' on ')
            modtype = targetdata[0]

            if (modtype=='Fe'):
                modtype = 'Fe[III]'
                affectedelement = 'E'

            else:
                modtype = targetdata[0]
                affectedelement = (targetdata[1].split(']')[0])
         

            keydatatwo[0] = modtype
            keydatatwo[1] = affectedelement
            
            datafive.append(keydatatwo)

    for lines in datasix:
        dat = lines.split('aa')
        tempdat.append(dat[1])
    

    for g in range(0,len(datafour)):
        datafour[g].append(tempdat[g])

            
    subsdata.append(datafour)
    numsubs.append(numsub)

    moddata.append(datafive)
    numothermods.append(numothermod)


validdata = []


for x in range(0,numaccessions):
    keyone = 'Alanine' # A
    keytwo = 'Arginine' # R
    keythree = 'Asparagine' # N
    keyfour = 'Aspartic acid' # D
    keyfive = 'Cysteine' # C
    keysix = 'Glutamic Acid' # E
    keyseven = 'Glutamine' # Q
    keyeight = 'Glycine' # G
    keynine = 'Histidine' # H
    keyten = 'Isoleucine' # I
    keyeleven = 'Leucine' # L
    keytwelve = 'Lysine' # K
    keythirteen = 'Methionine' # M
    keyfourteen = 'Phenylalanine' # F
    keyfifteen = 'Proline' # P
    keysixteen = 'Serine' #S
    keyseventeen = 'Threonine' # T
    keyeighteen = 'Tryptophan' # W
    keynineteen = 'Tyrosine' # Y
    keytwenty = 'Valine' # V
    

    tempdata = []


    if (numothermods[x]>0):
        for q in range(0,numothermods[x]):
            unitracker = ptmmodinfo(moddata[x][q][0],moddata[x][q][1])
            modtracker = moddata[x][q][0]
            unimodsearch = 'Unimod'

            if unitracker == None:
                continue
            
            else:
                if 'A' in moddata[x][q][1]:
                    subsmoddata = keyone
                if 'R' in moddata[x][q][1]:
                    subsmoddata = keytwo
                if 'N' in moddata[x][q][1]:
                    subsmoddata = keythree
                if 'D' in moddata[x][q][1]:
                    subsmoddata = keyfour
                if 'C' in moddata[x][q][1]:
                    subsmoddata = keyfive
                if 'E' in moddata[x][q][1]:
                    subsmoddata = keysix
                if 'Q' in moddata[x][q][1]:
                    subsmoddata = keyseven
                if 'G' in moddata[x][q][1]:
                    subsmoddata = keyeight
                if 'H' in moddata[x][q][1]:
                    subsmoddata = keynine
                if 'I' in moddata[x][q][1]:
                    subsmoddata = keyten
                if 'L' in moddata[x][q][1]:
                    subsmoddata = keyeleven
                if 'K' in moddata[x][q][1]:
                    subsmoddata = keytwelve
                if 'M' in moddata[x][q][1]:
                    subsmoddata = keythirteen
                if 'F' in moddata[x][q][1]:
                    subsmoddata = keyfourteen
                if 'P' in moddata[x][q][1]:
                    subsmoddata = keyfifteen
                if 'S' in moddata[x][q][1]:
                    subsmoddata = keysixteen
                if 'T' in moddata[x][q][1]:
                    subsmoddata = keyseventeen
                if 'W' in moddata[x][q][1]:
                    subsmoddata = keyeighteen
                if 'Y' in moddata[x][q][1]:
                    subsmoddata = keynineteen
                if 'V' in moddata[x][q][1]:
                    subsmoddata = keytwenty

                valid = 0

                for g in range(0,len(targetptmdata)):
                    if (valid == 0):
                        if unitracker in targetptmdata[g][3]:
                            if modtracker in targetptmdata[g][1]:
                                if subsmoddata in targetptmdata[g][2]:
                                    tempdata.append(targetptmdata[g][0])
                                    valid = 1
    
    validdata.append(tempdata)


os.chdir('..')

if os.path.isdir("appended_xml_files") == False:
    os.mkdir("appended_xml_files")
    os.chdir("appended_xml_files")
else:
    os.chdir("appended_xml_files")


if (len(sys.argv) < 2):
    totaloutput = open("combined_file.xml", "w")
    for q in range(0,numaccessions):
        addonsub = []
        addonmod = []

        accessionnum = accession[q]

        url = 'https://www.uniprot.org/uniprot/%s.xml'  % accessionnum
        result = requests.get(url)

        xmldata = result.text.splitlines()

        fileindex = -1
        counter = 0

        for lines in xmldata[0:]:
            keyword = '<evidence key'
            fileindex = fileindex + 1
            
            if (counter==0):
                if keyword in lines:
                    keyindex = fileindex
                    counter = 1
        
        if (numsubs[q]>0):
            for x in range(0,numsubs[q]):
                addonsub.append('<feature type="sequence variant" description="XXXX" id="XXXX" evidence="XX">')
                
                originalline = "<original>" + subsdata[q][x][0] + "</original>"
                addonsub.append(originalline)
                
                variationline = "<variation>" + subsdata[q][x][1] + "</variation>"
                addonsub.append(variationline)
                
                addonsub.append('<location>')

                positionline = '<position position="' + subsdata[q][x][2] + '"/>'
                addonsub.append(positionline)

                addonsub.append('</location>')
                addonsub.append('</feature>')
        
        
        if (len(validdata[q])>0):
            for x in range(0,len(validdata[q])):
                featuretype = '<feature type="modified residue" description="' + validdata[q][x] + '" evidence="XX">'
                addonmod.append(featuretype)

                addonmod.append('<location>')

                modpositionline = '<position position="' + moddata[q][x][2] + '"/>'
                addonmod.append(modpositionline)

                addonmod.append('</location>')
                addonmod.append('</feature>')

        finalfile = xmldata[0:keyindex] + addonsub + addonmod + xmldata[keyindex:]

        outputfilename = accession[q] + '.xml'
        output = open(outputfilename, "w")

        for x in range(0,len(finalfile)):
            output.write(finalfile[x] + "\n")
            totaloutput.write(finalfile[x] + "\n")
        output.close()

        totaloutput.write("\n")

    totaloutput.close()

else:
    totaloutput = open("combined_file.xml", "w")
    for x in range(1,len(sys.argv)):
        targetaccession = sys.argv[x].strip()

        for q in range(0,numaccessions):
            if (targetaccession == accession[q]):
                addonsub = []
                addonmod = []

                accessionnum = accession[q]

                url = 'https://www.uniprot.org/uniprot/%s.xml'  % accessionnum
                result = requests.get(url)

                xmldata = result.text.splitlines()

                fileindex = -1
                counter = 0

                for lines in xmldata[0:]:
                    keyword = '<evidence key'
                    fileindex = fileindex + 1
                    
                    if (counter==0):
                        if keyword in lines:
                            keyindex = fileindex
                            counter = 1
                
                if (numsubs[q]>0):
                    for x in range(0,numsubs[q]):
                        addonsub.append('<feature type="sequence variant" description="XXXX" id="XXXX" evidence="XX">')
                        
                        originalline = "<original>" + subsdata[q][x][0] + "</original>"
                        addonsub.append(originalline)
                        
                        variationline = "<variation>" + subsdata[q][x][1] + "</variation>"
                        addonsub.append(variationline)
                        
                        addonsub.append('<location>')

                        positionline = '<position position="' + subsdata[q][x][2] + '"/>'
                        addonsub.append(positionline)

                        addonsub.append('</location>')
                        addonsub.append('</feature>')
                
                
                if (len(validdata[q])>0):
                    for x in range(0,len(validdata[q])):
                        featuretype = '<feature type="modified residue" description="' + validdata[q][x] + '" evidence="XX">'
                        addonmod.append(featuretype)

                        addonmod.append('<location>')

                        modpositionline = '<position position="' + moddata[q][x][2] + '"/>'
                        addonmod.append(modpositionline)

                        addonmod.append('</location>')
                        addonmod.append('</feature>')

                finalfile = xmldata[0:keyindex] + addonsub + addonmod + xmldata[keyindex:]

                outputfilename = accession[q] + '.xml'
                output = open(outputfilename, "w")

                for x in range(0,len(finalfile)):
                    output.write(finalfile[x] + "\n")
                    totaloutput.write(finalfile[x] + "\n")
                output.close()
                
                totaloutput.write("\n")

    totaloutput.close()
