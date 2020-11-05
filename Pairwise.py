import numpy as np
import re

def main():
    # Open both genome fasta files
    # M-SARS: https://www.ncbi.nlm.nih.gov/nuccore/NC_004718.3?from=26398&to=27063&report=fasta
    # M-COVID: https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?from=26523&to=27191&report=fasta
    sarsfile = open('M-SARS.fasta', 'r')
    covidfile = open('M-COVID.fasta', 'r')
    # Skip the comment line
    sarsfile.readline()
    covidfile.readline()
    # Read files into strings
    msars = sarsfile.read()
    mcovid = covidfile.read()
    # Remove new line characters
    msars = msars.replace('\n', '')
    mcovid = mcovid.replace('\n', '')
    # Pairwise align both M gene sequences

    paired = PairwiseAlignment(msars, mcovid, [1, 0, -1])
    print(paired)

    alignscore = ScoreStrand(paired[0], paired[1])
    print(alignscore)

# Pairwise alignment function
def PairwiseAlignment(dna1, dna2, scoring):
    scores = np.zeros(shape=(len(dna1) + 1, len(dna2) + 1, 2))
    for x in range(0, len(dna1) + 1):
        for y in range(0, len(dna2) + 1):
            if x == 0:
                if y == 0:
                    scores[x, y, 0] = 0
                    scores[x, y, 1] = 0
                else:
                    scores[x, y, 0] = scores[x, y-1, 0] + scoring[2]
                    scores[x, y, 1] = 1
            else:
                if y == 0:
                    scores[x, y, 0] = scores[x-1, y, 0] + scoring[2]
                    scores[x, y, 1] = 2
                else:
                    upscore = scores[x, y-1, 0] + scoring[2]
                    leftscore = scores[x-1, y, 0] - 2
                    
                    score = 0
                    if dna1[x-1] == dna2[y-1]:
                        score = scoring[0]
                    else:
                        score = scoring[1]
                        
                    diagscore = scores[x-1,y-1, 0] + score
                    
                    maxscore = max(upscore, leftscore, diagscore)
                    dir = 0
                    if maxscore == upscore:
                        dir = 1
                    elif maxscore == leftscore:
                        dir = 2
                    else:
                        dir = 3
                        
                    scores[x, y, 0] = maxscore
                    scores[x, y, 1] = dir

    xstart = len(dna1) 
    ystart = len(dna2)
    dna1nuc = []
    dna1nuc[:] = dna1
    dna2nuc = []
    dna2nuc[:] = dna2
    strand1 = []
    strand2 = []

    while xstart > 0:
        if scores[xstart, ystart, 1] == 1:
            strand1.insert(0, "_")
            strand2.insert(0, dna2nuc.pop())
            ystart = ystart - 1
        elif scores[xstart, ystart, 1] == 2:
            strand1.insert(0, dna1nuc.pop())
            strand2.insert(0, "_")
            xstart = xstart - 1
        else:
            strand1.insert(0, dna1nuc.pop())
            strand2.insert(0, dna2nuc.pop())
            xstart = xstart - 1
            ystart = ystart - 1           

    # print(scores[:,:,0])
    # print(scores[:,:,1])

    return ["".join(strand1), "".join(strand2)]

def ScoreStrand(sarsstrand, covidstrand):
    score = {
        "SynCount": 0,
        "NonsynCount": 0,
        "IndelCount": 0
    }
    atgindexes = [m.start() for m in re.finditer('ATG', covidstrand)]
    
    protdict = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        }

    for indx in atgindexes:
        strindx = indx
        while True:
            sarstrip = sarsstrand[strindx:(strindx + 3)]
            covtrip = covidstrand[strindx:(strindx + 3)]

            if "_" not in covtrip:
                if protdict[covtrip] == '_':
                    break

            tripscore = ScoreTriplet(sarstrip, covtrip)
            score["SynCount"] = score["SynCount"] + tripscore["SynCount"]
            score["NonsynCount"] = score["NonsynCount"] + tripscore["NonsynCount"]
            score["IndelCount"] = score["IndelCount"] + tripscore["IndelCount"]
            strindx = strindx + 3

    return score

def ScoreTriplet(sarstriplet, covidtriplet):
    score = {
        "SynCount": 0,
        "NonsynCount": 0,
        "IndelCount": 0
    }
    if "_" in covidtriplet:
        gapindexes = [i for i, ltr in enumerate(covidtriplet) if ltr == "_"]
        if len(gapindexes) == 3:
            score["IndelCount"] = 1
        elif len(gapindexes) == 1:
            score["IndelCount"] = 1
        else:
            if 1 in gapindexes:
                score["IndelCount"] = 1
            else:
                score["IndelCount"] = 2
    elif "_" in sarstriplet:
        gapindexes = [i for i, ltr in enumerate(sarstriplet) if ltr == "_"]
        if len(gapindexes) == 3:
            score["IndelCount"] = 1
        elif len(gapindexes) == 1:
            score["IndelCount"] = 1
        else:
            if 1 in gapindexes:
                score["IndelCount"] = 1
            else:
                score["IndelCount"] = 2
    elif sarstriplet != covidtriplet:
        protdict = { 
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
        } 
        if protdict[sarstriplet] != protdict[covidtriplet]:
            score["SynCount"] = 1
        else:
            score["NonsynCount"] = 1

    return score

# Define entry point
if __name__ == "__main__":
    main()