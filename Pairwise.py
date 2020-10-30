import numpy as np

def main():
    paired = PairwiseAlignment("AGC", "AAAC")
    print(paired)
    # print(scr[:,:,0])
    # print(scr[:,:,1])

# Pairwise alignment function
def PairwiseAlignment(dna1, dna2):
    scores = np.zeros(shape=(len(dna1) + 1, len(dna2) + 1, 2))
    for x in range(0, len(dna1) + 1):
        for y in range(0, len(dna2) + 1):
            if x == 0:
                if y == 0:
                    scores[x, y, 0] = 0
                    scores[x, y, 1] = 0
                else:
                    scores[x, y, 0] = scores[x, y-1, 0] - 2
                    scores[x, y, 1] = 1
            else:
                if y == 0:
                    scores[x, y, 0] = scores[x-1, y, 0] - 2
                    scores[x, y, 1] = 2
                else:
                    upscore = scores[x, y-1, 0] - 2
                    leftscore = scores[x-1, y, 0] - 2
                    
                    score = 0
                    if dna1[x-1] == dna2[y-1]:
                        score = 1
                    else:
                        score = -1
                        
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

    return ["".join(strand1), "".join(strand2)]

# Define entry point
if __name__ == "__main__":
    main()