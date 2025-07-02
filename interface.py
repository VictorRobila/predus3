#! /usr/bin/perl 

## ----------------------------------------------------------------------------
## interface.pl
##   predict protein protein interface residues from the residues that
##   frequently interacting with other proteins among its structure 
##   neighors. The frequencies have been computed and recorded in the 
##   B-factor (BF) field. A logit function ( see below ) is designed to 
##   turn these frequencies into probabilities. One can supply its own 
##   probability cutoffs or use default ones to decide whether a residue 
##   is on interface or not.
##
## input:
##   pdb file with BF field recording contacting frequencies
## output:
##   interface residue list
##
## the function is: 
##   p = f(bf) = 1/(1+exp((-bf+max(BF)/2)/(max(BF)/10)))
##             = 1/(1+exp((5*max(BF)-bf)/max(BF))
##
## original author (perl):    Cliff@Honig Lab
## date:                      Nov 14, 2008
## 
## python Refactor & cleanup: Anthony@Honig Lab
## date:                      Jul 2, 2025
## ----------------------------------------------------------------------------

import argparse
import os
import sys
import math
import sys


PGM = os.path.basename(sys.argv[0])

# Create argument parser
parser = argparse.ArgumentParser(
    prog=PGM,
    usage=f"{PGM} -i PDBfile [-c cutoff] [-f interface] [-h]",
    description="""
where
  PDBfile        input protein structure with BF field recording contacting frequencies
  interface      output interface residue list
  scoringMethod  residue probabilities and bools decided by maximum or average of its
                 atom probabilities 
  cutoffSchema   cutoff applied to normalized scores or naive scores
  cutoff         cutoff of probabilities for making boolean predictions
""",
)

# Define options

#parser = argparse.ArgumentParser(add_help=False)


parser.add_argument('-i', required=True, help='Input protein structure (PDB file)')
parser.add_argument('-f', help='Output interface residue list')
parser.add_argument('-t', help='Cutoff schema')
parser.add_argument('-m', help='Scoring method')
parser.add_argument('-s', help='(corresponding option for -s in Perl)')
parser.add_argument('-c', type=float, default=0.075, help='Cutoff of probabilities for making boolean predictions')
#parser.add_argument('-h', action='help', help='Show this help message and exit')



# Parse args
args = parser.parse_args()

pwd = os.getcwd()

# Build full path if needed
pdbInput = args.i
if not os.path.isabs(pdbInput):
    pdbInput = os.path.join(pwd, pdbInput)

# Extract structID from file name (strip directory + extension)
structID = os.path.basename(pdbInput).split('.')[0]

# Check if input file can be opened
if not os.path.isfile(pdbInput):
    print(f"ERROR! cannot open PDB file {pdbInput}.", file=sys.stderr)
    sys.exit(1)

# Check if output file can be created if -f was provided
if args.f:
    try:
        intf_handle = open(args.f, "w")
        intf_handle.close()  # Just testing if we can write
    except IOError:
        print(f"ERROR! cannot write to file {args.interface}.", file=sys.stderr)
        sys.exit(1)

# Use cutoff value
cutoff = args.c

# OUTPUT TO SCREEN AND OUTPUT FILE (if provided)
print(f"#  cutoff: {cutoff}")
print("#--------------------------------")
print(structID)

if args.f:
    with open(args.f, "a") as intf_handle:
        intf_handle.write(f"#  cutoff: {cutoff}\n")
        intf_handle.write("#--------------------------------\n")
        intf_handle.write(f"{structID}\n")

# COMPUTE RESIDUE CONTACTING FREQUENCIES
# Read all lines from input file
with open(pdbInput, "r") as infile:
    pdbText = infile.readlines()

pdbLength = len(pdbText)

# Initialize variables
maxAtomBF = 0.0          # maximum atom BF
maxResBF_avg = 0.0       # max residue BF by atom averages

bfFieldWth = [0]*pdbLength

resNum = []
resBF_max = []           # residue BF defined by max atom BF
resBF_avg = []           # residue BF defined by avg atom BF

# Placeholder for loop index and current residue number
newResNum = None


for idx in range(pdbLength):
    pdbLine = pdbText[idx]
    if pdbLine.startswith(('ATOM', 'HETATM', 'SIGATM')):
        # Columns 23-26 in PDB format are 22-25 in 0-indexed Python
        resNum.append(pdbLine[22:26].strip())
        resBF_max = [0.0]
        break

resCount=0;
resAtomCount=0;

for idx in range(pdbLength):
    pdbLine=pdbText[idx]
    if not pdbLine.startswith(('ATOM', 'HETATM', 'SIGATM')):
       continue
    newResNum = pdbLine[22:26].strip()
    bfFieldWth[idx]=6
    atomBF = pdbLine[60:60+bfFieldWth[idx]].strip()

    while not atomBF.endswith('.00'):
        bfFieldWth[idx]+=1
        atomBF = pdbLine[60:60+bfFieldWth[idx]].strip()
    
    atomBF_val = float(atomBF)
    if atomBF_val>maxAtomBF:
        maxAtomBF = atomBF_val
    if newResNum != resNum[resCount] or len(resNum)==0:
        if resCount<len(resBF_avg):
            resBF_avg[resCount]/=resAtomCount
            if resBF_avg[resCount] > maxResBF_avg:
                maxResBF_avg = resBF_avg[resCount]


            resNum.append(newResNum)
            resBF_avg.append(atomBF_val)
            resBF_max.append(atomBF_val)
            resAtomCount=1
            resCount = len(resNum)-1
    else:
        resAtomCount +=1
        while resCount >=len(resBF_avg):
            resBF_avg.append(0.0) # account for perl auto resizing arrays
        resBF_avg[resCount] += atomBF_val      
        
        while resCount >= len(resBF_max):
            resBF_max.append(0.0)
        if(atomBF_val>resBF_max[resCount]):
            resBF_max[resCount] = atomBF_val


resBF_avg[resCount]/=resAtomCount
if resBF_avg[resCount] > maxResBF_avg:
    maxResBF_avg = resBF_avg[resCount]

resBF = resBF_max.copy()
maxBF = maxAtomBF

if maxBF < 1:
    print("Warning! Too few structure neighbors, not predicted.", file=sys.stderr)
    # check if INTF needs to be closed
    sys.exit()

resProb = []
for idx in range(resCount + 1):
    prob = 1 / (1 + math.exp((5 * maxBF - 10 * resBF[idx]) / maxBF))
    resProb.append(prob)


intf = ""
resCount=0
for idx in range(pdbLength):
    pdbLine=pdbText[idx]
    if not pdbLine.startswith(('ATOM', 'HETATM', 'SIGATM')):
       continue

    newResNum = pdbLine[22:26].strip()

    if newResNum != resNum[resCount] or len(resNum)==0:
        res=resNum[resCount].replace(" ", "")
        intf += f"{res}:{resProb[resCount]} " # this is 1 more s.f. than the perl representation of the code

        resCount+=1

        if resCount < len(resNum):
            resNum[resCount] = newResNum
        else:
            resNum.append(newResNum)

    if idx == pdbLength-1:
        res=resNum[resCount].replace(" ","")
        intf += f"{res}:{resProb[resCount]} "

        resCount +=1
        if resCount < len(resNum):
            resNum[resCount] = newResNum
        else:
            resNum.append(newResNum)
    #isInt = 1 if resProb[resCount] > cutoff else 0 <-- unsure when used 

if intf.endswith(" "):
    intf = intf[:-1]

print(intf)


if args.f:
    with open(args.f, "a") as intf_handle:
        intf_handle.write(intf+"\n")

