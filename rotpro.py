import os
import sys
import argparse
import subprocess
import Settings_Predus
import re


def checkPDBClst(ID1, ID2):
    """
    Check if ID1 and ID2 are in the same cluster mapping.
    Returns 1 if found, 0 otherwise.
    """
    # Reformat IDs
    ID1 = ID1.upper()
    ID1 = f"{ID1[:4]}:{ID1[-1]}"
    ID2 = ID2.upper()
    ID2 = f"{ID2[:4]}:{ID2[-1]}"

    # Define mapping file path
    mapping_file = "/ifs/home/c2b2/bh_lab/qz2126/data/databases/uniprot/idmapping/pdb/uni_pdb.capri.map"

    isTheSame = 0

    # Check if both IDs are on the same line
    try:
        with open(mapping_file) as f:
            for line in f:
                if ID1 in line and ID2 in line:
                    isTheSame = 1
                    break
    except FileNotFoundError:
        print(f"Mapping file {mapping_file} not found.", file=sys.stderr)

    return isTheSame


def getStructure(PDB_chainID, predus_root):
    """
    Returns a list of PQS models for a given PDB chain ID.
    If no PQS models are found, returns [PDB_chainID].
    """
    structureList = []

    pdb2pqs = os.path.join(predus_root, "dat/pdb2pqs.map")

    try:
        with open(pdb2pqs) as f:
            linesPDB2PQS = [
                line.split("\t")[1].strip()
                for line in f
                if line.startswith(f"{PDB_chainID}\t")
            ]
    except FileNotFoundError:
        print(f"Mapping file {pdb2pqs} not found.", file=sys.stderr)
        return [PDB_chainID]

    if linesPDB2PQS:
        PQS_chainIDs = linesPDB2PQS[0].split("|")
        for PQSID_chainID in PQS_chainIDs:
            # Replace .X at end with .mmol.X
            parts = PQSID_chainID.rsplit(".", 1)
            if len(parts) == 2:
                PQSID_chainID = f"{parts[0]}.mmol.{parts[1]}"
            structureList.append(PQSID_chainID)
    else:
        structureList = [PDB_chainID]

    return structureList




def getChains(structureFile, predus_root, timeout=10):
    """
    Returns (chainNumber, chainList) for a given structure file using pdbstat.
    If pdbstat fails or times out, returns (0, "").
    """
    chainNumber = 0
    chainList = ""

    try:
        result = subprocess.run(
            [f"{predus_root}/bin/pdbstat", structureFile],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=True
        )
        pdbStatLines = result.stdout.strip().splitlines()

    except subprocess.TimeoutExpired:
        print(f"pdbstat timed out for {structureFile}", file=sys.stderr)
        return (0, "")
    except subprocess.CalledProcessError as e:
        print(f"pdbstat failed for {structureFile}: {e}", file=sys.stderr)
        return (0, "")
    except FileNotFoundError:
        print(f"pdbstat not found at {predus_root}/bin/pdbstat", file=sys.stderr)
        return (0, "")

    # Parse chain IDs
    for pdbStatLine in pdbStatLines:
        if "Chain " in pdbStatLine and ":" in pdbStatLine:
            parts = pdbStatLine.split()
            try:
                idx = parts.index("Chain")
                chainID = parts[idx + 1].rstrip(":")
            except (ValueError, IndexError):
                continue

            if chainID not in chainList:
                chainNumber += 1
                chainList += chainID

    return (chainNumber, chainList)

def resetBF(pdbFile):
    """
    Resets B-factor (columns 61-66) to 0.00 in ATOM, HETATM, or SIGATM lines of a PDB file.
    Overwrites the input file.
    """
    try:
        with open(pdbFile, "r") as f:
            pdbText = f.readlines()
    except FileNotFoundError:
        print(f"PDB file {pdbFile} not found.", file=sys.stderr)
        return

    with open(pdbFile, "w") as f:
        for line in pdbText:
            if line.startswith(("ATOM", "HETATM", "SIGATM")):
                # Replace columns 61-66 (0-based index 60:66) with "  0.00"
                line = line[:60] + "  0.00" + line[66:]
            f.write(line)


def combineMapping(listFile, mapCombined, seqsim):
    """
    Combines mapping files listed in listFile into mapCombined.
    Updates B-factor fields with counts and weighted sequence similarities.
    """
    # Read map file list
    try:
        with open(listFile) as f:
            mapList = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"List file {listFile} not found.", file=sys.stderr)
        return

    if not mapList:
        print(f"No maps listed in {listFile}.", file=sys.stderr)
        return

    # Read first map to initialize
    try:
        with open(mapList[0]) as f:
            outText = f.readlines()
    except FileNotFoundError:
        print(f"Map file {mapList[0]} not found.", file=sys.stderr)
        return

    outLength = len(outText)
    BFactor = [0] * outLength
    BFactorW = [0.0] * outLength

    mapindex = 0

    for mapFile in mapList:
        try:
            with open(mapFile) as f:
                mapText = f.readlines()
        except FileNotFoundError:
            print(f"Map file {mapFile} not found.", file=sys.stderr)
            continue

        if len(mapText) != outLength:
            print(f"ERROR! map {mapFile} does not have equal length ...skipped", file=sys.stderr)
            continue

        resNum = []
        resBF = []
        resStart = 0
        resCount = -1

        # Initialize first residue
        for idx, line in enumerate(mapText):
            if line.startswith(("ATOM", "HETATM", "SIGATM")):
                resNum.append(line[22:26])
                resBF.append(0)
                resCount = 0
                break

        if resCount == -1:
            continue  # No ATOM lines found

        resEnd = 0

        for idx, mapLine in enumerate(mapText):
            if not mapLine.startswith(("ATOM", "HETATM", "SIGATM")):
                continue

            newResNum = mapLine[22:26]
            atomBF = int(float(mapLine[60:66]))

            if newResNum != resNum[resCount]:
                resEnd = idx

                if resBF[resCount] == 1:
                    for idxRes in range(resStart, resEnd):
                        BFactor[idxRes] += 1
                        BFactorW[idxRes] += seqsim[mapindex]

                resCount += 1
                resNum.append(newResNum)
                resBF.append(1 if atomBF > 0 else 0)
                resStart = idx
            else:
                if atomBF > 0:
                    resBF[resCount] = 1

        # Last residue block
        if resBF[resCount] == 1:
            for idxRes in range(resStart, outLength):
                BFactor[idxRes] += 1
                BFactorW[idxRes] += seqsim[mapindex]

        mapindex += 1

    # Update output lines with BFactors
    with open(mapCombined, "w") as OUT:
        for idx, line in enumerate(outText):
            if line.startswith(("ATOM", "HETATM", "SIGATM")):
                bf_str = f"{BFactor[idx]:6.2f}"
                bfw_str = f"{BFactorW[idx]:6.2f}"
                new_line = line[:60] + bf_str + line[66:66] + bfw_str + "\n"
                OUT.write(new_line)
            else:
                OUT.write(line)

    os.chmod(mapCombined, 0o777)



def main():
    QSUBNUM=1
    PDBSTATTIMEOUT = 60
    PREDUSROOT = "/path/to/our/root"
    PWD = os.getcwd()
    PGM = os.path.basename(sys.argv[0])
    usage = f"""
USAGE:
\t{PGM} -i structID -q qPDB -o outPDB -n FILE [-c #] [-d #] [-w working_directory] [-s] [-r] [-m] [-e] [-V] [-D] [-h] 

This script predicts the protein binding interface for PDB structure qPDB from its structural
neighbors (-nbr file) by superimposing the neighbors with qPDB and mapping their interfacial 
residues onto the surface of qPDB.
"""

    help_text = """
  where
\t-i structID is the query structure ID

\t-q qPDB is the query structure. Could be a PDB ID, or a PDB file. 

\t-o outPDB is the output file, with B-factor field denoting the contacting frequency.

\t-n FILE contains the structural neighbor list of qPDB, usually generated from skan.

\t-c # to turn on check collision and specify the minimum percentage of atoms in protein 
\tpartner not in collision for protein partner to not be rejected by collision detector.

\t-d # distance cut-off to rotated partner chains

\t-s to suppress including homologs

\t-r to remove redundancy of the templates of structural neighbors

\t-e runs neatly

\t-h prints more help information.
"""

    # ## GET ARGUMENTS -----------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Predicts protein binding interface by superimposing structural neighbors onto query structure.",
        epilog=help_text,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Adding flags based on Perl getopts string 'i:q:o:n:t:p:b:c:d:m:w:sreVDh'
    parser.add_argument('-i', required=True, help='Query structure ID')
    parser.add_argument('-q', required=True, help='Query structure (PDB ID or file)')
    parser.add_argument('-o', required=True, help='Output PDB file with B-factor field recording contacting frequencies')
    parser.add_argument('-n', help='File containing structural neighbor list of qPDB, usually generated from skan')
    parser.add_argument('-t', help='Option t (not documented)')
    parser.add_argument('-p', help='Option p (not documented)')
    parser.add_argument('-b', help='Option b (not documented)')
    parser.add_argument('-c', type=float, help='Check collision with minimum non-collision atom percentage')
    parser.add_argument('-d', type=float, help='Distance cutoff to rotated partner chains')
    parser.add_argument('-m', help='Option m (not documented)')
    parser.add_argument('-w', help='Working directory')
    parser.add_argument('-s', action='store_true', help='Suppress including homologs')
    parser.add_argument('-r', action='store_true', help='Remove redundancy of templates')
    parser.add_argument('-e', action='store_true', help='Runs neatly')
    parser.add_argument('-V', action='store_true', help='Option V (not documented)')
    parser.add_argument('-D', action='store_true', help='Option D (not documented)')
    parser.add_argument('-h', action='help', help='Prints more help information and exits')

    args = parser.parse_args()
    # Perl: my $debug = $options{D}? 1 : 0;
    debug = 1 if args.D else 0

    # Perl: my $verbose = $options{V}? 1 : 0;
    verbose = 1 if args.V else 0

    # Perl: $verbose = 1 if ( $debug );
    if debug:
        verbose = 1

    # Perl: my $structID = $options{i};
    structID = args.i  # query structure ID, an ID without path and file extension

    # Perl: my $pdbInput = $options{q};
    pdbInput = args.q  # query PDB

    # Perl: my $pdbOutput = $options{o};
    pdbOutput = args.o  # output PDB file

    # Perl: my $nbrF = $options{n};
    nbrF = args.n  # file with structural neighbors

    # Perl: my $cutoff = $options{c} ? $options{c} : 0;
    cutoff = args.c if args.c is not None else 0

    # Perl: my $distance = $options{d} ? $options{d} : 5;
    distance = args.d if args.d is not None else 5

    # Perl: my $maximumNumber = $options{m} ? $options{m} : 50;
    maximumNumber = args.m if args.m is not None else 50

    # Perl: my $removeRedundancy = $options{r} ? 1 : 0;
    removeRedundancy = 1 if args.r else 0

    # Perl: my $neatRun = $options{e} ? 1 : 0;
    neatRun = 1 if args.e else 0

    # Perl: my $wrkDir = defined ( $options{w} )? $options{w} : $PWD;
    wrkDir = args.w if args.w is not None else PWD

    # Perl: print STDERR `/bin/mkdir $wrkDir` if ( not -e $wrkDir );
    if not os.path.exists(wrkDir):
        os.makedirs(wrkDir)
        if verbose:
            print(f"[VERBOSE] Created working directory: {wrkDir}")

    # Perl: chdir $wrkDir;
    os.chdir(wrkDir)

    # Perl: print STDERR $PGM, " started\n";
    print(f"{PGM} started")

    # Perl: print STDERR "\tTime: ", `/bin/date`;
    from datetime import datetime
    print(f"\tTime: {datetime.now()}")

    # Perl: print STDERR "Running directory: $wrkDir\n\n";
    print(f"Running directory: {wrkDir}\n")

    ####################################
    # Perl: my $finalnbrF = "$wrkDir/$structID".".final.nbr";
    finalnbrF = f"{wrkDir}/{structID}.final.nbr"

    # Perl: open ( OUT, ">$finalnbrF");
    try:
        OUT = open(finalnbrF, "w")
    except IOError as e:
        print(f"ERROR: Cannot open file {finalnbrF} for writing: {e}", file=sys.stderr)
        sys.exit(1)
    ####################################

    # ## BACKUP INPUT FILES AND CREATE TEMPORARY OUTPUT FILE --------------------------------
    # Perl: print STDERR `/bin/mkdir model rotated map`;
    for folder in ["model", "rotated", "map"]:
        if not os.path.exists(folder):
            os.makedirs(folder)
            if verbose:
                print(f"[VERBOSE] Created folder: {folder}")

    # Perl: print STDERR `/bin/mkdir collision ncRotated` if ( $cutoff );
    if cutoff:
        for folder in ["collision", "ncRotated"]:
            if not os.path.exists(folder):
                os.makedirs(folder)
                if verbose:
                    print(f"[VERBOSE] Created folder: {folder}")

    # Perl: print STDERR `/bin/mkdir qsub` if ( not -e "qsub");
    if not os.path.exists("qsub"):
        os.makedirs("qsub")
        if verbose:
            print("[VERBOSE] Created folder: qsub")

    # Perl: print STDERR `/bin/cp $pdbInput $wrkDir/tgrt.pdb`;
    import shutil
    shutil.copy(pdbInput, f"{wrkDir}/tgrt.pdb")
    if verbose:
        print(f"[VERBOSE] Copied {pdbInput} to {wrkDir}/tgrt.pdb")

    # Perl: $pdbInput = "$wrkDir/tgrt.pdb";
    pdbInput = f"{wrkDir}/tgrt.pdb"

    # Perl: resetBF ( $pdbInput );  ## reset B-factor for $pdbInput
    # Placeholder for resetBF function call
    # TODO: Implement resetBF equivalent
    resetBF(pdbInput)

    # ## ROTATE THE BINDING PROTEIN PARTNER AND MAP THE SURFACE ATOMS IN CONTACT ---------------
    # generate job label
    os.chdir("qsub")

    # Perl: my $joblabel = `/bin/date +%w%H%S`;
    # chomp($joblabel);
    from datetime import datetime
    joblabel = datetime.now().strftime("%w%H%S")

    # Perl: $joblabel = "RP" . $joblabel . "." . $structID;
    joblabel = f"RP{joblabel}.{structID}"

    # Perl: my $allNbrHomoModelCount = 0;
    allNbrHomoModelCount = 0

    # Perl: my $qsubCount = int( $allNbrHomoModelCount/$QSUBNUM + 1 );
    qsubCount = int(allNbrHomoModelCount / QSUBNUM + 1)

    # Perl: my $qsubF = "$joblabel.$qsubCount.sh"; 
    qsubF = f"{joblabel}.{qsubCount}.sh"

    # Perl: my $qoutF = "$joblabel.$qsubCount.out"; 
    qoutF = f"{joblabel}.{qsubCount}.out"

    # Perl: my $qerrF = "$joblabel.$qsubCount.err"; 
    qerrF = f"{joblabel}.{qsubCount}.err"

    # ## generate structure neighbor list: structures to be rotated
    # Perl: my $mapListF = "$wrkDir/map.lst";
    mapListF = f"{wrkDir}/map.lst"

    # Perl: my $rotEnough = 0;
    rotEnough = 0

    # Perl: my @allNbrLst;
    allNbrLst = []

    # Perl: my $allNbrCount = 0;
    allNbrCount = 0

    # Perl: my @allNbrHomoModelLst;
    allNbrHomoModelLst = []

    # Perl: my %allNbrHomoModelPartners = ();
    allNbrHomoModelPartners = dict()
    
    allNbrHomoModelLstF = os.path.join(wrkDir, "tgrt.anbl")

    with open(allNbrHomoModelLstF, "w") as NBRLST:

        # Use subprocess to replicate `/bin/cut -f N $nbrF`
        def cut_field(n):
            result = subprocess.run(["cut", "-f", str(n), nbrF], capture_output=True, text=True, check=True)
            return result.stdout.strip().splitlines()

        nbrLst = cut_field(1)
        SASLst = cut_field(2)
        PSDLst = cut_field(3)
        RMSDLst = cut_field(4)

        nbrNum = len(nbrLst)
        count = 0

        for idxNbr in range(nbrNum):
            count += 1
            if rotEnough:
                break

            nbr = nbrLst[idxNbr].strip()
            sas = SASLst[idxNbr].strip()
            psd = PSDLst[idxNbr].strip()
            rmsd = RMSDLst[idxNbr].strip()

            # Perl regex: if ($nbr =~ /^>(\w+) /)
            if nbr.startswith(">") and " " in nbr:
                nbr = nbr.split(" ")[0][1:]

            tmpnbr = nbr
            if len(tmpnbr) == 6:
                tmpnbr = tmpnbr[:4] + tmpnbr[5]

            tmpnbr = tmpnbr.lower()

            if verbose:
                print(f"\nstructural neighbor No.{count}: {nbr}", file=sys.stderr)
            
    # New code block starts here
            try:
               
                result = subprocess.run([os.path.join(Settings_Predus.PREDUSROOT, "getNbr.py"), tmpnbr, "60"], capture_output=True, text=True, check=True)
                nbrHomos = result.stdout.strip().split()
            except subprocess.CalledProcessError as e:
                nbrHomos = []
            
            for nbrHomo in nbrHomos:
                if rotEnough:
                    break
                if not nbrHomo.strip():
                    continue

                nbrHomo = nbrHomo.strip()
                if nbrHomo.startswith("d"):
                    nbrHomo = nbrHomo[1:6]

                if verbose:
                    print(f"    sequence homolog of {nbr}: {nbrHomo}", file=sys.stderr)

                found = nbrHomo in allNbrLst

                if found:
                    if verbose:
                        print(f"        structural neighbor {nbrHomo} has been processed before. ...skipped.", file=sys.stderr)
                    continue
                else:
                    allNbrLst.append(nbrHomo)
                    allNbrCount += 1

                nbrHomoPDBID = nbrHomo[0:4]
                nbrHomoChain = nbrHomo[4:5].upper()

                if nbrHomoChain == "_":
                    continue

                nbrHomoStructID = f"{nbrHomoPDBID}.{nbrHomoChain}"
                nbrHomoModels = getStructure(nbrHomoStructID)

                for nbrHomoModel in nbrHomoModels:
                    if rotEnough:
                        break

                    if verbose:
                        print(f"        {nbrHomoModel}", file=sys.stderr)

                    # Extract structure ID and chain ID using regex
                    match = re.match(r"^(\S+).(\S)$", nbrHomoModel)
                    if match:
                        nbrHomoModelStructID, nbrHomoModelChainID = match.groups()
                    else:
                        continue

                    nbrHomoModelF = None

                    if '.mmol.' in nbrHomoModel:
                        # e.g. 1ti8_1.mmol.B -> /dat/pqs/1ti8_1.mmol.pdb
                        nbrHomoModelF = os.path.join(PREDUSROOT, "dat/pqs", f"{nbrHomoModelStructID}.pdb")
                        if not os.path.exists(nbrHomoModelF):
                            # Try alternative file naming
                            nbrHomoModelF = nbrHomoModelF.replace(".mmol.pdb", ".pdb")
                    else:
                        # e.g. 1ti8.B -> /dat/pdb/1ti8.pdb
                        nbrHomoModelF = os.path.join(PREDUSROOT, "dat/pdb", f"{nbrHomoModelStructID}.pdb")

                    if not os.path.exists(nbrHomoModelF):
                        if verbose:
                            print(f"\tERROR: cannot find structure file: {nbrHomoModelF}. ...skipped.", file=sys.stderr)
                        continue

                    nbrHomoModelChainNum, nbrHomoModelChainIDs = getChains(nbrHomoModelF)

                    if nbrHomoModelChainNum == 0:
                        if verbose:
                            print(f"\tERROR: pdbstat failed structure {nbrHomoModelF}. ...skipped.", file=sys.stderr)
                        continue

                    if nbrHomoModelChainNum == 1:
                        if verbose:
                            print(f"\tStructure {nbrHomoModelF} has only one chain. ...skipped.", file=sys.stderr)
                        continue

                    nbrHomoModelPartnerChainIDs = nbrHomoModelChainIDs.replace(nbrHomoModelChainID, "")
                    nbrHomoModelPartnerNum = len(nbrHomoModelPartnerChainIDs)

                    if removeRedundancy == 1:
                        try:
                            clusterF = os.path.join(PREDUSROOT, "dat/pdb_pqs.pred.clstr")
                            # Equivalent to grep in Perl
                            with open(clusterF) as f:
                                modelCluster = [line.strip() for line in f if nbrHomoModel in line]
                        except FileNotFoundError:
                            modelCluster = []

                        if modelCluster:
                            models = modelCluster[0].split("|")

                            for model in models:
                                if nbrHomoModelPartnerNum == 0:
                                    break

                                for idxAllNbr in range(allNbrHomoModelCount):
                                    if nbrHomoModelPartnerNum == 0:
                                        break

                                    if model == allNbrHomoModelLst[idxAllNbr]:
                                        if debug:
                                            print(f"{nbrHomoModel} is in the same cluster as {model}", file=sys.stderr)

                                        allNbrHomoModelPartnerChainIDs = allNbrHomoModelPartners.get(model, "")

                                        # Additional processing for redundancy removal continues here...
                                        if debug:
                                            print(f"\t{nbrHomoModel} partner chains: {nbrHomoModelPartnerChainIDs}", file=sys.stderr)
                                            print(f"\t{model} partner chains: {allNbrHomoModelPartnerChainIDs}", file=sys.stderr)

                                        allNbrHomoModelPartnerChainIDs_list = list(allNbrHomoModelPartnerChainIDs)
                                        nbrHomoModelPartnerChainIDs_list = list(nbrHomoModelPartnerChainIDs)

                                        # Load entire cluster file once for efficiency if repeatedly used
                                        try:
                                            with open(clusterF) as f:
                                                cluster_lines = f.readlines()
                                        except FileNotFoundError:
                                            cluster_lines = []

                                        for nbrHomoModelPartnerChainID in nbrHomoModelPartnerChainIDs_list.copy():  # use .copy() to modify list during iteration
                                            nbrHomoModelPartner = re.sub(r"\S$", nbrHomoModelPartnerChainID, nbrHomoModel)

                                            if debug:
                                                print(f"\t{nbrHomoModel} partner {nbrHomoModelPartner}", file=sys.stderr)

                                            for allNbrHomoModelParterChainID in allNbrHomoModelPartnerChainIDs_list:
                                                modelPartner = re.sub(r"\S$", allNbrHomoModelParterChainID, model)

                                                if debug:
                                                    print(f"\t{model} partner {modelPartner}", file=sys.stderr)

                                                # Emulate `grep -F "$nbrHomoModelPartner" $clusterF | grep -F "$modelPartner"`
                                                inSameCluster = [
                                                    line for line in cluster_lines
                                                    if nbrHomoModelPartner in line and modelPartner in line
                                                ]

                                                if inSameCluster:
                                                    if debug:
                                                        print(f"\t\t{nbrHomoModel} partner {nbrHomoModelPartner} is in the same cluster as:", file=sys.stderr)
                                                        print(f"\t\t\t{model} partner {modelPartner}.", file=sys.stderr)

                                                    # Remove this chain from partner chain IDs
                                                    nbrHomoModelPartnerChainIDs = nbrHomoModelPartnerChainIDs.replace(nbrHomoModelPartnerChainID, "")
                                                    nbrHomoModelPartnerChainIDs_list.remove(nbrHomoModelPartnerChainID)
                                                    nbrHomoModelPartnerNum -= 1

                                                    if debug:
                                                        print(f"\t\t{nbrHomoModel} partner chains: {nbrHomoModelPartnerChainIDs}", file=sys.stderr)

                                                    break  # exit inner loop if match found

                                            if nbrHomoModelPartnerNum == 0:
                                                break  # exit outer loop if no partners left
                    if (allNbrHomoModelCount == 0) or (nbrHomoModelPartnerNum != 0):
                        if rotEnough:
                            break

                        # Rotate HERE!!
                        # Generate scripts, copy model structure files here, rename it
                        # (Implementation depends on your specific pipeline requirements)
                        # For example:
                        # generate_scripts_and_copy_models()

                        new_qsubCount = int(allNbrHomoModelCount / QSUBNUM + 1)
                        if qsubCount != new_qsubCount:
                            # Print current time (equivalent to `/bin/date`)
                            print(f"\tTime: {datetime.now()}", file=sys.stderr)

                            # Run bash script, capturing stdout and stderr
                            try:
                                subprocess.run(
                                    ["bash", qsubF],
                                    stdout=open(qoutF, "w"),
                                    stderr=open(qerrF, "w"),
                                    check=True
                                )
                            except subprocess.CalledProcessError as e:
                                print(f"\tERROR running {qsubF}: {e}", file=sys.stderr)

                            qsubCount = new_qsubCount
                            qsubF = f"{joblabel}.{qsubCount}.sh"
                            qoutF = f"{joblabel}.{qsubCount}.out"
                            qerrF = f"{joblabel}.{qsubCount}.err"

                        allNbrHomoModelLst.append(nbrHomoModel)
                        allNbrHomoModelPartners[nbrHomoModel] = nbrHomoModelPartnerChainIDs
                        allNbrHomoModelCount += 1

                        ##########################################################################
                        pdbid = nbr
                        # if pdbid.startswith("d"): pdbid = pdbid[1:6]

                        OUT.write(f"{nbrHomoPDBID.upper()}\t{nbrHomoChain}\t{allNbrHomoModelPartners[nbrHomoModel]}\t{pdbid}\t{structID}\t{sas}\t{psd}\t{rmsd}\n")
                        ##########################################################################

                        print(f"all model list: {allNbrHomoModelCount}\t{nbrHomoModel}\t{allNbrHomoModelPartners[nbrHomoModel]}", file=sys.stderr)
                        print(f"\tTime: {datetime.now()}", file=sys.stderr)

                        NBRLST.write(f"{allNbrHomoModelCount}\t{nbrHomoModel}\t{allNbrHomoModelPartners[nbrHomoModel]}\n")

                        tmpNbrHomoModel = f"{allNbrHomoModelCount:04d}"
                        tmpNbrHomoModelF = os.path.join(wrkDir, "model", f"{tmpNbrHomoModel}.pdb")

                        try:
                            shutil.copy(nbrHomoModelF, tmpNbrHomoModelF)
                            print(f"Copied {nbrHomoModelF} to {tmpNbrHomoModelF}", file=sys.stderr)
                        except Exception as e:
                            print(f"Error copying file: {e}", file=sys.stderr)

                        rotatedF = os.path.join(wrkDir, "rotated", f"{tmpNbrHomoModel}.pdb")
                        mapF = os.path.join(wrkDir, "map", f"{tmpNbrHomoModel}.pdb")

                        # Append to map list file
                        with open(mapListF, "a") as f:
                            f.write(f"{mapF}\n")

                        includeOpt = f" -include {allNbrHomoModelPartners[nbrHomoModel]} "
                        distanceOpt = f" -distance {distance}"

                        # Write qsub script
                        with open(qsubF, "w") as QSUB:
                            QSUB.write(f"{PREDUSROOT}/bin/rotpro {pdbInput} {tmpNbrHomoModelF} {nbrHomoModelChainID}{includeOpt}{distanceOpt} > {rotatedF}\n")
                            
                            if cutoff:
                                collisionF = os.path.join(wrkDir, "collision", f"{tmpNbrHomoModel}.txt")
                                collisionOpt = f" -cc {collisionF} -cutoff {cutoff} "
                                rotatedF_nc = os.path.join(wrkDir, "ncRotated", f"{tmpNbrHomoModel}.pdb")
                                
                                QSUB.write(f"{PREDUSROOT}/bin/surfaceExtractor -surfLigTest {pdbInput} {rotatedF_nc} 1.4 0.1 | grep 'non-collision' > {collisionF}\n")
                                QSUB.write(f"{PREDUSROOT}/bin/rotpro {pdbInput} {tmpNbrHomoModelF} {nbrHomoModelChainID}{includeOpt}{distanceOpt}{collisionOpt} > {rotatedF_nc}\n")
                            
                            QSUB.write(f"grep ATOM {rotatedF} > {rotatedF}.model\n")
                            QSUB.write(f"{PREDUSROOT}/bin/mappro {pdbInput} {rotatedF}.model > {mapF}\n")

                        # Check if rotation limit reached
                        if maximumNumber and (allNbrHomoModelCount >= maximumNumber):
                            rotEnough = True
                    else:
                        print(f"\tAll similar templates of the target chain {nbrHomoModel} and all of its partners have been processed.", file=sys.stderr)


    # Assuming qsubF, qoutF, qerrF, NBRLST, wrkDir, mapListF, OUT, finalnbrF,
    # structID, pdbOutput, cutoff, neatRun, PWD, PGM are defined earlier.

    # Run bash script and capture stdout/stderr
    try:
        with open(qoutF, "w") as out_f, open(qerrF, "w") as err_f:
            subprocess.run(["bash", qsubF], stdout=out_f, stderr=err_f, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {qsubF}: {e}", file=sys.stderr)

    # Close NBRLST if it was opened with open()
    NBRLST.close()

    # Set permissions recursively for map directory
    for root, dirs, files in os.walk(os.path.join(wrkDir, "map")):
        for momo in dirs + files:
            os.chmod(os.path.join(root, momo), 0o777)

    # Set permission for mapListF
    os.chmod(mapListF, 0o777)

    ####################################
    OUT.close()
    os.chmod(finalnbrF, 0o777)

    print(f"{wrkDir}/{structID}.08.skan.......{finalnbrF}......{wrkDir}......")

    # Call HomoAlign.nbrAlign equivalent (assuming you have implemented or imported it)
    # Example placeholder:
    # Alignment.HomoAlign.nbrAlign(f"{wrkDir}/{structID}.08.skan", finalnbrF, wrkDir)
    # Replace with actual function call in your pipeline

    finalnbraln = f"{finalnbrF}.aln"
    os.chmod(finalnbraln, 0o777)

    # Parse alignment output
    seqsim = []
    try:
        with open(finalnbraln) as f:
            for line in f:
                terms = line.strip().split("\t")
                sim = terms[2]
                seqsim.append(sim)
    except Exception as e:
        print(f"Cannot open {finalnbraln}: {e}", file=sys.stderr)

    ####################################

    # Combine mapping if map list exists and is not empty
    if os.path.exists(mapListF) and os.path.getsize(mapListF) > 0:
        combineMapping(mapListF, pdbOutput)  # replace with your actual function
    else:
        print("\nNo job running. No structural neighbors?", file=sys.stderr)

    ## CLEAN --------------------------------------------------------------------------------
    if neatRun:
        os.chdir(wrkDir)
        print("Cleaning up model, rotated, qsub directories...", file=sys.stderr)
        for folder in ["model", "rotated", "qsub"]:
            shutil.rmtree(folder, ignore_errors=True)

        if cutoff:
            for folder in ["collision", "ncRotated"]:
                shutil.rmtree(folder, ignore_errors=True)

    # Return to original directory
    os.chdir(PWD)

    print(f"{PGM} ended", file=sys.stderr)
    print(f"\tTime: {datetime.now()}", file=sys.stderr)




if __name__ == "__main__":
    main()




