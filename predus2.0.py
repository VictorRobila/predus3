#!/usr/bin/env python3

import os
import sys
import argparse
import tempfile
import shutil
import subprocess
from datetime import datetime
import glob

def main():
    PREDUSROOT = "/path/to/predusroot"
    PREDUSCRATCH = "/path/to/predusscratch" # we might not want to use the global perl txt file
    PSDcutoff = 0.6

    PWD = os.getcwd();
    PGM= os.path.basename(sys.argv[0])

    usage_text = f"""
    USAGE:
            {PGM} -s structure_ID -f structure_file -o structure_interface_file -k skan_file
            options:
            -t  server_type
            -V  verbose mode
            -D  debug mode
            -h  this help
            -k  skan file
    """

    parser = argparse.ArgumentParser(description="PredUs script", usage=usage_text)
    parser.add_argument('-s', required=True, help='Structure ID')
    parser.add_argument('-f', required=True, help='Structure file')
    parser.add_argument('-o', required=True, help='Output file for structure interface')
    parser.add_argument('-k', required=True, help='SKAN file')
    parser.add_argument('-t', help='Server type', default=None)
    parser.add_argument('-V', action='store_true', help='Verbose mode')
    parser.add_argument('-D', action='store_true', help='Debug mode')
    args = parser.parse_args()

    verbose = 1 if args.V else 0
    debug = 1 if args.D else 0
    if debug: 
        verbose = 1

    inputStructLstF = args.s
    structFileLstF = args.f
    predus_emp_input_str = os.path.basename(structFileLstF)
    predus_emp_temp_input = structFileLstF

    result_file_location = "./"
    usesvm = 'yes' if args.V else 'no'
    if args.f and not (structFileLstF.startswith('/') or structFileLstF.startswith('~')):
        structFileLstF = os.path.join(PWD, structFileLstF)

    outputIntfF = args.o
    predus_svm_output = "./"
    svm_outfile_location = "./"

    if args.o and not (outputIntfF.startswith('/') or outputIntfF.startswith('~')):
        ouptutIntfF = os.path.join(PWD, outputIntfF)

    skanFile="";
    if args.k is not None:
        skanFile=args.k

    clean = True
    if debug:
        clean=False

    wrkDir = tempfile.mkdtemp(prefix="scratch", dir = PWD) # this may be in conflict with the perl code, but the line $wrkDir = "$PWD/$wrkDir"; seems to be in error 
    try:
        os.chmod(wrkDir,0o777)
        print(f"chmod 666 {wrkDir}", file=sys.stderr)
    except Exception as e:
        print(f"Failed chmod on {wrkDir}: {e}", file = sys.stderr)


    if clean:
        import atexit
        atexit.register(shutil.rmtree, wrkDir)

    print(f"{PGM} started...")
    print(f"\tTime: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    if skanFile != "":
        try:
            shutil.copy(skanFile,wrkDir)
            print(f"Copied {skanFile} to {wrkDir}", file=sys.stderr)
        except Exception as e:
            print(f"Failed to copy {skanFile} to {wrkDir}: {e}", file=sys.stderr)

    os.chdir(wrkDir)

    # Ensure 'log' directory exists and set permissions
    if not os.path.exists("log"):
        os.mkdir("log")
        try:
            os.chmod("log", 0o777)
            print(f"chmod 777 log", file=sys.stderr)
        except Exception as e:
            print(f"Failed chmod on log: {e}", file=sys.stderr)

    # Ensure 'intf' directory exists and set permissions
    if not os.path.exists("intf"):
        os.mkdir("intf")
        try:
            os.chmod("intf", 0o777)
            print(f"chmod 777 intf", file=sys.stderr)
        except Exception as e:
            print(f"Failed chmod on intf: {e}", file=sys.stderr)



    structID = args.s
    intfRes = runPredUS (structID, structFileLstF) # defined at the bottom
    # Print to stderr
    sys.stderr.write("**** SVM\n\n")

    # List files matching 'map' (as ls -l . map would do)
    try:
        ls_output = subprocess.check_output("ls -l . map", shell=True, stderr=subprocess.STDOUT)
        sys.stderr.write(ls_output.decode())
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.output.decode())
    

    temp_svm_out = os.path.join(PWD, "svm_temp.xx")
    with open(temp_svm_out, "w") as OUTPUT:
        OUTPUT.write(f"{struct_id}\t{' '.join(intf_res)}\n")

    # The line below was originally in the perl code. However, there are no loops. Pasting here just in case. 
    # next if $intfRes[0] eq "not predicted";
    # ---------------------------------------------------------------------
    SvmPredict.predict(structID,workDir, *intfRes) #TODO: svmPredict has not yet been written
    
    sys.stderr.write("**** SVM\n\n")

    # List files matching 'map' (as ls -l . map would do)
    try:
        ls_output = subprocess.check_output("ls -l . map", shell=True, stderr=subprocess.STDOUT)
        sys.stderr.write(ls_output.decode())
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.output.decode())
    

    temp_svm_out = os.path.join(PWD, "svm_temp.xx")
    with open(temp_svm_out, "w") as OUTPUT:
        OUTPUT.write(f"{struct_id}\t{' '.join(intf_res)}\n")



    # Build file paths
    svmTestFile = os.path.join(wrk_dir, f"{structID}.svm.t")
    svmOutFile = os.path.join(wrk_dir, f"{structID}.svm.t.out")

    # Empty list (equivalent to @svmpred=())
    svmpred = []

    # Try to open test file
    try:
        svmt = open(svmTestFile, "r")
    except IOError as e:
        sys.stderr.write(f"Can not open the svm test file!\n{e}\n")
        sys.exit(1)

    # Try to open out file
    try:
        svmo = open(svmOutFile, "r")
    except IOError as e:
        sys.stderr.write(f"Can not open the svm out file!\n{e}\n")
        sys.exit(1)

    # Example of writing SVM_value line
    # Assuming OUTPUT is some file object you opened earlier:
    OUTPUT.write(f"SVM_value:{struct_id}\t")


    for tline, oline in zip(svmt, svmo):
        terms = tline.strip().split("\t")
        resid = terms[0]
        svmvalue = oline.strip()
        OUTPUT.write(f"{resid}:{svmvalue}\t")
        
        try:
            if float(svmvalue) > 0:
                svmpred.append(resid)
        except ValueError:
            # handle case where svmvalue is not a number
            pass

    if usesvm == 'yes':
        for resid in svmpred:
            print(f"PD2: {resid}")
        sys.exit(0)

    OUTPUT.write("\n")
    svmt.close()
    svmo.close()
    OUTPUT.close()
    
    os.chdir(svm_outfile_location)

    # Print working directory and file info to stderr
    sys.stderr.write("Here pwd:\n")
    sys.stderr.write(f"{os.getcwd()}\n")
    sys.stderr.write(f"{predus_emp_temp_input}\n")

    # Copy predus_emp_temp_input to current directory
    try:
        shutil.copy(predus_emp_temp_input, "./")
    except IOError as e:
        sys.stderr.write(f"Error copying {predus_emp_temp_input}: {e}\n")

    # Copy temp_svm_out to ./svm_out.temp
    try:
        shutil.copy(temp_svm_out, "./svm_out.temp")
    except IOError as e:
        sys.stderr.write(f"Error copying {temp_svm_out}: {e}\n")

    cmd = f"//ifs/data/c2b2/bh_lab/shares/bin/predus2_emp.exe {predus_emp_input_str} svm_out.temp"
    subprocess.run(cmd, shell=True, check=True)

    # Log header
    sys.stderr.write("**** predus2.0\n\n")

    # List files (like `ls -l . map`)
    try:
        ls_output = subprocess.check_output("ls -l . map", shell=True, stderr=subprocess.STDOUT)
        sys.stderr.write(ls_output.decode())
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.output.decode())

    # Copy temp.txt to options_o if defined
    if args.o is not None:
        try:
            shutil.copy("temp.txt", args.o)
        except IOError as e:
            sys.stderr.write(f"Error copying temp.txt to {args.o}: {e}\n")

    # Read temp.text
    try:
        with open("temp.text", "r") as f:
            predifc = f.readlines()
    except IOError as e:
        sys.stderr.write(f"Error reading temp.text: {e}\n")
        predifc = []

    # Process each line
    for line in predifc:
        line = line.strip()
        data = line.split("\t")
        if len(data) > 1:
            try:
                if float(data[1]) > 0:
                    print(f"PD2: {data[0]}")
            except ValueError:
                # Handle non-numeric data[1] gracefully
                pass
    try:
        shutil.copy("svm_format.cutoff", predus_svm_output)
    except IOError as e:
        sys.stderr.write(f"Error copying svm_format.cutoff to {predus_svm_output}: {e}\n")

    # Copy PD2.* files to result_file_location
    for file in glob.glob("PD2.*"):
        try:
            shutil.copy(file, result_file_location)
        except IOError as e:
            sys.stderr.write(f"Error copying {file} to {result_file_location}: {e}\n")

    # Change directory back to PWD
    os.chdir(PWD)

    # Change file and directory permissions
    try:
        os.chmod(outputIntfF, 0o777)
    except Exception as e:
        sys.stderr.write(f"Error changing permissions on {outputIntfF}: {e}\n")

    # Apply chmod -R 777 equivalent on wrkDir
    for root, dirs, files in os.walk(wrkDir):
        try:
            os.chmod(root, 0o777)
            for d in dirs:
                os.chmod(os.path.join(root, d), 0o777)
            for f in files:
                os.chmod(os.path.join(root, f), 0o777)
        except Exception as e:
            sys.stderr.write(f"Error changing permissions in {wrkDir}: {e}\n")

    # Move temp.text and log the action
    dst = f"/ifs/scratch/c2b2/bh_lab/shares/tmp/{args.s}.pd2.txt"
    try:
        shutil.move(os.path.join(wrkDir, "temp.text"), dst)
        sys.stderr.write(f"Moved temp.text to {dst}\n")
    except IOError as e:
        sys.stderr.write(f"Error moving temp.text: {e}\n")

def runPredUS(structID, structFileLstF):
    # code 
    structureID = structID
    structureFile = structFileLstF
    logFile = os.path.join(wrkDir, "log", f"{structureID}.PredUS.log")
    if(verbose):
        print(f"\t...PredUS, log file at {logFile}")
        print(f"\t\tTime Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    neighborFile = os.path.join(wrkDir, f"{structureID}.nbr")
    neighbor_SAS={}
    neighbor_PSD={}
    neighbor_RMSD={}
    structID_nbrInfoDB = structureID[:4] + "." + structureID[-1]
    first_char_ascii = ord(structID_nbrInfoDB[0])
    if first_char_ascii < 49 or first_char_ascii > 57:
        structID_nbrInfoDB = structID_nbrInfoDB[:4]
    neighbors = []

    neighborNum = len(neighbors) # should always return 0
    if verbose:
        print(f"structure neighbors: {neighborNum}.", file = sys.stderr)
    if not neighborNum:
        skanTmpFile = os.path.join(wrkDir,f"{structureID}.08.skan.fa")
        #if(not os.path.exists(skanFile)) or (os.path.getsize(skanFile) ==0):
            #TODO: print some sort of error, im not sure what htey want 
    else:
        skanTmpFile = skanFile
    if os.path.exists(skanTmpFile) and os.path.getsize(skanTmpFile) >0:
        with open(logFile,"w") as logf:
            try:
                result - subprocess.run(
                    ["/path/to/compileNbrTable.py",skanTmpFile,structureID],
                    stdout=subprocess.PIPE,
                    stderr=logf,
                    text=True,
                    check=True
                )
                neighbors=result.stdout.strip().splitlines()
            except subprocess.Called.ProcessError as e:
                print(f"Error running compileNbrTable.py: {e}", file=sys.stderr)
        print(f"n_neigh: {len(neighbors)}", file=sys.stderr)

    for neighbor in neighbors:
        fields = neighbor.split('\t')
        structure, neighborID = fields[0].split('-', 1)  # split only once (just like Perl)

        PSD = float(fields[1])
        if PSD > PSDcutoff:
            continue

        neighbor_SAS[neighborID] = fields[3]
        neighbor_PSD[neighborID] = fields[1]
        neighbor_RMSD[neighborID] = fields[2]

    closestNeighbors = sorted(neighbor_SAS.keys(),key=lambda k: float(neighbor_SAS[k]))
    print(f"n_close: {len(closestNeighbors)}", file=sys.stderr)

    # Open neighbor file for writing
    with open(neighborFile, "w") as nbrfile:
        for neighborID in closestNeighbors:
            nbrfile.write(
                f"{neighborID}\t{neighbor_SAS[neighborID]}\t{neighbor_PSD[neighborID]}\t{neighbor_RMSD[neighborID]}\n"
            )         
        # Copy neighbor file to xfer location and list it
    try:
        shutil.copy(neighborFile, "/ifs/scratch/c2b2/bh_lab/shares/xfer") #TODO: change to new directory
        print(f"Copied {neighborFile} to /ifs/scratch/c2b2/bh_lab/shares/xfer", file=sys.stderr)
    except Exception as e:
        print(f"Failed to copy {neighborFile}: {e}", file=sys.stderr)

    try:
        subprocess.run(["ls", "-l", neighborFile], check=True, text=True, stderr=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"ls failed: {e}", file=sys.stderr)

    # Check if neighborFile exists and is not empty
    if not os.path.exists(neighborFile) or os.path.getsize(neighborFile) == 0:
        print("ERROR! PredUS did not generate prediction: no structure neighbors.", file=sys.stderr)
        print(f"\tstructure:\t{structureID}", file=sys.stderr)
        print(f"\tfile:\t{structureFile}", file=sys.stderr)
        print(f"\tlog:\t{logFile}\n", file=sys.stderr)

        if verbose:
            print(f"\t\tTime: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stderr)

        if os.path.exists(neighborFile):
            try:
                os.remove(neighborFile)
                print(f"Removed {neighborFile}", file=sys.stderr)
            except Exception as e:
                print(f"Failed to remove {neighborFile}: {e}", file=sys.stderr)

        if verbose:
            print(f"\t\tTime Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stderr)

        # Signal not predicted
        result = "not predicted from_predus_hh.pl (neighbor not found)"
    else:
        result = "predicted"

    rotOutFile= os.path.join(wrkDir, f"{structureID}.out")
    
    print ("**** before rotpro\n", file=sys.stderr)

    try:
        subprocess.run(["ls", "-l", ".", "map"], check=True, text=True, stderr=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"ls failed: {e}", file=sys.stderr)
    rotpro_cmd = [ # TODO: change predus_ / the rotpro program
        f"{PREDUSROOT}/scr/rotpro.pl",
        "-i", structureID,
        "-qgit add ", structureFile,
        "-n", neighborFile,
        "-o", rotOutFile,
        "-w", wrkDir,
        "-m", "30",
        "-c", "0.5",
        "-e",
        "-r"
    ]
    # Print the command to STDERR like Perl does
    print(' '.join(rotpro_cmd), file=sys.stderr)

    # Run rotpro and redirect both stdout and stderr to logFile
    with open(logFile, "a") as logf:
        try:
            subprocess.run(rotpro_cmd, stdout=logf, stderr=logf, check=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"rotpro.pl failed: {e}", file=sys.stderr)

    print("**** rotpro\n", file=sys.stderr)

    # List files again
    try:
        subprocess.run(["ls", "-l", ".", "map"], check=True, text=True, stderr=sys.stderr)
    except subprocess.CalledProcessError as e:
        print(f"ls failed: {e}", file=sys.stderr)


    if (not os.path.exists(rot_out_file)) or (os.path.getsize(rot_out_file) == 0):
        sys.stderr.write("ERROR! PredUS did not generate prediction: no combined map.\n")
        sys.stderr.write(f"\tstructure:\t{structure_id}\n")
        sys.stderr.write(f"\tfile:\t{structure_file}\n")
        sys.stderr.write(f"\tlog:\t{log_file}\n\n")

        if verbose:
            sys.stderr.write(f"\t\tTime: {datetime.now()}\n")

        if verbose:
            sys.stderr.write(f"\t\tTime Finished: {datetime.now()}\n")

        return "not predicted !! from_predus2.0.py (from rotfile)"

    intf_file = os.path.join(wrk_dir, f"{structure_id}.intf")

    # Log the command
    cmd = f"{predus_root}/scr/interface.py -i {rot_out_file}" # TODO: change interface.pl into interface.py
    with open(os.path.join(wrk_dir, "logFile"), "a") as log_f:
        log_f.write(f"{cmd}\n")

    # Run the command, capture stderr in logFile, get last line of stdout
    try:
        proc = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()

        # Append stderr to log file
        with open(os.path.join(wrk_dir, "logFile"), "a") as log_f:
            log_f.write(stderr.decode())

        # Get last line of output
        intf_res_string = stdout.decode().strip().split("\n")[-1]

    except Exception as e:
        sys.stderr.write(f"ERROR running interface.py: {e}\n")
        return "NULL"

    if not intf_res_string:
        return "NULL"

    # Move rotOutFile
    os.rename(rot_out_file, os.path.join(wrk_dir, "intf", f"{structure_id}.predus.pdb"))

    # Example for removing files (commented out)
    # os.remove(neighbor_file)
    # for f in glob.glob(os.path.join(wrk_dir, "tgrt.*")):
    #     os.remove(f)
    # skan_file = os.path.join(wrk_dir, f"{structure_id}.08.skan")
    # if os.path.exists(skan_file):
    #     os.remove(skan_file)

    if verbose:
        print(f"\t\tTime Finished: {datetime.now()}")

    # Set permissions
    os.chmod(output_intf_f, 0o777)

    # Return list of residues
    intf_residue = intf_res_string.split()
    return intf_residue

if __name__ == "__main__":
    main()
