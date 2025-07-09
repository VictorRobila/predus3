#!/usr/bin/env python3

import os
import sys
import subprocess
import Settings_Predus

# ==========================
# Equivalent to "use Settings"
# (Assume you have a settings.py with relevant variables)
# ==========================

#import settings  # If your Settings.pm was converted to settings.py

# ==========================
# Determine root directory
# ==========================

# our $PUDGEROOT = '/ifs/home/c2b2/bh_lab/shares/pudge';
#
#cwd = os.getcwd()

#sys.path.append(os.path.join(cwd, "scr"))

#root_dir = os.environ.get("PUDGE", os.environ.get("PUDGEROOT"))
#if root_dir is None:
#    print("Error: PUDGE or PUDGEROOT environment variable must be set.")
#    sys.exit(1)

# ==========================
# Load pudgeUtil module
# (Assume it is converted to pudgeUtil.py in your scr folder)
# ==========================

# pudgeUtil = /groups/bh6_gp/home/shares/pudge/scr/pudgeUtil.pm

#sys.path.append(os.path.join(root_dir, "scr"))
#import pudgeUtil

# ==========================
# Usage info
# ==========================

PGM = os.path.basename(sys.argv[0])
usage = f"""
USAGE:
    {PGM} pdbid similarity

    Get the representative of a particular structure in the given database
"""

if len(sys.argv) != 3:
    print(usage, file=sys.stderr)
    sys.exit(1)

# ==========================
# Run pipeline data setup
# ==========================

#pudgeUtil.getPipelineData()

# ==========================
# Parse arguments
# ==========================

pdbid = sys.argv[1]
simlev = sys.argv[2]

# ==========================
# Get representative at 100% similarity
# ==========================

rep100_cmd = os.path.join(Settings_Predus.pdg_scr_dir, "getRep.pl") + f" {pdbid} 100"
rep100 = subprocess.getoutput(rep100_cmd).strip()

if not rep100:
    print(f"Representative not found for {pdbid}. Exiting.", file=sys.stderr)
    sys.exit(-1)

# ==========================
# Check similarity file exists
# ==========================

simfile = os.path.join(Settings_Predus.pdg_seqdb_dir, f"chain_templates.{simlev}.nbr")
if not os.path.exists(simfile):
    print("Similarity level must be 100, 90 or 60.", file=sys.stderr)
    sys.exit(-1)

# ==========================
# Grep for rep100 in similarity file
# ==========================

with open(simfile) as f:
    nbrs = [line.strip() for line in f if rep100 in line]

if not nbrs:
    print(f"{pdbid} not found in {simlev}% database. Exiting.", file=sys.stderr)
    sys.exit(-1)

nbrs = nbrs[0].split("\t")

# ==========================
# Load all neighbors at 100%
# ==========================

nbrs100_file = os.path.join(Settings_Predus.pdg_seqdb_dir, "chain_templates.100.nbr")
with open(nbrs100_file) as f:
    nbrs100_lines = [line.strip() for line in f]

# ==========================
# Build list of all neighbors
# ==========================

allnbrs = list(nbrs)
for id in nbrs:
    if not id:
        continue
    temp = [line for line in nbrs100_lines if id in line]
    if temp:
        temp_split = temp[0].split("\t")
        allnbrs.extend(temp_split)

# ==========================
# Sort and deduplicate
# ==========================

allnbrs_sorted = sorted(set(filter(None, allnbrs)))

# ==========================
# Print results
# ==========================

for id in allnbrs_sorted:
    print(id)

