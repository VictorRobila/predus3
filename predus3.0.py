#!/usr/bin/env python3
# -----------------------------------------------------------------------------#
#  PredUs 2.0  –  Pure-Python front-end (driver)                                #
#  --------------------------------------------------------------------------- #
#  This file is a **direct port** of the original Perl script predus2.0.pl      #
#  used in the Honig-lab PredUs / PrePPI pipeline.  All the heavy science       #
#  (structure alignment, neighbour compilation, interface extraction,          #
#  SVM scoring, empirical reweighting) is farmed out to **placeholder_*.py**    #
#  modules so the logic can be filled-in incrementally.                         #
#                                                                              #
#  Road-map of the pipeline as reproduced here:                                #
#                                                                              #
#      ┌──────────────┐   1   ┌────────────┐   2   ┌────────────┐              #
#      │   SKAN       │──────▶│ neighbour  │──────▶│  rotpro    │              #
#      │  alignment   │       │  table     │       │ (map merge)│              #
#      └──────────────┘       └────────────┘       └────────────┘              #
#            |                                    │                            #
#            |  3 (nbr file)                      │ 4 ( .out )                 #
#            ▼                                    ▼                            #
#      ┌──────────────┐   5   ┌────────────┐   6   ┌───────────────┐           #
#      │ interface.pl │──────▶│  SVM       │──────▶│ predus2_emp   │           #
#      │  residue list│       │ rescoring  │       │  final scores │           #
#      └──────────────┘       └────────────┘       └───────────────┘           #
#                                                                              #
#  Only steps 1-3 are executed in `run_predus()`.  Steps 4-6 happen in `main`. #
#  Each numbered rectangle corresponds to a placeholder module you will        #
#  implement later.                                                            #
# -----------------------------------------------------------------------------#

"""
WAS predus2.0.pl  –  Python translation of the original Perl PredUs front-end.

Placeholder modules you must eventually implement 1-for-1:
    ▸ placeholder_skan.run_skan()               (step-1)
    ▸ placeholder_compileNbrTable.compile_nbr_table()   (step-2)
    ▸ placeholder_rotpro.rotpro()               (step-3/4)
    ▸ placeholder_interface.extract_interface() (step-5)
    ▸ placeholder_svmPredict.svm_predict()      (step-6a)
    ▸ placeholder_predus2_emp.combine_empirical()       (step-6b)
"""

# -------------------------- std-lib imports ---------------------------------#
import argparse
import shutil
import sys
import tempfile
from pathlib import Path

# ------------------------- placeholder imports ------------------------------#
# each import MUST exist (even if empty) so that the driver runs.             #
# replace the `raise NotImplementedError` lines inside them as you port code. #
import placeholder_skan                 as skan_mod          # step-1
import placeholder_compileNbrTable      as nbr_mod           # step-2
import placeholder_rotpro               as rotpro_mod        # step-3/4
import placeholder_interface            as iface_mod         # step-5
import placeholder_svmPredict           as svm_mod           # step-6a
import placeholder_predus2_emp          as emp_mod           # step-6b

# -------------------------- static settings -------------------------------- #
# In the original Perl these came from Settings.pm
PREDUSROOT    = Path("/ifs/data/c2b2/bh_lab/shares/hfpd/web/PredUs_new")
PREDUSSCRATCH = Path("/ifs/scratch/c2b2/bh_lab/shares")
PSD_CUTOFF    = 0.6     # identical numeric threshold to Perl

# --------------------------------------------------------------------------- #
# utility: gated stderr print (Verbosity flag)                                #
def vprint(flag: bool, *msg):
    """print to STDERR only if flag is true (mimics Perl $verbose)"""
    if flag:
        print(*msg, file=sys.stderr)

# --------------------------------------------------------------------------- #
# ----------------------------- core subroutine ----------------------------- #
# direct Python replacement of Perl runPredUS()                               #
# --------------------------------------------------------------------------- #
def run_predus(struct_id: str,
               pdb_file: Path,
               skan_override: Path | None,
               wrk: Path,
               verbose=False) -> list[str]:
    """
    • Ensures SKAN output exists   (.08.skan.fa)
    • Builds neighbour table       (.nbr)
    • Calls rotpro & interface.pl  (.out → residue list)
    • Returns list[str]  of residue IDs predicted by the template-mapping step
    """
    # --- create log/ intf/ sub-dirs inside scratch -------------------------#
    (wrk / "log").mkdir(exist_ok=True)
    (wrk / "intf").mkdir(exist_ok=True)

    # ----------------------------------------------------------------------#
    # 1) SKAN alignment                                                     #
    # ----------------------------------------------------------------------#
    skan_fa = wrk / f"{struct_id}.08.skan.fa"
    if skan_override:
        # user supplied -k ⇒ copy verbatim
        shutil.copy2(skan_override, skan_fa)
        vprint(verbose, f"Using pre-computed SKAN file: {skan_fa}")
    else:
        vprint(verbose, "Running SKAN …")
        skan_mod.run_skan(
            query_pdb=pdb_file,
            skads=PREDUSSCRATCH / "databases" / "skads" / "templates.60.skads",
            out_fasta=skan_fa
        )

    # ----------------------------------------------------------------------#
    # 2) compileNbrTable → parse neighbours                                 #
    # ----------------------------------------------------------------------#
    neighbor_lines = nbr_mod.compile_nbr_table(skan_fa, struct_id)

    # filter neighbours by PSD and stash SAS/RMSD for ranking
    info: dict[str, tuple[float, float, float]] = {}
    for line in neighbor_lines:
        pair, psd, rmsd, sas = line.rstrip().split("\t")
        psd = float(psd)
        if psd > PSD_CUTOFF:
            continue
        neighbour_id = pair.split("-")[1]              # same parse as Perl
        info[neighbour_id] = (float(sas), psd, float(rmsd))

    if not info:
        raise RuntimeError("No structural neighbours within PSD cutoff")

    # write <ID>.nbr exactly like the Perl script (sorted by SAS ascending)
    nbr_file = wrk / f"{struct_id}.nbr"
    with open(nbr_file, "w") as nf:
        for nid in sorted(info, key=lambda k: info[k][0]):
            sas, psd, rmsd = info[nid]
            nf.write(f"{nid}\t{sas}\t{psd}\t{rmsd}\n")
    vprint(verbose, f"Wrote neighbour file: {nbr_file}")

    # ----------------------------------------------------------------------#
    # 3) rotpro + interface.pl                                              #
    # ----------------------------------------------------------------------#
    rot_out = wrk / f"{struct_id}.out"
    rotpro_mod.rotpro(struct_id, pdb_file, nbr_file, rot_out, wrk)

    residue_string = iface_mod.extract_interface(rot_out)
    if not residue_string:
        raise RuntimeError("interface.pl returned empty string")

    # move .out to intf/<ID>.predus.pdb to mirror Perl
    final_pdb = (wrk / "intf" / f"{struct_id}.predus.pdb")
    shutil.move(rot_out, final_pdb)

    return residue_string.split()

# --------------------------------------------------------------------------- #
# ----------------------------- command-line wrap --------------------------- #
# --------------------------------------------------------------------------- #
def main():
    # ------- CLI identical to Perl driver ---------------------------------#
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("-s", required=True, help="structure ID (e.g. 6re0_E)")
    ap.add_argument("-f", required=True, type=Path, help="query PDB file")
    ap.add_argument("-o", required=True, type=Path, help="output .intf/.txt")
    ap.add_argument("-k", type=Path, help="pre-computed SKAN file (.fa)")
    ap.add_argument("-t")                       # optional, ignored
    ap.add_argument("-V", action="store_true")  # verbose
    ap.add_argument("-D", action="store_true")  # debug (keep scratch)
    ap.add_argument("-v", action="store_true")  # SVM-only mode (early exit)
    ap.add_argument("-h", action="help")
    args = ap.parse_args()

    # ------- create scratch directory (Perl tempdir("scratchXXXXXX") ) ----#
    workdir = Path(tempfile.mkdtemp(prefix="scratch"))
    if args.D:  # override to local cwd for easy inspection
        workdir = Path(tempfile.mkdtemp(prefix="scratch", dir="."))
    vprint(args.V or args.D, f"Scratch directory: {workdir}")

    try:
        # ---------------- stage-1: template mapping -----------------------#
        residues = run_predus(args.s, args.f, args.k, workdir, args.V)

        # ---------------- replicate Perl's svm_temp.xx --------------------#
        svm_temp = workdir / "svm_temp.xx"
        with open(svm_temp, "w") as fh:
            fh.write(f"{args.s}\t" + " ".join(residues) + "\n")

        # ---------------- stage-2: SVM rescoring --------------------------#
        svm_test, svm_out = svm_mod.svm_predict(args.s, workdir, residues)

        # optional early exit (`-v` in Perl prints only residues>0)
        if args.v:
            with open(svm_out) as s_out:
                for res, score in zip(residues, s_out):
                    if float(score) > 0:
                        print("PD2:", res)
            return

        # ---------------- stage-3: empirical combiner ---------------------#
        emp_mod.combine_empirical(args.s, svm_temp, workdir / "temp.text")

        shutil.copy2(workdir / "temp.text", args.o)
        print(f"Interface written to {args.o}")

    finally:
        if args.D:
            vprint(True, f"Debug: leaving scratch dir {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
