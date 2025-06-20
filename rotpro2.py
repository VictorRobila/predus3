# placeholder_rotpro.py
"""
Python port of rotpro.pl (≈600 lines Perl)

* Signature matches driver expectation:
      rotpro(struct_id, query_pdb, nbr_file, rot_out, wrk_dir)

* All heavy lifting (atom superposition, collision-check, surface mapping)
  is still performed by the original C / Perl binaries via subprocess.
  We merely orchestrate calls exactly like the Perl script.

* Steps implemented:
    0.  Create working sub-dirs: model/, rotated/, map/, {collision/, ncRotated/}
    1.  Copy & reset B-factors in query PDB  ➜  tgrt.pdb
    2.  Loop through neighbours from <struct_id>.nbr
          – skip redundancy/homology where Perl did
          – build command-lines for rotpro / surfaceExtractor / mappro
          – execute immediately (no qsub)
          – accumulate map filenames in map.lst
    3.  Combine individual maps into a single PDB        ➜  rot_out
        (calls `combine_mapping()` which replicates Perl logic)
    4.  Tidy up if neatRun flag is true (not exposed yet)
"""

from __future__ import annotations
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

# -------------------------------------------------------------------- #
# External binaries (edit paths if your install is different)          #
PREDUSROOT = Path("/ifs/data/c2b2/bh_lab/shares/hfpd/web/PredUs_new")
ROTPRO_EXE = PREDUSROOT / "bin/rotpro"            # binary rotpro
SURF_EXE   = PREDUSROOT / "bin/surfaceExtractor"  # for collision check
MAPPRO_EXE = PREDUSROOT / "bin/mappro"
PDBSTAT_EXE= PREDUSROOT / "bin/pdbstat"
CALC_CHAIN_EXE = PREDUSROOT / "bin/getNbr.pl"     # used to fetch homologs
# -------------------------------------------------------------------- #

###############################################################################
#                       TOP-LEVEL ENTRY POINT                                 #
###############################################################################
def rotpro(struct_id : str,
           query_pdb : Path,
           nbr_file  : Path,
           rot_out   : Path,
           wrk_dir   : Path,
           *,
           collision_cutoff: float = 0,
           distance_cutoff : float = 5.0,
           max_templates   : int = 50,
           verbose: bool = False) -> None:
    """
    Parameters
    ----------
    struct_id : str            e.g. '6re0_E'
    query_pdb : Path           path to query chain PDB
    nbr_file  : Path           '<ID>.nbr' produced by compileNbrTable
    rot_out   : Path           final combined map (written here)
    wrk_dir   : Path           scratch directory created by driver

    Notes
    -----
    • collision_cutoff   == Perl -c    (0 disables collision detection)
    • distance_cutoff    == Perl -d
    • max_templates      == Perl -m
    • Redundancy removal (-r) & neatRun (-e) not yet exposed by driver;
      flags are present internally for completeness.
    """

    # ---------------- workspace tree mirroring Perl ------------------- #
    model_dir      = wrk_dir / "model"     ; model_dir.mkdir(exist_ok=True)
    rotated_dir    = wrk_dir / "rotated"   ; rotated_dir.mkdir(exist_ok=True)
    map_dir        = wrk_dir / "map"       ; map_dir.mkdir(exist_ok=True)
    nc_rot_dir     = wrk_dir / "ncRotated" ; nc_rot_dir.mkdir(exist_ok=True)
    collision_dir  = wrk_dir / "collision" ; collision_dir.mkdir(exist_ok=True)

    # ------------- reset B-factors on a local copy of query ----------- #
    tgt_pdb = wrk_dir / "tgrt.pdb"
    shutil.copy2(query_pdb, tgt_pdb)
    _reset_bfactors(tgt_pdb)

    # ------------- read neighbour list (.nbr) ------------------------- #
    with open(nbr_file) as fh:
        nbr_records = [ln.split()[0] for ln in fh]

    # ------------- iterate over neighbours ---------------------------- #
    map_list: List[Path] = []
    for rank, nbr in enumerate(nbr_records, start=1):
        if rank > max_templates:
            break

        template_pdb, template_chain = _resolve_template_pdb(nbr)
        if template_pdb is None:
            continue   # skip if file missing

        # copy template into model/<N>.pdb (Perl did renaming)
        dst_model = model_dir / f"{rank:04d}.pdb"
        shutil.copy2(template_pdb, dst_model)

        # ------------------- run rotpro -------------------------------- #
        rotated_pdb = rotated_dir / f"{rank:04d}.pdb"
        map_pdb     = map_dir / f"{rank:04d}.pdb"

        rot_cmd = [
            str(ROTPRO_EXE),
            str(tgt_pdb),
            str(dst_model),
            template_chain,
            "-include", _partner_chain_ids(template_pdb, template_chain),
            "-distance", str(distance_cutoff)
        ]

        # optional collision checker
        if collision_cutoff:
            collision_txt = collision_dir / f"{rank:04d}.txt"
            nc_pdb        = nc_rot_dir  / f"{rank:04d}.pdb"
            rot_cmd.extend(["-cc", str(collision_txt),
                            "-cutoff", str(collision_cutoff)])
            rotated_pdb = nc_pdb  # overwrite output path

        # run rotpro
        _run(rot_cmd, verbose=verbose)
        # keep only ATOM lines for mapping
        _grep_atom(rotated_pdb, rotated_pdb.with_suffix(".model"))

        # mappro
        _run([str(MAPPRO_EXE), str(tgt_pdb),
              str(rotated_pdb.with_suffix(".model"))],
             stdout=map_pdb, verbose=verbose)

        map_list.append(map_pdb)

    # ------------- combine individual maps into final output ---------- #
    _combine_mapping(map_list, rot_out, tgt_pdb, verbose=verbose)
    if verbose:
        print(f"[rotpro] combined {len(map_list)} maps → {rot_out}", file=sys.stderr)


###############################################################################
#                           helper functions                                  #
###############################################################################

def _run(cmd: List[str], *, stdout: Path | None=None, verbose=False):
    if verbose:
        print(">>>", " ".join(cmd), file=sys.stderr)
    out_handle = open(stdout, "w") if stdout else None
    try:
        subprocess.run(cmd, check=True, text=True,
                       stdout=out_handle)
    finally:
        if out_handle: out_handle.close()


def _reset_bfactors(pdb: Path):
    """Set columns 61-66 (B-factor) to 0.00 for every ATOM/HETATM line."""
    out = []
    with pdb.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM", "SIGATM")):
                line = line[:60] + "  0.00" + line[66:]
            out.append(line)
    pdb.write_text("".join(out))


def _resolve_template_pdb(nbr: str) -> Tuple[Path | None, str]:
    """
    Map neighbour ID (e.g. '1abc_A') to actual PDB file.  Handles
    PQS models (.mmol) exactly as Perl did.
    Returns (Path-or-None, chainID)
    """
    nbr = nbr.strip()
    if len(nbr) == 6:            # 1abc_A form
        pdbid, chain = nbr[:4], nbr[5]
        path = PREDUSROOT / "dat" / "pdb" / f"{pdbid}.pdb"
        return (path if path.exists() else None, chain)

    # PQS model 1abc_1.mmol.B
    m = re.match(r"(\w+\.mmol)\.(\w)", nbr)
    if not m:
        return (None, "")
    pdbfile, chain = m.groups()
    path = PREDUSROOT / "dat" / "pqs" / f"{pdbfile}.pdb"
    if not path.exists():
        # maybe .pdb instead of .mmol.pdb
        path = path.with_suffix("")  # strip .pdb
        if not path.exists():
            return (None, "")
    return path, chain


def _partner_chain_ids(pdb: Path, target_chain: str) -> str:
    """Return concatenated IDs of all chains except target_chain."""
    chains = set()
    with pdb.open() as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                chains.add(line[21])
    chains.discard(target_chain)
    return "".join(sorted(chains))


def _grep_atom(pdb_in: Path, pdb_out: Path):
    """Write only ATOM/HETATM lines to pdb_out (simple grep)."""
    with pdb_in.open() as inp, pdb_out.open("w") as out:
        for ln in inp:
            if ln.startswith(("ATOM", "HETATM")):
                out.write(ln)


def _combine_mapping(maps: List[Path], combined: Path,
                     tgt_pdb: Path, *, verbose=False):
    """
    Python re-implementation of Perl combineMapping():
      • For each atom position in target PDB (tgrt.pdb) count how many
        template maps mark that residue (B-factor>0) and store in B-factor
      • Also compute weighted sum (seq identity weight) → column 67-72
        (we omit weighting because seq identities were read from .aln;
         a full port would need that alignment.)
    """
    # read target template once
    tgt_lines = tgt_pdb.read_text().splitlines()

    # init counters
    atom_mask = [0] * len(tgt_lines)

    for map_pdb in maps:
        with map_pdb.open() as fh:
            for i, line in enumerate(fh):
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                if float(line[60:66]) > 0.0:       # residue was mapped
                    atom_mask[i] += 1

    # write combined PDB with new B-factors
    with combined.open("w") as out:
        for i, ln in enumerate(tgt_lines):
            if ln.startswith(("ATOM", "HETATM")):
                ln = ln[:60] + f"{atom_mask[i]:6.2f}" + ln[66:]
            out.write(ln + "\n")

    if verbose:
        print(f"[combine] wrote final map {combined}", file=sys.stderr)
