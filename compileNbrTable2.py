# placeholder_compileNbrTable.py
"""
Python clone of compileNbrTable.pl

Perl original logic
-------------------
1. Run  `cali <skanFile>`   (from $HFPD_BIN) — this re-formats *.08.skan.fa
   into groups of 8 lines per SKAN hit:
      0  >tc...      template header with PSD/RMSD
      1  query-seq
      2  hit-seq
      3  query-align-string  (with ‘-’ for gaps)
      4  hit-align-string    (with ‘-’ for gaps)
      5  ...
      6  ...
      7  (blank / newline)

2. For every hit:
   • Count positions where BOTH query & hit have a residue (no ‘-’) ➜ alignLen
   • Extract from header line:
         P:PSD  (field $2)   e.g. 0.08
         R:RMSD (field $3)   e.g. 1.70
         % identity is ignored
   • Skip hits with RMSD == 0   (same test as Perl’s ‘next if $3==0’)
   • SAS =  (RMSD * 100 / alignLen)   then formatted to 0.1f
   • Print:  "<queryID>-<templateID>\tPSD\tRMSD\tSAS"
"""

from pathlib import Path
import re
import subprocess
from typing import List


# ---------------------------------------------------------------------- #
# CONFIG – where the `cali` executable lives (edit if your path differs) #
HFPD_BIN = Path("/ifs/home/c2b2/bh_lab/petrey/research/software/hfpd/bin")
CALI_EXE = HFPD_BIN / "cali"
# ---------------------------------------------------------------------- #


def _run_cali(skan_file: Path) -> list[str]:
    """Call external 'cali' and capture its stdout as a list of lines."""
    result = subprocess.run(
        [str(CALI_EXE), str(skan_file)],
        check=True,
        text=True,
        capture_output=True
    )
    return result.stdout.splitlines()


def compile_nbr_table(skan_fa: Path, struct_id: str) -> List[str]:
    """
    Parameters
    ----------
    skan_fa   : Path   – the <ID>.08.skan.fa alignment file
    struct_id : str    – query chain ID (e.g. '6re0_E')

    Returns
    -------
    List[str] – each entry already formatted:
                '<queryID>-<templateID>\tPSD\tRMSD\tSAS\n'
    """
    lines = _run_cali(skan_fa)

    out_records: list[str] = []
    idx = 0
    while idx < len(lines):
        # 0-based indices mirror the Perl i, i+3, i+4, i+7 logic
        header      = lines[idx]         # >tc_sse:XXXX ...
        query_align = lines[idx + 3]     # aligned query string
        hit_align   = lines[idx + 7]     # aligned hit string
        idx += 8                         # advance to next block

        # Count positions where BOTH have a residue (no gap)
        align_len = sum(
            1 for q, h in zip(query_align, hit_align) if q != "-" and h != "-"
        )
        if align_len == 0:
            continue

        # Parse the header line for template ID, PSD, RMSD
        # Example header:
        #   >tc_sse:1abc_A (P:0.08 R:1.70 I:0.25 A:115)
        m = re.match(r">tc_sse:(\w+)\s+\(P:([\d.]+)\s+R:([\d.]+)", header)
        if not m:
            continue  # malformed line; skip safely
        template_id, psd_str, rmsd_str = m.groups()
        rmsd_val = float(rmsd_str)
        if rmsd_val == 0.0:              # Perl: next if $3==0
            continue

        # Compute SAS = (RMSD * 100 / alignLen), one decimal place
        sas = f"{rmsd_val * 100.0 / align_len:.1f}"

        out_records.append(
            f"{struct_id}-{template_id}\t{psd_str}\t{rmsd_str}\t{sas}\n"
        )

    return out_records
