"""Module to align protein sequences using Clustal Omega"""

import os
import subprocess
from pathlib import Path

from Bio import AlignIO, SeqIO


def get_alignment(input_file: Path | str, out_dir: Path | str) -> AlignIO.MultipleSeqAlignment:
    if not isinstance(input_file, Path):
        input_file = Path(input_file)
    if not isinstance(out_dir, Path):
        out_dir = Path(out_dir)
    if not out_dir.exists():
        os.makedirs(out_dir, exist_ok=True)

    """Align protein sequences using Clustal Omega."""
    records = list(SeqIO.parse(input_file, 'fasta'))

    padded_sequences_path = out_dir / f'{input_file.stem}_padded{input_file.suffix}'
    aligned_fasta_path = out_dir / f'{input_file.stem}_aligned{input_file.suffix}'
    # write to temporary file and do alignment
    with padded_sequences_path.open('w', encoding="utf-8") as f:
        SeqIO.write(records, f, 'fasta')

    if len(records) == 1:
        align = records
    else:
        # TODO: magic path for now, but this won't work with a proper library installation
        #   - maybe do relative to __file__?
        clustal_omega_path = '/home/hendraet/stud_sync/Studium/phd/proteomics/PROTzilla/backend/ptm-visualization/clustal-omega/clustalo-1.2.4-Ubuntu-x86_64'
        cmd = f"{clustal_omega_path} \
                --infile={padded_sequences_path} \
                --outfile={aligned_fasta_path} \
                --outfmt=fasta \
                --iter=0 \
                --force"
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True, check=True)
        align = AlignIO.read(f"{aligned_fasta_path}", "fasta")

    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / 'aligned.fasta').open('w', encoding="utf-8") as f:
        SeqIO.write(align, f, 'fasta')

    # clean up temporary files
    if padded_sequences_path.exists():
        padded_sequences_path.unlink()
    if aligned_fasta_path.exists():
        aligned_fasta_path.unlink()

    return align
