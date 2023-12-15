# ahmp
**A**ntigen-**H**LA **M**odeling **P**ipeline

rahmp is a Rosetta-based modeling pipeline that creates templated peptide/HLA models for use as starting structures in Rosetta FlexPepDocking
refinement decoy generation.

## Dependencies
### Python v3.11.3
- Python 3.11.3
- argparse
- sys
- os
- multiprocessing
- math
- pathlib.Path
- zipfile.ZipFile
- shutil
- re
- xml.etree.ElementTree
- subprocess
### Biopython v1.81
- SeqIO
- Seq
- SeqRecord
- pairwise2
- AlignIO
- Align.substitution_matrices
- PDB.PDBParser
- PDB.PDBIO
### Pandas v2.1
### numpy v1.26
### requests 2.31
### Rosetta v3.13
- Rosetta 3.13
