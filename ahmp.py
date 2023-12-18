'''
            ahmp: the antigen-hla modeling pipeline          

========================| DESCRIPTION |=======================
    
    ahmp is a multipurpose utility for creating antigen-      
    hla models and extracting useful data from those
    models for downstream scientific applications.

    ahmp searches a locally available template database to
    identify an HLA-peptide template on which to base the
    model's starting spatial configuration. This search is
    performed using a sequence alignment to both the HLA
    and the peptide of interest. The templated structure
    is used as input for a RosettaScripts protocol that
    performs the following operations:

    1. Thread desired peptide sequence onto template
        <SimpleThreadingMover>
   1b. [optional] Mutate residue(s) to non-canonical AA
        <MutateRes>
    2. Add MHC template structure to the pose
        <DeleteRegionMover>
        <AddChain>
    3. Prepack the structure in preparation for a
       FlexPepDock refinement run.
        <FlexPepDock>

    In addition, ahmp creates all files necessary to run
    FlexPepDocking refinement on the prepacked output
    structure.

    ahmp contains additional functionalities in its
    associated modules. To access these help files, use
    the command:

        | $ python ahmp.py --module <module_name>

    To view a list of modules with help documentation,
    pass the command --module without arguments.

========================| REFERENCES |========================

    Bloodworth N. et al. PLOS ONE; 17(2)e0275759.

===========================| USAGE |==========================

===============| DEPENDENCIES AND INSTALLATION |==============

'''
## IMPORTS ##
# Python standard library
import argparse
# ahmp modules

# BODY OF EXECUTION
def main(args):
    success=False

    if args.model:
        # Initiate modeling protocol
        success=True

    if args.query:
        # Initiate query protocol
        success=True

    return success
## INPUT ARGUMENTS ##
# Argument parser
parser=argparse.ArgumentParser(
        description="ahmp: The Antigen-HLA Modeling Pipeline",
        usage=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
# === Protocols: Specify Behavior ===
parser.add_argument("-m","--model",
        action="store_true",
        help="Create prepacked peptide-HLA models for FlexPepDock refinement"
        )
parser.add_argument("-q","--query",
        actions="store_true",
        help="Query the Immune Epitope Database or local template database"
        )
# === Inputs: Per protocol '--model'  ===
## Mandatory
parser.add_argument("--peptides",
        nargs="*",
        default=[],
        help="Peptides to thread onto HLA templates. Command-line sequence(s) and/or file(s)"
        )
parser.add_argument("--allele",
        nargs=1,
        default=[],
        help="HLA allele. For human alleles, use the format <GENE>*<ALLELE_GROUP>:<PROTEIN>. For example: A*02:01"
        )
## Optional
parser.add_argument("--rosetta_options",
        nargs="*",
        help="[OPTIONAL] Rosetta FlexPepDocking options. Flag/value pairs are formatted as flag=value. Pass 'help' to view defaults"
        )
parser.add_argument("--slurm",
        nargs="*",
        help="[OPTIONAL] SLURM batch #SBATCH options. Flag/value pairs are formatted as flag=value. Pass 'help' to view defaults."
        )
parser.add_argument("--rosetta_main",
        nargs=1,
        default=["/dors/meilerlab/apps/rosetta/rosetta-3.13/main"],
        help="[OPTIONAL] Path to Rosetta/main"
        )
parser.add_argument("--params",
        nargs=2,
        default=["","0"],
        help="[OPTIONAL] Rosetta .params file for a non-canonical AA, and optional position at which to make the NCAA substitution"
        )
parser.add_argument("--hla_template",
        nargs=1,
        default=[""],
        help="[OPTIONAL] PDB to use as the HLA template (rather than automatically selecting from available structures)"
        )
parser.add_argument("--verbose",
        action="store_true",
        help="[OPTIONAL] Toggle detailed command-line feedback"
        )
parser.add_argument("--ignore_epitope_match",
        action="store_true",
        help="[OPTIONAL] [BENCHMARKING ONLY] Ignore exact sequence matches when selecting a peptide structural template"
        )
parser.add_argument("--worst_epitope_template",
        action="store_true",
        help="[OPTIONAL] [BENCHMARKING ONLY] Select the peptide structural template by worst alignment score"
        )
# === Inputs: Per protocol '--query' ===
## Mandatory
## Optional

# Parse and execute
args=parser.parse_args()

if __name__=="__main__":
    succes=main(args)
    if success:

