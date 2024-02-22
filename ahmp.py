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
    structure. File structure:
    /<record_name>/
        input.pdb 
                        : Prepacked input for FlexPepDock
                          Refinement
        <record_name>.fasta
                        : Peptide fasta sequence modified
                          to include NCAA if applicable
        <record_name>.slurm
                        : SLURM file with options specified
                          by user
        <record_name>.xml
                        : RosettaScripts protocol used to
                          produce <record_name>_input.pdb
        score.sc        : Rosetta scorefile for
                          <record_name>_input.pdb
                            
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
import os

# ahmp modules
from modules import modeling

# BODY OF EXECUTION
def main(args):
    if args.model:
        allele_list=modeling.create_allele_list(
                args.allele,
                )
        for allele in allele_list:
            modeling.modeling_pipeline(
                    args.peptides,
                    allele,
                    args.rosetta_options,
                    args.slurm,
                    args.batch_size,
                    args.mhcdb_location,
                    args.rosetta_main,
                    args.params,
                    args.hla_template,
                    args.verbosity,
                    args.ignore_epitope_match,
                    args.worst_epitope_template,
                    args.root_dir_prod,
                    args.root_dir_input,
                    args.prod_output,
                    args.multitask
                    )
    #if args.query:
        # Initiate query protocol

    return
## INPUT ARGUMENTS ##
# Argument parser
parser=argparse.ArgumentParser(
        description="AHMP: The Antigen-HLA Modeling Pipeline",
        usage=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
# === Protocols: Specify Behavior ===
parser.add_argument("-m","--model",
        action="store_true",
        help="Create prepacked peptide-HLA models for FlexPepDock refinement"
        )
parser.add_argument("-q","--query",
        action="store_true",
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
        nargs="*",
        default=[],
        help="HLA allele, list of alleles, or file containing allele names onto which peptides will be threaded. For human alleles, use the format <GENE>*<ALLELE_GROUP>:<PROTEIN>. For example: A*02:01"
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
parser.add_argument("--mhcdb_location",
        type=str,
        default=os.getcwd(),
        help="[OPTIONAL] Location of MHC database"
        )
parser.add_argument("--rosetta_main",
        type=str,
        default="/dors/meilerlab/apps/rosetta/rosetta-3.13/main",
        help="[OPTIONAL] Path to Rosetta/main"
        )
parser.add_argument("--params",
        nargs=2,
        default=["","0"],
        help="[OPTIONAL] Rosetta .params file for a non-canonical AA and the peptide position (or one letter AA code) at which to make the NCAA substitution"
        )
parser.add_argument("--hla_template",
        type=str,
        default="",
        help="[OPTIONAL] PDB to use as the HLA template (rather than automatically selecting from available structures)"
        )
parser.add_argument("--verbosity",
        type=str,
        default="all",
        choices=["all","warnings","errors","none"],
        help="[OPTIONAL] Toggle detailed command-line feedback"
        )
parser.add_argument("--ignore_epitope_match",
        action="store_true",
        default=False,
        help="[OPTIONAL] [BENCHMARKING ONLY] Ignore exact sequence matches when selecting a peptide structural template"
        )
parser.add_argument("--worst_epitope_template",
        action="store_true",
        default=False,
        help="[OPTIONAL] [BENCHMARKING ONLY] Select the peptide structural template by worst alignment score"
        )
parser.add_argument("--multitask",
        action="store_true",
        default=False,
        help="[OPTIONAL] Use python's multiprocessing to create inputs simultaneously"
        )
parser.add_argument("--batch_size",
        type=int,
        default=1,
        help="[OPTIONAL] Number of inputs to run simultaneously on supercomputing cluster"
        )
parser.add_argument("--root_dir_prod",
        type=str,
        default=os.getcwd(),
        help="[OPTIONAL] Root directory for the production run on supercomputing cluster"
        )
parser.add_argument("--root_dir_input",
        type=str,
        default=os.getcwd(),
        help="[OPTIONAL] Directory containing input for production run on supercomputing cluster"
        )
parser.add_argument("--prod_output",
        type=str,
        choices=["silent","score"],
        default="score",
        help="[OPTIONAL] Specify Rosetta output type (silent or scorefile) for production run"
        )
# Parse and execute
arguments=parser.parse_args()

if __name__=="__main__":
    main(arguments)

