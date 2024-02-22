'''
ahmp.modules.modeling

ahmp module for generating FlexPepDocking input.

A brief word on directory structure and nomenclature created by
ahmp:
    > root_dir_input: provided by the user, the root directory
    containing the input for Rosetta during the production run
    > root_dir_prod: provided by the user, the root directory
    from which Rosetta will be launched during the production
    run
    > parent_dir: the parent directory containing peptide
    directories with input files; also contains the Rosetta
    params file and PDB rotamers file (if using). Direct child
    of root_dir_input. Named after the HLA allele.
    > peptide_dir: Directory containing input files for
    production run. Child of parent_dir.
    > output: Optional folder set as an argument for 
    -out:path:all in Rosetta, if argument not already 
    specified by the user. Must be created at time of production
    run by user.

    ./root_dir_input/
        parent_dir/
            rosetta.params
            rotamers.pdb
            peptide_dir/
                input.pdb
    
    ./root_dir_prod/
        parent_dir/
            output/

'''
# Python standard libraries
from pathlib import Path
import multiprocessing
import shutil
import os
import xml.etree.ElementTree as ET
import subprocess
import random
import math
import sys

# ahmp
from modules.utilities import Notify
from modules.mhcdb import MHCdatabase

## Biopython
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO

def create_allele_list(alleles):
    '''
    Takes in a list of HLA alleles and/or files containing a list of
    alleles and returns a compiled list to of allels to iterate through
    during processing.

    Arguments:
    > alleles [list]: list of alleles or files containing allele names
    
    Positional returns:
    [0] [list]: List of strings to validate as HLA alleles downstream
    '''
    allele_list=[]
    for item in alleles:
        if Path(item).is_file():
            # Assume each line represents an allele name in the provided file
            with open(item,"r") as f:
                lines=f.readlines()
            allele_list.extend(lines)
        else:
            allele_list.append(item)
    
    if not len(allele_list):
        Notify().error(
                "No HLA alleles found! Use the --allele flag to provide "+
                "HLA alleles or a file containing HLA alleles."
                )
    return allele_list

def create_peptide_list(peptides, messenger):
    '''
    Takes in a list of peptides and .fasta files containing peptides and 
    returns a list of Bio.Seq objects.

    Arguments:
    >peptides [list]: List of sequences and/or files containing peptide sequences
    to model
    >messenger [obj]: instance of utilities.Notify class 

    Positional returns:
    [0] [list]: List of validated Bio.Seq objects 
    '''
    valid_AA_codes=set(
            ["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","X"]
            )
    
    peptide_seq_records=[]
    peptide_seq_records_curated=[]
    peptide_seq_names=[]
    for item in peptides:
        # Check if the item is a file. If so, assume it is of the fasta type
        if Path(item).is_file():
            try:
                peptide_data=SeqIO.parse(item,"fasta")
            except:
                messenger.warning(f"Unable to read peptides from file {item}")
                continue
            # add to sequence records
            for peptide in peptide_data:
                peptide_seq_records.append(peptide)
        else:
        # Otherwise, 
            peptide_seq_records.append(SeqRecord(item))

    # After consolidating the sequences and records, QC sequence values and assign unique IDs 
    for peptide in peptide_seq_records:
        # QC the sequence, ensuring it has no invalid AA codes (X is a generic stand-in for NCAAs)
        if len(set(list(peptide.seq.__str__()))-set(valid_AA_codes))==0:
            # Ensure the sequence has a unique name; if not, create one:
            if peptide.id=="" or peptide.id=="<unknown id>":
                peptide.id=str(random.randint(1000000,9999999))
            if peptide.id in peptide_seq_names:
                while peptide.id in peptide_seq_names:
                    peptide.id=str(random.randint(1000000,9999999))
            peptide_seq_records_curated.append(peptide)
            peptide_seq_names.append(peptide.id)
        else:
            messenger.warning(f"Sequence {peptide.seq.__str__()} with invalid AA codes")
    
    return peptide_seq_records_curated

def add_ncaa(
        peptide,
        ncaa_code,
        toreplace,
        messenger
        ):
    '''
    Adds the NCAA to the input peptide sequence.
    Given a BioPython SeqRecord object, searches the sequence for
    either (1) a numeric index correponsing to AA position or (2) a
    one-letter AA code. Subsitutes a placeholder AA ("X") at that
    position(s) and returns the SeqRecord object with the updated
    sequence.

    Arguments:
    >peptide [obj]: Biopython SeqRecord object
    >ncaa_code [str]: 3-letter NCAA code
    >toreplace [str]: Position or AA code to replace with NCAA
    >messenger [obj]: Instance of Notice class

    Positional returns:
    [0] [obj]: Modified Biopython SeqRecord object (unchanged if no NCAA sub made)
    '''
    # First check: ensure that a change is requested
    if ncaa_code=="":
        return peptide
    
    # Make the substitution (if possible)
    pep_seq=peptide.seq.__str__()
    pep_seq_bak=pep_seq
    aa_or_idx="position"
    # The substitution can be either an AA or a position in the peptide
    # Use conditionals to test for both
    try:
        pep_seq[toreplace]="X"
    except:
        if toreplace in pep_seq:
            pep_seq=pep_seq.replace(toreplace,"X")
            aa_or_idx="amino acid"
        else:
            messenger.warning(
                    f"Unable to substitute NCAA {ncaa_code} at {aa_or_idx} {toreplace} "+
                    f"in peptide {peptide.id} with sequence {pep_seq}"
                    )
            return peptide

    peptide.seq=Seq(pep_seq)
    messenger.notice(
            f"Substituted NCAA {ncaa_code} at {aa_or_idx} {toreplace} in peptide {peptide.id}"
            )
    messenger.message(
            f"Old sequence: {pep_seq_bak}\nNew sequence: {pep_seq}"
            )
    return peptide

def params_file_validation(
        params_file,
        parent_dir,
        messenger
        ):
    
    '''
    Validates a Rosetta .params file and associated .pdb rotamers file
    exist. If they do, moves them to the parent_dir for prepacking and
    production run.

    Arguments:
    > params_file [str][path like]: Rosetta .params file
    > parent_dir [str][path like]: Directory to save .params and .pdb
    rotamers file to for use in production run
    > messenger [obj]: Instance of Notify class
    
    Positional returns:
    [0] [str][path like]: New .params file location
    [1] [str]: 3-letter NCAA code. Empty string denotes failure.
    [2] [bool]: Check for operation success. Useful for downstream tasks.
    '''
    # Case 1: Default of no params file passed by user in ahmp.py
    if params_file=="":
        return "","",True

    # Otherwise, make sure file and associated requirements are present.
    params_file=os.path.abspath(params_file)
    if Path(params_file).is_file():
        # Get the rotamer filename and validate it exists:
        params_dir=os.path.dirname(params_file)
        with open(params_file,"r") as pf:
            params_data=pf.readlines()
        rotamers_file=os.path.join(
                params_dir,
                params_data[-1].strip("\n").split(" ")[-1]
                )
        # Attempt to locate the rotamers file. If found, move it
        # to the parent directory and retrieve the 3-letter ncaa
        # code
        if Path(rotamers_file).is_file():
            shutil.copy(
                    params_file,
                    os.path.join(
                        parent_dir,
                        os.path.basename(params_file)
                        )
                    )
            shutil.copy(
                    rotamers_file,
                    os.path.join(
                        parent_dir,
                        os.path.basename(rotamers_file)
                        )
                    )
            # Get the 3 letter code for the NCAA
            ncaa_code=params_data[0].strip("\n").split(" ")[-1]
            new_params_loc=os.path.join(
                    parent_dir,
                    os.path.basename(params_file)
                    )
            messenger.notice(f"Successfully located params file {params_file} with "+
                    f"non-canonical AA {ncaa_code}")
            return new_params_loc, ncaa_code, True
        else:
            messenger.error(f"PDB rotamers file {rotamers_file} not located")
    else:
        messenger.error(f"Rosetta .params file {params_file} not located")

    return params_file,"",False


def make_flexpepdock_protocol(
        peptide,
        mhcdb,
        messenger,
        peptide_dir,
        params_file,
        ncaa_code,
        mhci_model,
        ignore_epitope_match,
        worst_epitope_template
        ):
    '''
    Make RosettaScripts .xml file that will create a prepacked pose for
    a given peptide bound to an HLA given a template backbone shape.

    Arguments:
    >peptide [obj]: Biopython SeqRes object with sequence of interest
    >mhcdb [obj]: Instance of MHCdatabase with assigned HLA allele
    >messenger [obj]: Instance of utilities.Notify class
    >peptide_dir [str][path like]: Directory containing input files
    >params_file [str][path like]: Path to Rosetta .params file, if using
    >ncaa_code [str]: 3-letter NCAA code
    >mhci_model [str][path like]: PDB to use as an MHC-I, overrides database search
    >ignore_epitope_match [bool]: If True, will ignore perfect template sequence matches
    >worst_epitope_template [bool]: If True, finds a structure template with the worst sequence match

    Positional returns:
    [0] [str][path like]: RosettaScripts .xml filename, or empty string if failed
    '''

    # Step 1: Find a template for our peptide
    messenger.message(f"Building {peptide.id}.xml...")
    # Remove the [ncaa_code] insert prior to sequence alignment
    # Set **kwargs values for get_peptide_template()
    query_sequence=peptide.seq.__str__()
    method="best"
    omit=[]
    if worst_epitope_template:
        method="worst"
    if ignore_epitope_match:
        omit=["self"]
    
    pdb_template=mhcdb.get_peptide_template(
            query_sequence,
            method=method,
            omit=omit
            )
    if not pdb_template:
        messenger.error(
                f"Failed to locate suitable template for peptide {peptide.id}"
                )
        return ""
    
    messenger.message(
            f"Using PDB template {pdb_template}"
            )

    # Copy the template PDB with a standardized file name (fpd_ppk_input.pdb) into
    # the peptide directory
    shutil.copy(
            os.path.join(
                mhcdb.location,
                "templates",
                pdb_template+".pdb"
                ),
            os.path.join(
                peptide_dir,
                "fpd_ppk_input.pdb"
                )
            )

    # Query the database to get information on chain assignments for the template structure
    template_data=mhcdb.query(
            action="get_templates",
            template=pdb_template
            )[0]
    mhc_chain=template_data["MHC_PDB_Chain1"]
    pep_chain=template_data["Antigen_PDB_Chain(s)"]
    template_sequence=template_data["Epitope_Description"]

    # Step 2: Make our RosettaScripts file
    root=ET.Element("ROSETTASCRIPTS")
    movers=ET.SubElement(root,"MOVERS")
    protocols=ET.SubElement(root,"PROTOCOLS")


    # SIMPLE THREADING MOVER
    ET.SubElement(
            movers,
            "SimpleThreadingMover",
            name="thread_peptide",
            thread_sequence=query_sequence,
            start_position=f"1{pep_chain}",
            pack_neighbors="1",
            skip_unknown_mutant="1",
            pack_rounds="1"
            )
    ET.SubElement(
            protocols,
            "Add",
            mover="thread_peptide"
            )

    # MUTATE RESIDUE MOVER
    # Add a mutate res mover and corresponding step in the protocol for each
    # ncaa in the peptide
    if "X" in query_sequence:
        for i,AA in enumerate(query_sequence):
            if AA=="X":
                ncaa_pos=i+1
                ET.SubElement(
                        movers,
                        "MutateResidue",
                        name=f"mutate_res_{ncaa_pos}",
                        target=f"{ncaa_pos}{pep_chain}",
                        new_res=ncaa_code
                        )
                ET.SubElement(
                        protocols,
                        "Add",
                        mover=f"mutate_res_{ncaa_pos}"
                        )
                messenger.message(
                        f"Added NCAA {ncaa_code} at position {ncaa_pos}"
                        )

    # DELETE OVERHANGING RESIDUES
    # A necessary step if the template peptide is longer than the query sequence
    if len(template_sequence)>len(query_sequence):
        ET.SubElement(
                movers,
                "DeleteRegionMover",
                name="delete_res_peptide",
                start=f"{len(query_sequence)+1}{pep_chain}",
                end=f"{len(template_sequence)}{pep_chain}"
                )
        ET.SubElement(
                protocols,
                "Add",
                mover="delete_res_peptide"
                )
        messenger.message(
                f"Deleting residual unmutated residues {len(query_sequence)+1}-{len(template_sequence)}"
                )
    # ADD MHC-I TO POSE
    # Can either be user-defined (mhci_model tag) or from the mhc database
    if not Path(mhci_model).is_file():
        # Checks for a file in the "models" directory of the mhc database matching
        # the allele name if no mhc-i model is provided by the user
        mhci_model=os.path.join(
                mhcdb.location,
                "models",
                mhcdb.allele+".pdb"
                )
        if not Path(mhci_model).is_file():
            messenger.error(f"MHC-I model {mhci_model} not located!")
            return ""
    
    # In order to use the AddChain mover, we must first ensure the MHC-I model we wish to add to
    # the pose has the same number of residues as the existing MHC-I from the template.
    parser=PDBParser(QUIET=True)
    template_struct=parser.get_structure(
            "template",
            os.path.join(
                mhcdb.location,
                "templates",
                pdb_template+".pdb"
                )
            )
    mhcimodel_struct=parser.get_structure(
            "mhci_model",
            mhci_model
            )
    # Quick check: There should only be 1 chain in our MHC-I model, and it should be labeled
    # 'A'. If otherwise, notify user of exceptional circumstances and adjust chain assignment
    # accordingly
    mhcimodel_chains=[c.id for c in mhcimodel_struct[0].get_chains()]
    mhcimodel_chainid=mhcimodel_chains[0]
    if len(mhcimodel_chains)>1:
        messenger.warning(
                f"MHC-I model {mhci_model} has {len(mhcimodel_chains)-1} unexpected chains!"
                )
    
    template_struct_length=len(
            [r.resname for r in template_struct[0][mhc_chain].get_residues()]
            )
    mhcimodel_struct_length=len(
            [r.resname for r in mhcimodel_struct[0][mhcimodel_chainid].get_residues()]
            )
    chain_len_diff=template_struct_length-mhcimodel_struct_length
    # If the mhc-i model we want to align to the pose has more residues than the existing
    # template structure, delete the overhanging residues to Rosetta doesn't crash
    if chain_len_diff>0:
        del_start=mhcimodel_struct_length+1
        del_end=template_struct_length
        ET.SubElement(
                movers,
                "DeleteRegionMover",
                name="delete_res_mhci",
                start=f"{del_start}{mhc_chain}",
                end=f"{del_end}{mhc_chain}"
                )
        ET.SubElement(
                protocols,
                "Add",
                mover="delete_res_mhci"
                )
        messenger.message(
                f"Removing {chain_len_diff} residues from PDB template {pdb_template} chain {mhc_chain}"
                )

    # Now we can add the chain of interest
    chain_to_swap=str(ord(mhc_chain)-64)
    ET.SubElement(
            movers,
            "AddChain",
            name="add_mhci_model_to_pose",
            file_name=mhci_model,
            swap_chain_number=chain_to_swap
            )
    ET.SubElement(
            protocols,
            "Add",
            mover="add_mhci_model_to_pose"
            )
    messenger.message(
            f"Added MHC-I model from file {mhci_model} to pose"
            )

    # RUN FLEXPEPDOCK PREPACK
    ET.SubElement(
            movers,
            "FlexPepDock",
            name="fpd_prepack",
            ppk_only="1"
            )
    ET.SubElement(
            protocols,
            "Add",
            mover="fpd_prepack"
            )
    
    # Make the .xml file
    xml_filename=os.path.join(
            peptide_dir,
            peptide.id+".xml"
            )
    xmltree=ET.ElementTree(root)
    ET.indent(xmltree, space="    ")
    xmltree.write(xml_filename)
    
    return xml_filename

def run_flexpepdock_protocol(
        xml_filename,
        peptide_dir,
        params_file,
        messenger,
        rosetta_main
        ):
    '''
    Subroutine to run the RosettaScripts .xml protocol generated
    by make_flexpepdock_protocol.

    Arguments:
    >xml_filename [str][path like]: RosettaScripts .xml file to run
    >peptide_dir [str][path like]: Directory to place output .pdb file
    >params_file [str][path like]: Location of params file, if using.
    >messenger [obj]: Instance of the Notify class
    >rosetta_main [str][path like]: Location of rosetta/main
    
    Positional returns:
    [0] [str][path like]
    '''
    # First, we check to make sure the rosetta executable exists:
    rosetta_scripts=os.path.join(
            rosetta_main,
            "source/bin/rosetta_scripts.linuxgccrelease"
            )
    if not Path(rosetta_scripts).is_file():
        messenger.error(
                f"Rosetta Scripts executable at {rosetta_scripts} not located!"
                )
        return ""
    
    # Set the arguments
    args=[
            rosetta_scripts,
            "-s", os.path.join(peptide_dir, "fpd_ppk_input.pdb"),
            "-parser:protocol", xml_filename,
            "-extra_res_fa", params_file,
            "-database", os.path.join(rosetta_main,"database"),
            "-out:path:all", peptide_dir,
            "-overwrite",
            "-mute", "all"
            ]

    messenger.notice(
            f"Running FlexPepDock prepack protocol"
            )
    process=subprocess.run(args)

    # Rename the output
    shutil.move(
            os.path.join(
                peptide_dir,
                "fpd_ppk_input_0001.pdb"
                ),
            os.path.join(
                peptide_dir,
                "input.pdb"
                )
            )

    return

def prep_production_run(
        peptide_dir,
        parent_dir,
        messenger,
        rosetta_main,
        params_file,
        rosetta_options,
        slurm,
        root_dir_prod,
        root_dir_input,
        prod_output
        ):
    '''
    Prepares Rosetta .options file and SLURM batch options
    file for FlexPepDocking production run on cluster.

    Arguments:
    > peptide_dir [str][path like]: Directory containing input files
    > parent_dir [str][path like]: Parent directory name
    > messenger [obj]: Instance of Notify class
    > rostta_main [str][path like]: Path to /rosetta/main/
    > params_file [str][path like]: Rosetta .params file (if using)
    > rosetta_options [list]: Rosetta options in flag/value pairings delimited by "="
    > slurm [list]: SLURM #SBATCH options in flag/value pairings delimited by "="
    > params_file [str][path like]: Path to Rosetta params file (if using)
    > root_dir_prod [str][path like]: Directory Rosetta is launched from during production
    > root_dir_input [str][path like]: Directory containing generated input files
    > prod_output [str]: Specify Rosetta output file type for production ("silent" or "score")

    Positional returns:
    '''
    # Define functions for parsing options into key:value pairs
    def parse_options(user_options, default_options):
        '''
        Description:
        Utility to add command line input to default rosetta/SLURM options

        Arguments:
        >user_options: options specified by user
        >default_options: default options dictated by this utility.

        Positional returns:
        [0] [dict]: Dictionary where key:value pairs correspond to command line argument and
        value(s)
        '''

        updated_options=default_options
        if user_options:
            for option in user_options:
                if "=" in option:
                    tmp_argval=option.split("=")
                    tmp_arg=tmp_argval[0]
                    tmp_val=tmp_argval[-1]
                else:
                    tmp_arg=option
                    tmp_val=""
                updated_options[tmp_arg]=tmp_val
        return updated_options

    # Some filepath manipulation to ensure Rosetta can find input files at production
    prod_parent_file_path=os.path.join(
            os.path.abspath(root_dir_input),
            os.path.basename(parent_dir)
            )
    prod_input_file_path=os.path.join(
            prod_parent_file_path,
            os.path.basename(peptide_dir)
            )

    # Define our default options for a FlexPepDock refinement run
    rosetta_options_default = {
            "database":os.path.join(rosetta_main,"database"),
            "nstruct":"250",
            "pep_refine":"",
            "ex1":"",
            "ex2aro":"",
            "in:file:s":os.path.join(
                prod_input_file_path,
                "input.pdb"
                )
            }

    # Add location of Rosetta .parms file (if using)
    if Path(params_file).is_file():
        rosetta_options_default[
                "extra_res_fa"
                ]=os.path.join(
                        prod_parent_file_path,
                        os.path.basename(params_file)
                        )

    # Specify our output type for Rosetta
    if prod_output=="score":
        rosetta_options_default[
                "out:file:scorefile"
                ]=os.path.basename(peptide_dir)+".sc"
    elif prod_output=="silent":
        rosetta_options_default[
                "out:file:silent"
                ]=os.path.basename(peptide_dir)+".silent"

    slurm_default = {
        "mail-user":"",
        "mail-type":"ALL",
        "ntasks":"1",
        "time":"24:00:00",
        "mem":"1G",
        "output":""
        }


    # Create our updated rosetta and SLURM batch options
    rosetta_options=parse_options(
            rosetta_options,
            rosetta_options_default
            )
    slurm=parse_options(
            slurm,
            slurm_default
            )

    if "out:path:all" not in rosetta_options.keys():
        default_outpath=os.path.join(
                root_dir_prod,
                os.path.basename(parent_dir),
                "output"
                )
        rosetta_options["out:path:all"]=default_outpath
        messenger.warning(
                f"Rosetta option -out:path:all not specified by user! "+
                f"Default directory {default_outpath} will be used. "+
                "Please ensure it exists before production run!"
                )

    if not slurm["mail-user"]:
        messenger.warning(
                f"#SBATCH --mail-user not set!"
                )
    
    if not slurm["output"]:
        messenger.warning(
                f"#SBATCH --output not set!"
                )
    
    rosetta_options_file=os.path.join(
            peptide_dir,
            "rosetta.options"
            )
    slurm_file=os.path.join(
            peptide_dir,
            os.path.basename(peptide_dir)+".slurm"
            )
    # Write our Rosetta options file:
    with open(rosetta_options_file,"w") as rof:
        for option,arg in rosetta_options.items():
            rof.write(f"-{option}")
            if arg:
                rof.write(f" {arg}")
            rof.write("\n")
    
    # Redefine our rosetta.options file path based on the root_dir_input argument:
    rosetta_options_file=os.path.join(
            prod_input_file_path,
            "rosetta.options"
            )

    # Write our .sbatch file:
    with open(slurm_file,"w") as sf:
        sf.write("#!/bin/batch\n")
        for option,arg in slurm.items():
            sf.write(f"#SBATCH --{option}")
            if arg:
                sf.write(f"={arg}")
            sf.write("\n")
        sf.write("module load GCC/6.4.0-2.28\n")
        flexpepdock=os.path.join(
                rosetta_main,
                "source/bin/FlexPepDocking.default.linuxgccrelease"
               )
        sf.write(f"{flexpepdock} @{rosetta_options_file}\n")

    return

def make_batch_file_list( 
        successful_peptide_inputs,
        messenger,
        parent_dir,
        root_dir_input
        ):
    '''
    Creates a file in root_dir_input/parent_dir that can be used in
    conjunction with the provided runbatch.sh shell script to
    sequentially send jobs to production on a server that uses
    SLURM.
    Assumes that SLURM batch files are named using the convention:
    <peptide_dir>/<peptide_name>.slurm

    Arguments:
    > successful_peptide_inputs [list]: List of BioPython SeqRecord objs
    > messenger [obj]: Instance of the Notify class
    > parent_dir [str][path like]: Path to directory containing peptide dirs
    > root_dir_input [str][path like]: Path to input files for production run
    '''

    # Define our script here:
    run_batch_shell_script='''
#!/bin/bash

# Takes in two positional inputs:
#   $1: .batch file containing list of directories to peptides
#   $2: [optional] SLURM job ID to set as a dependency
# -Iterates through each line in $1 and looks for a file named
# after the peptide with a ".slurm" extension in that directory
# -Runs the sbatch command with that file as input.
# -Saves the name of the .batch file to a record called
# "completedbatches.txt"

BATCHFILE=$1
while read line; do
    BASENAME="$(basename "${line}")"
    if [[ $2 -eq 0 ]]
    then
        sbatch "$line/$BASENAME.slurm"
    else
         sbatch "--dependency=afterany:$2" "$line/$BASENAME.slurm"
    fi
done < $1
printf "$BATCHFILE\\n" >> "completedvhts.batches"
'''
    batch_file_list_name=os.path.join(
            root_dir_input,
            os.path.basename(parent_dir),
            os.path.basename(parent_dir)+".batch"
            )
    with open(batch_file_list_name,"w") as f:
        for peptide_name in successful_peptide_inputs:
            f.write(os.path.join(
                    root_dir_input,
                    os.path.basename(parent_dir),
                    peptide_name
                    )
                )
            f.write("\n")

    run_batch_shell_script_name=os.path.join(
            parent_dir,
            "runbatch.sh"
            )
    if not Path(run_batch_shell_script_name).is_file():
        messenger.message(
                f"Creating {run_batch_shell_script_name} for production run"
                )
        with open(run_batch_shell_script_name,"w") as f:
            for line in run_batch_shell_script.split("\n"):
                f.write(line)
                f.write("\n")

    return

def prepare_input_files(
        peptide_list,
        mhcdb,
        parent_dir,
        messenger,
        rosetta_options,
        slurm,
        rosetta_main,
        params,
        mhci_model,
        ignore_epitope_match,
        worst_epitope_template,
        root_dir_prod,
        root_dir_input,
        prod_output
        ):
    '''
    Subroutine to handle preparation of input files for modeling_pipeline.
    Can be called as parallel, independent processes.

    Arguments:
    > peptide_list [list]: List of Bio.SeqRes objects, each with sequence and unique id
    > mhcdb [obj]: Instance of MHCdatabase class
    > parent_dir [str][path like]: Parent directory name
    > messenger [obj]: Instance of Notify class
    > rosetta_options [list]: Rosetta options in flag/value pairings delimited by "="
    > slurm [list]: SLURM #SBATCH options in flag/value pairings delimited by "="
    > rostta_main [str][path like]: Path to /rosetta/main/
    > params [list]: Rosetta .params file [0] and residue number(s) or AA code(s) to replace with NCAA [1:]
    > mhci_model [str][path like]: Path to custom model to use for modeling the HLA
    > ignore_epitope_match [bool]: Ignore perfect sequence matches when searching for templates
    > worst_epitope_template [bool]: Find the worst epitope template (instead of the best)
    > root_dir_prod [str][path like]: Directory Rosetta is launched from during production
    > root_dir_input [str][path like]: Directory containing generated input files
    > prod_output [str]: Specify Rosetta output file type for production ("silent" or "score")

    Positional returns:
    None
    '''
    # Iterate through the peptides and perform the following operations:
    #   1. Generate parent directory
    #   2. Create a new .fasta for the peptide of interest, modified with NCAA if applicable
    #   3. Create a RosettaScripts .xml file and used to create prepacked input file
    #   4a.Create a Rosetta .options file with user-specified options
    #   4b.Create a SLURM batchfile for the production run with user-specified options
    
    # Prework: validate the Rosetta .params file exists if passed by the user:
    params_file=params[0]
    params_validated=params_file_validation(
            params_file,
            parent_dir,
            messenger
            )
    params_file=params_validated[0]
    ncaa_code=params_validated[1]
    validation_success=params_validated[2]
    toreplace=params[1]

    if not validation_success:
        return
    # Iterate through the peptides and perform steps as described:
    successful_peptide_inputs=[]
    for peptide in peptide_list:
        messenger.notice(
                f"Peptide: {peptide.id:<15}Sequence: {peptide.seq}"
                )
        # Step 1: generate parent directory
        peptide_dir=os.path.join(
                parent_dir,
                peptide.id
                )
        Path(peptide_dir).mkdir(parents=True,exist_ok=True)
        
        # Step 2b: modify .seq record with NCAA if applicable
        peptide=add_ncaa(
                peptide,
                ncaa_code,
                toreplace,
                messenger
                )
        
        #Step 2c: generate new .fasta file, with NCAA if applicable
        # This is for record keeping purposes only. The peptide .fasta
        # file can be used for a FlexPepDock ab-initio run if so
        # desired.
        fasta_fn=os.path.join(
                peptide_dir,
                peptide.id+".fasta"
                )
        with open(fasta_fn,"w") as f:
            new_peptide_seq=peptide.seq.__str__().replace("X","["+ncaa_code+"]")
            f.write(f">{peptide.id}\n{new_peptide_seq}")

        # Step 3: create RosttaScripts .xml file
        xml_filename=make_flexpepdock_protocol(
                peptide,
                mhcdb,
                messenger,
                peptide_dir,
                params_file,
                ncaa_code,
                mhci_model,
                ignore_epitope_match,
                worst_epitope_template
                )
        if not Path(xml_filename).is_file():
            continue

        # Run the protocol to generate the input file
        run_flexpepdock_protocol(
                xml_filename,
                peptide_dir,
                params_file,
                messenger,
                rosetta_main
                )

        # Keep a record of peptide directories that we sucessfully generated
        # inputs for. We will consolidate these into a file that can be used
        # as input for a shell script that will submit them all to a SLURM
        # scheduler during a production run.
        successful_peptide_inputs.append(peptide_dir)

        # Prepare rosetta options and slurm files for production
        prep_production_run(
                peptide_dir,
                parent_dir,
                messenger,
                rosetta_main,
                params_file,
                rosetta_options,
                slurm,
                root_dir_prod,
                root_dir_input,
                prod_output
                )

    
    # Now prepare a list of all sbatch files to run from root_dir_input
    messenger.notice(
            f"Created input files for {len(successful_peptide_inputs)} peptides"
            )
    return
      
def modeling_pipeline(
        peptides,
        allele,
        rosetta_options,
        slurm,
        batch_size,
        mhcdb_location,
        rosetta_main,
        params,
        mhci_model,
        verbosity,
        ignore_epitope_match,
        worst_epitope_template,
        root_dir_prod,
        root_dir_input,
        prod_output,
        multitask
        ):
    '''
    Pipeline for threading peptides onto PDB models to serve
    as starting structures for FlexPepDock refinement.

    For the given HLA, iterates through each peptide and produces a folder
    named after the peptide in the provided (or generated) .fasta file.
    Generates the following filetree structure:

    ./HLA_name/
        peptide_1/
            file_1
            file_2
            ...
            file_n
        ...
        peptide_n/

    Arguments:
    > peptides [list]: Peptide sequence(s) or file(s) containing peptide sequences
    > allele [str]: HLA allele in the format <GENE>*<ALLELE_GROUP>:<PROTEIN>
    > rosetta_options [list]: Rosetta options in flag/value pairings delimited by "="
    > slurm [list]: SLURM #SBATCH options in flag/value pairings delimited by "="
    > batch_size [int]: Number of simultaneous production runs
    > mhcdb_location [str][path like]: Location of the MHC database
    > rostta_main [str][path like]: Path to /rosetta/main/
    > params [list]: Rosetta .params file [0] and residue number to replace with NCAA [1]
    > mhci_model [str][path like]: Path to custom model to use for modeling the HLA
    > verbosity [str]: Toggle output logged to console ("all","warnings","errors",or "none")
    > ignore_epitope_match [bool]: Ignore perfect sequence matches when searching for templates
    > worst_epitope_template [bool]: Find the worst epitope template (instead of the best)
    > root_dir_prod [str][path like]: Directory Rosetta is launched from during production
    > root_dir_input [str][path like]: Directory containing generated input files
    > prod_output [str]: Specify Rosetta output file type for production ("silent" or "score")
    > multitask [bool]: If True, run thread_model in parallel using available resources

    Positional returns:
    None
    '''
    if multitask:
        verbosity="notifications"
    
    # Set up our MHC database object and validate our allele:
    mhcdb=MHCdatabase(
            location=mhcdb_location,
            toprint=verbosity,
            build=True
            )
    check_success=mhcdb.set_allele(allele)
    if not check_success:
        return
    
    # Set up our notifications object:
    messenger=Notify(toprint=verbosity)

    # Obtain and QC our list of peptides:
    peptide_list=create_peptide_list(peptides, messenger)
    if len(peptide_list)==0:
        messenger.error(f"No valid peptide sequences found")
        return
    
    # Create our parent directory for the allele of interest
    # While the directory will be made in the current working directory,
    # the user can specify a different parent directory for the production
    # run using the --root_dir_input flag
    dirname=mhcdb.allele.replace("*","").replace(":","")
    parent_dir=os.path.join(
            os.getcwd(),
            dirname
            )
    Path(parent_dir).mkdir(parents=True,exist_ok=True)
    
    # Thread peptides and create input files (with or without multitasking):
    if multitask:
        # Default CPUs to use is all available minus 1
        available_cpus=multiprocessing.cpu_count()-1
        batch_size=math.ceil(len(peptide_list)/available_cpus)
        peptide_batchlist=[
                peptide_list[i:i+batch_size] for i in range(0,len(peptide_list),batch_size)
                ]
        # For each list of peptides in peptide_batchlist, prepare input files
        pool=multiprocessing.Pool(available_cpus)
        for peptide_batch in peptide_batchlist:
            pool.apply_async(
                   prepare_input_files,(
                       peptide_batch,
                       mhcdb,
                       parent_dir,
                       messenger,
                       rosetta_options,
                       slurm,
                       rosetta_main,
                       params,
                       mhci_model,
                       ignore_epitope_match,
                       worst_epitope_template,
                       root_dir_prod,
                       root_dir_input,
                       prod_output
                       )
                    )
        pool.close()
        pool.join()
    else:
    # Otherwise, prepare input files as normal for each peptide in sequence.
        prepare_input_files(
                peptide_list,
                mhcdb,
                parent_dir,
                messenger,
                rosetta_options,
                slurm,
                rosetta_main,
                params,
                mhci_model,
                ignore_epitope_match,
                worst_epitope_template,
                root_dir_prod,
                root_dir_input,
                prod_output
                )
    # Make our batch file list:
    successful_peptide_inputs=[
            x for x in os.listdir(parent_dir) if Path(x).is_dir()
            ]
    make_batch_file_list(
        successful_peptide_inputs,
        messenger,
        parent_dir,
        root_dir_input
        )

    return
