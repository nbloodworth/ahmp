'''
ahmp.modules.modeling

ahmp module for generating FlexPepDocking input.
'''
# Python standard libraries

# ahmp
from utilities import Notify

def build_fasta():
    return

def thread_template(
        input_seq, 
        mhc_db, 
        rosetta_main,
        params="",
        ):
    '''
    Description:
    Builds a RosettaScripts protocol for producing input to
    FlexPepDocking refinement. Optionally adds an NCAA at
    user-determined postion.
    
    Arguments:
    >input_seq [obj]: BioPython sequence record 

    Outputs:
    
    '''
    allele=db.allele
    if user_defined_receptor!=None and Path(user_defined_receptor).is_file():
        chk_success=db.set_allele(allele, build_receptor=False) 
    else:
        chk_success=db.set_allele(allele)# This will automatically build the receptor if it is not present and create a list of alleles in the structure database to use as peptide backbone templates. Building a receptor with AlphaFold2 could take upwards of ~1h
    if not chk_success:
        notify().error(f"Failed to set allele to {allele}")
        return False

 # Retrieve the peptide sequences in the structural database for either the allele or (if unavailable) its next nearest sequence neighbor struct_data_fn=os.path.join(db.locate,"database.info")
    if Path(struct_data_fn).is_file(): 
        struct_data=pd.read_csv(struct_data_fn)
    else:
        notify().error(f"Unable to locate {struct_data_fn}! Does it exist?")
        return False

     # Retrieve the default receptor location
     receptor_default=db.receptor

     # Import the fasta file for the sequence to thread
     input_record = SeqIO.read(input_fasta, "fasta")
     input_record_originalSeq=input_record.seq # Save the original sequence as we will be modifying the record sequence later.
     
     # Clean the sequence for alignment if it is formatted for Rosetta with a NCAA placeholder, and save the 3-letter code for the NCAA for later use
     if "X" in input_record.seq:
         if not params:
             notify().error("NCAA placeholder character 'X' found in {}, but no params file provided\033[0m".format(input_fasta))
             return False
         else:
             pfile=open(os.path.abspath(params),"r")
             NCAA_name=pfile.readline().strip()[-3::] #Assumes that the three letter code for the NCAA is identical to what is in the provided params file (it should be, but this might miss edge cases where it is not, which would also cause Rosetta to crash later in the pipeline)
             pfile.close()
             input_record.seq=Seq(str(input_record.seq).replace("["+NCAA_name+"]",""))

     # Retrieve the peptide backbone template:
     if ignore_perfect_match:
         to_omit=["self"]
     else:
         to_omit=[]
     if find_worst_template:
         template_id_method="worst"
     else:
         template_id_method="best"

     # Get the PDB to use as a template for peptide threading (based on sequence similarity)
     if "X" in input_record.seq:
         queryseq=str(input_record.seq)
     else:
         queryseq=str(input_record_originalSeq)
     
     best_pdb=db.get_peptide_template(queryseq, method=template_id_method, omit=to_omit)
     # Retrieve the sequence for the best template identified
     best_sequence_data=db.query_mhc_database(gettemplate=best_pdb, quiet=True)[0]
     best_sequence=best_sequence_data["Epitope_Description"]

     best_pdb_fn=os.path.abspath(os.path.join(db.locate,"templates",best_pdb+".pdb"))
     # Now we build an rosetta_scripts.xml file to do the following:
     #     (1) Sequentially mutate residues using the simple threading mover to match our desired sequence
     #     (2) If the template sequence is longer than our input sequence by N residues, delete the last N residues from the threaded model
     #     (3) Add the template (or user defined) MHC receptor to the threaded peptide from the database of aligned structures
     #    (4) Call FlexPepDock prepack on the output structure in preparation for refinement run 

     parent_dir=os.path.abspath(os.path.dirname(input_fasta))
     XML_filename = "{}/flexpepdock_input_{}.xml".format(parent_dir, os.path.basename(input_fasta)[0:-6])
     args=["{}/source/bin/rosetta_scripts.linuxgccrelease".format(ROSETTA_PATH), "-s", best_pdb_fn, "-parser:protocol", XML_filename, "-database", "{}/database".format(ROSETTA_PATH), "-overwrite","-mute","all"]
     N=len(best_sequence)-len(input_record.seq) #for resolving length disparities between input sequence and selected template

     # Get the correct peptide chain from the file
     db_records=pd.read_csv(os.path.join(db.locate,"database.info"))
     pChain=db_records["Antigen_PDB_Chain(s)"].loc[db_records["PDB_ID"]==best_pdb].values[0]
     mhcChain=db_records["MHC_PDB_Chain1"].loc[db_records["PDB_ID"]==best_pdb].values[0]

     # STEP 1: SIMPLE THREADING MOVER
     print("Building .xml for input model generation...")
     root=ET.Element("ROSETTASCRIPTS")
     movers=ET.SubElement(root,"MOVERS")

     movers_threadModel=ET.SubElement(movers, "SimpleThreadingMover", name="threadModel", pack_neighbors="1", start_position="1{}".format(pChain), thread_sequence=input_record.seq, skip_unknown_mutant="1", pack_rounds="1")
     protocols=ET.SubElement(root,"PROTOCOLS")
     protocols_add_threadModel=ET.SubElement(protocols,"Add", mover="threadModel")

     #    Step 1b: Mutate NCAA if present
     if params:
         for i,AA in enumerate(input_record.seq):
             if AA=="X":
                 tmp_pos=i+1
                 protocols_add_mutateRes=ET.SubElement(protocols,"Add", mover="mutateRes_{}".format(tmp_pos))
                 movers_mutateRes=ET.SubElement(movers, "MutateResidue", name="mutateRes_{}".format(tmp_pos), target="{}{}".format(tmp_pos,pChain), new_res=NCAA_name)
         args.extend(["-extra_res_fa", os.path.abspath(params)])

     # STEP 2: DELETE OVERHANGING RESIDUES
     if N>0:
         protocols_add_deleteRes=ET.SubElement(protocols,"Add", mover="deleteRes")
         movers_deleteRes=ET.SubElement(movers, "DeleteRegionMover", name="deleteRes", start="{}{}".format(len(input_record.seq)+1,pChain), end="{}{}".format(len(best_sequence),pChain))

     # STEP 3: ADD TEMPLATE (OR USER DEFINED) MHC TO POSE

     # Determine if we should use a user-defined receptor or the default receptor for this allele, and make sure the file(s) exist
     if user_defined_receptor!=None:
         if Path(os.path.abspath(user_defined_receptor)).is_file():
             print(f"Adding user defined receptor: {os.path.basename(user_defined_receptor)}")
             replace_receptor_fn=os.path.abspath(user_defined_receptor)
         else:
             notify().warn(f"User defined MHC model {os.path.abspath(user_defined_receptor)} not located.")
             print(f"Adding default receptor for {db.allele}")
             replace_receptor_fn=receptor_default
     else:
         print(f"Adding default receptor {db.receptor} for {db.allele}")
         replace_receptor_fn=receptor_default

     # Read in the file and copy it to a temporary file - IMPORTANT - assumes that the default receptor will have a single chain named "A"
     pdb_parser=PDBParser(QUIET=True)
     new_receptor_fn = os.path.join(parent_dir,"tmp_receptor.pdb") #Create temporary receptor file to graft into pose
     new_receptor = pdb_parser.get_structure("defaultReceptor", replace_receptor_fn)
     io=PDBIO()
     
     for chain in new_receptor.get_chains():
         if chain.get_id()=="A":
             io.set_structure(chain)
             io.save(new_receptor_fn)

     print("Temporary file {} created for adding to pose".format(new_receptor_fn))
     # First make sure the template and template receptor have the same number of residues.
     structure=pdb_parser.get_structure("template",best_pdb_fn)
     model=structure[0]
     chain=model[mhcChain]
     template_rescount=0
     for residue in chain:
         template_rescount+=1

     structure=pdb_parser.get_structure("receptor",new_receptor_fn)
     model=structure[0]
     chain=model["A"]
     receptor_rescount=0
     for residue in chain:
         receptor_rescount+=1     
     # Use delete region to remove excess residues from the template
     resdiff=template_rescount-receptor_rescount
     if resdiff>0:
         print(f"Deleting residual residues {receptor_rescount}-{template_rescount} from {best_pdb}")
         protocols_add_deleteRes2=ET.SubElement(protocols,"Add", mover="deleteRes_receptor")
         movers_deleteRes2=ET.SubElement(movers, "DeleteRegionMover", name="deleteRes_receptor", start="{}{}".format(receptor_rescount+1,mhcChain), end="{}{}".format(template_rescount,mhcChain),rechain="1")

     # Use addchain mover to add new receptor to pose
     # Determine chain to swap
     chaindict={"A":1,"B":2,"C":3,"D":4,"E":5,"F":6,"G":7,"H":8,"I":9,"J":10,"K":11,"L":12,"M":13,"N":14,"O":15,"P":16,"Q":17,"R":18,"S":19,"T":20,"U":21,"V":22,"W":23,"X":24,"Y":25,"Z":26}
     chain_toSwap=chaindict[db_records["MHC_PDB_Chain1"].loc[db_records["PDB_ID"]==best_pdb].values[0]]

     # Add the new receptor. This will also (conveniently) align the new receptor to the old receptor.
     protocols_add_addReceptor=ET.SubElement(protocols,"Add", mover="addReceptor")
     movers_addReceptor=ET.SubElement(movers, "AddChain", name="addReceptor", file_name="{}".format(new_receptor_fn), swap_chain_number=f"{chain_toSwap}")

     # STEP 4: RUN FLEXPEPDOCK PREPACK
     protocols_add_prepack=ET.SubElement(protocols,"Add", mover="prepack")
     movers_prepack=ET.SubElement(movers, "FlexPepDock", name="prepack", ppk_only="1")

     # Write the final xml file
     xmltree=ET.ElementTree(root)
     ET.indent(xmltree,space='    ')
     xmltree.write(XML_filename)

     # Run the script in Rosetta
     notify().notice("Now running Rosetta FlexPepDock prepack protocol...")
     p=subprocess.run(args)

     # Rename the output .pdb file
     output_fn_old = os.path.join(parent_dir, "{}_0001.pdb".format(best_pdb))
     output_fn_new = os.path.join(parent_dir, "{}_input.pdb".format(os.path.basename(input_fasta)[0:-6]))

     shutil.move(output_fn_old, output_fn_new)

     print("\n{} saved and ready for production run.".format(output_fn_new))

     return output_fn_new # return the output filename for use in thread_all

    return

def thread_all():
    return
