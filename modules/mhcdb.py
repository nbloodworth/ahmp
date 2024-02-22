'''
 MHCdb is a class with a collection of methods and attributes used
 create and access information contained within a database of all HLA
 allele sequences and known structures. The class is designed to facilitate
 computational modeling of peptide-MHC structures with downstream
 applications (such as Rosetta or py.Rosetta)

 ============================PYTHON DEPENDENCIES============================
 BioPython
 Pandas
 Requests

==================================MODULES===================================

 ==========================DATABASE FILE STRUCTURE==========================
 *All .info files are .csv files organized by column:value format and can
 be read into a pandas dataframe

 Database folder/file------->Description
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 >MHCdb/------------->Parent folder
 |-->build/----------------->Folder containing files important for database
    |                        build
    |-->allele_list.info---->File containing full sequences for each allele,
    |                        the allele the sequence belongs to, and a list
    |                        of alleles with at least one PDB to use as
    |                        templates (in order of sequence identity). The
    |                        templates list is time consuming to create and
    |                        is made "as needed", and afterwards the
    |                        templates list is saved to this file so it
    |                        won't need to be created again if necessary.
    |-->allele_seq.info ---->File with all unique alpha1/alpha2 domain
    |                        sequences for every HLA allele, and a list of
    |                        HLA alleles with that sequence. Used for
    |                        building receptor models with AlphaFold2.
    |-->mhc_3d_assays.csv--->CSV file from IEDB with list of PDBs and PDB
    |                        info for all MHC models
    |-->tcr_3d_assays.csv--->CSV file from IEDB with list of PDBs and PDB
    |                        info for all TCR/MHC models
 |-->templates/------------->Folder containing all PDB templates for all
 |                           MHC alleles
 |-->models/---------------->Folder containing MHC-I PDBs to use as models
 |                           for peptide docking
 |-->allele_struct.info----->File containing a summary of all MHC/peptide
                             PDB structures in database.
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ==========================CLASSES AND SUBROUTINES==========================

 [class] chainSelect
 Overloads default Select class in Biopython's PDBIO module to allow
 selection of chains with MHC-I and peptide antigen. Also removes HETATMs.
 Called by self.build() method

 [class] MHCdatabase
    Creates an instance of the MHC database. From this instance further
    subroutines can be called, and properties inferred.
 ---------------------------------------------------------------------------
    [attribute] locate
       Location of database; can be set with constructor.
    [attribute] [hidden] __LOCATE__
       Default location of database. Defaults to
       /home/<user>/Data/MHCdb.
    [attribute] allele
       The allele that the MHC database object is currently pointed towards.
    [attribute] template
       The allele, present in the database, that will be used as a template
       for peptide backbone threading.
  ---------------------------------------------------------------------------
    [method] __init__
      Constructor for a database object. Optionally takes in database
      location if different from default (/home/usr/Data/MHCdb/), an
      option to set the allele and/or template, and an option to build the
      database if it does not exist.
  ---------------------------------------------------------------------------
     [method] build(allele_list=None)
       Builds the database at location specified by self.location. If
       self.location does not exist, will make the folder and parents.
  ---------------------------------------------------------------------------
     [method] (check_allele, get_templates=True,
                 build_receptor=True)
        Sets an allele for the database instance. Prints a warning and
        returns False if allele requested does not exist. Also sets allele
        structural templates if no structure for the allele is available,
        and builds a receptor model with AlphaFold2 if one is not already
        created (optional).
  ---------------------------------------------------------------------------
     [method] query(allele="all", quiet=False, getpdb=False,
                 getepitopes=False, gettemplate=False)
        Method for retrieving information from the database. Variable
        behavior depending on which option is set by the user, but always
        returns a tuple with the first element corresponding to a data
        structure (list or dict) and the second a boolean. Outputs for each
        option are below in priority order:

        allele="all": Returns a list of all available alleles (default
        behavior). Returns a True value if the value for
        MHCdatabase().allele is present in 

        allele="<some HLA here>": Returns a list of all available alleles
        and a boolean value indicating whether or not requested allele is
        present (True) or absent (False) (ignores the value for
        MHCdatabase().allele for the purpose of this output)

        getpdb=True: Prints and returns a list of PDBs with templates for
        MHCdatabase().allele. If no allele value set for database instance,
        will return all PDBs belonging to the variable allele. If that
        variable is set to the default ("all"), then will return all PDBs
        in database. If no PDBs are found, also returns False.

        getepitopes=True: Returns a dictionary where each key corresponds
        to a PDB ID in the database belonging to the allele assigned to
        the database instance (MHCdatabase().allele), and the value is the
        epitope sequence. If MHCdatabase().allele is set to None, returns
        all PDB:epitope sequences in the database. If no PDBs for the
        allele in quesiton are found, returns False.

        gettemplate=<4-letter-PDB-code>: Returns a dictionary where each
        key corresponds to  column values for a given PDB in
        the database, and a boolean indicating whether or not the PDB was
        found. Keys include:
        'Epitope_IRI', 'Epitope_Description', 'MHC_Allele',
        'MHC_PDB_Chain1', 'Antigen_PDB_Chain(s)', 'PDB_ID',
        'Resolution_(Angstrom)', 'Epitope_ID'
        If the PDB is not found in the database, returns an empty dict and
        False.

        quiet=True: Supress the printing of execution info to screen by
        setting quiet=True
  ---------------------------------------------------------------------------
     [method] __fetch_mhc_allele_info(allele_list=None, update=False)
        Fetches a list of all MHC alleles and sequences available in either
        a provided fasta file (optional argument allele_list=<filename or
        URL>) or the IMGT database. Writes the allele_list.info and
        allele_seq.info files in MHCdb/build/.
  ---------------------------------------------------------------------------
     [method] __set_mhc_templates()
        Method for setting the allele templates, if self.allele is not
        present in the MHC structure database. Otherwise,
        self.template=self.allele.
  ---------------------------------------------------------------------------
     [method] __trim_seq(seq, HLA_allele)
        Method that takes in an HLA sequence and allele and returns the
        alpha 1 and alpha 2 domains only to build a receptor model.
  ---------------------------------------------------------------------------
     [method] build_receptor(overwrite=False,
                 AF2_MINICONDA="/sb/apps/alphafold211/miniconda3",
                 AF2_DATADIR="/sb/apps/alphafold-data",
                 AF2_REPO="/sb/apps/alphafold211/alphafold")
        Wrapper for a call to AlphaFold2 to build the MHC receptor based on
        the trimmed alpha 1 and alpha 2 domain(s). AF2_MINICONDA refers to
        the location of the local installation of the AlphaFold2 conda
        environment; AF2_DATADIR refers to the local AlphaFold2 database,
        and AF2_REPO refers to the local AlphaFold2 repository.
 ---------------------------------------------------------------------------
     [method] get_peptide_template(self, query_sequence, method="best",
                 cutoff=0, omit=["self"]):
        Method to retrieve a structure from the database whose peptide
        backbone sequence matchest that of query_sequence. Returns the PDB
        code that corresponds to the matching structure in
        <database>/templates. If execution fails, returns False.

        Description of method variables:
        query_sequence: the sequence we are looking to find a template for
        method:   The method to perform the match (string). Options include:
          "best":      [default] Finds the best matching template by
                       alignment similarity score using a BLOSUM62 matrix for
                       AA evolutionary similarity
          "worst":     As above, but returns the template with the worst
                       score
          "bestabove": Finds the best matching template by %% identity above
                       a threshold given by the variable 'cutoff'. Does not
                       consider evolutionary similarity between AA
          "worstbelow":Finds the worst matching template by %% identity
                       below a threshold given by the variable 'cutoff'
        cutoff:   Set a value between 0 and 1 that indicates a threshold for
                  sequence identity for "bestabove" and "worstbelow"
                  methods.
        omit:     A list of PDB templates to omit from template selection
                  (4-character PDB codes). If "self" is included in the
                  list, any template with a 100% matching sequence will be
                  ignored.
 ---------------------------------------------------------------------------
     [method] help()
        Prints this help message and exits.

 =============================USE AND EXAMPLES==============================
  Use examples:

 #Example 0: Print these usage instructions:

    from HLA_db import MHCdatabase
    MHCdatabase().help()

    >MHCdb is a class with a collection of methods and attributes
    >used create a and access information contained within a database of
    >all HLA allele sequences and known structures. The class is designed
    ...

 #Example 1: Create an instance of the database: (note: will NOT build
             database or assign an allele to the instance unless specified
             by user)

     from HLA_db import MHCdatabase
     db=MHCdatabase()

 #Example 2: Build the database:

     db.locate="/home/user/database-location" #Optionaly specify location
    db.build()

 #Example 3: List all alleles with structural templates in database:

    allele_list=db.query(quiet=True)
    print(allele_list[0])

    >['A*01:01', 'A*02:01', 'A*02:03', 'A*02:06'...]

 #Example 4: Assign an allele to the database instance and initiate query:

    print(db.allele) #Before assigning an allele, the default value is None

    >None

    db.("A*02:01", get_templates=False, build_receptor=False)
    print(db.allele)

    >A*02:01 #Value is a string with format <gene>*<allele group>:<protein>

 #Example 5: Retrieve all PDBs in the database matching this allele:

    pdb_list=db.query(quiet=True, getpdb=True)
    print(pdb_list[0])

    >['3MRE', '3MRG', '3D25', '7KGQ', ...]

 #Example 6: Retrieve peptide template from database based on sequence
             similarity closest match.

     query_sequence="ALAKIMKAP"
     pdb_template=db.get_peptide_template()
     print(pdb_template)

     >6R2L
'''

# ===============================================================================================
# To do:
# [X]    Add a step to self.build() that renumbers the receptor chain and peptide chain when building the database. There are a lot of messy PDBs that confuse Rosetta later in the pipeline. Fixed them manually for now, but this will help with most of them (they throw an "assertion "). The main problem is that the residues listed in  are not all present in the actual PDB.
#  [ ] Add additional functionality to remove_bad_pdbs that considers residues with single atoms at the N-terminus (a surprising number of failure cases are caused by this). Can try running simple Rosetta thread peptide implementation, or double checking that each residue has the expected number of atoms
# [ ] Add functionality to include MHC-I from other species (mouse, primarily)
# [ ] Add a cluster_templates method that assess template similarity by backbone structure rather than sequence similarity, and benchmark.
# [ ] Implement a "remote query" function, that can retrieve templates directly from the PDB rather than storing them locally.


## IMPORTS ##
# Python standard library
import os
import sys
from pathlib import Path
from zipfile import ZipFile
import subprocess
import warnings

# Biopython
from Bio import BiopythonWarning
from Bio import BiopythonDeprecationWarning
warnings.simplefilter('ignore',BiopythonWarning)
warnings.simplefilter('ignore',BiopythonDeprecationWarning)

from Bio import SeqIO
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select
from Bio import Align
from Bio.Align import substitution_matrices

# requests
import requests

# pandas
import pandas as pd
warnings.simplefilter(action='ignore', category=FutureWarning)

# ahmp modules
from modules.utilities import Notify

## BEGIN FUNCTION DEFINITIONS##
def remove_bad_pdbs(templates_dir, struct_data, notifications):
    '''
    Description:
    A method to quickly error check all of our templates
    and get rid of bad ones. Based on the presumption 
    that most errors are due to low peptide electron 
    density and resultant missing residues in the PDB.

    Arguments:
    >templates_dir: [str] Directory to the downloaded PDB templates
    >struct_data: [pd.DataFrame obj] data from the  file
    >notificatoins: [Notify obj] Notify object for controlling output urgency
    '''
    bad_pdbs=[]
    fix_pdbs=[]
    bad_chain=[]
    bad_res=[]
    reported_seq=[]
    pdb_seq=[]
    pdbtot=len(struct_data)
    pdbnum=0
    for pdb in struct_data["PDB_ID"].tolist():
        # First, ensure the sequence of the PDB template matches what is recorded in the database
        pdbnum+=1
        notifications.message(f"Now error checking PDB {pdbnum:<5} of {pdbtot} in {templates_dir}...",endline="\r")
        tmp_reported_seq=struct_data.loc[
                struct_data["PDB_ID"]==pdb,["Epitope_Description"]
                ].iloc[0,0]
        epitope_chain_id=struct_data.loc[
                struct_data["PDB_ID"]==pdb,["Antigen_PDB_Chain(s)"]
                ].iloc[0,0]

        for record in SeqIO.parse(os.path.join(templates_dir,pdb)+".pdb","pdb-atom"):
            # print(record.id)
            if record.annotations["chain"]==epitope_chain_id:
                if tmp_reported_seq != str(record.seq):
                    bad_pdbs.append(pdb)
                    reported_seq.append(tmp_reported_seq)
                    pdb_seq.append(str(record.seq))

        # If the sequence is OK, now we iterate through the residues and make sure
        # the expected number of heavy atoms are present in each residue
        tmp_struct=PDBParser(
                PERMISSIVE=1,
                QUIET=1
                ).get_structure(
                        pdb,
                        os.path.join(templates_dir,pdb+".pdb")
                        )[0]
        rescount=1
        for chain in tmp_struct.get_chains():
            for residue in chain.get_residues():
                atomcount=0
                for atom in residue.get_atoms():
                    atomcount+=1
        # If the number of atoms is unexpectedly low for a given residue, 
        # remove it from the template (a more or less arbitrary fix which 
        # takes care of edge cases where the residue at the N-terminus of 
        # the MHC-I receptor has no atoms). Check only the terminus atoms
                if atomcount<3:
                    bad_chain.append(chain)
                    bad_res.append(rescount)
                    fix_pdbs.append(pdb)
                rescount+=1
    notifications.message("",endline="\n")

    # Remove the bad terminal residues from PDBs identified earlier
    if len(fix_pdbs)>0:
        notifications.message("Removing bad residues from: {}".format(", ".join(fix_pdbs)))
        for pdb,chain,res in list(zip(fix_pdbs,bad_chain,bad_res)):
            tmp_fn=os.path.abspath(os.path.join(templates_dir,pdb+".pdb"))
            tmp_struct=PDBParser(PERMISSIVE=1, QUIET=1).get_structure(pdb,tmp_fn)[0]
            io=PDBIO()
            io.set_structure(tmp_struct)
            io.save(tmp_fn,residueSelect(chain,res))

    # Print information on PDBs that were culled because their epitope sequence did not match the database sequence
    notifications.notice(f"Error checking complete. {len(bad_pdbs)} bad pdbs found and removed)")
    
    # Return the culled data 
    return struct_data

## BEGIN CLASS DEFINITIONS##
class residueSelect(Select):
    '''
    Description:
    Subclass for overloading default selector class in Bio.PDB.PDBIO.
    Defines "valid" residues for selection by testing if they
    belong to the desired chain(s).
    '''
    def __init__(self,rem_chain,rem_res):
        self.rem_chain=rem_chain
        self.rem_res=rem_res

    def accept_residue(self,residue):
        if residue.get_id()[1]==self.rem_res:
            return False
        else:
            return True

class chainSelect(Select):
    '''
    Class to help us select chains from the PDBs (both MHC and peptide chains)
    Overloads the default Select class in Bio.PDB.PDBIO. Also includes
    a method to omit HETATMs from the template structure.
    '''
    def __init__(self,MHC,pep):
        self.MHC_chain=MHC
        self.pep_chain=pep

    def accept_chain(self,chain):
        if chain.get_id()==self.MHC_chain or chain.get_id()==self.pep_chain:
            return True
        else:
            return False

    # Omit HETATMs
    def accept_atom(self, atom):
        # In BioPDB terms, a space at the 3,0 position of get_full_id 
        # corresponds to a HETATM
        if atom.get_full_id()[3][0]==' ':
            return True
        else:
            return False

## MHCdatabase CLASS ##
class MHCdatabase:
    '''
    Description:
    Class defining the attributes and methods of an MHCdatabase object.

    Attributes:
    >allele: Current assigned allele
    >templates: Assigned allele templates
    >sequence: Assigned allele alpha1/alpha2 domain sequence
    >model: Assigned allele PDB to use in downstream modeling
    >location: Root directory of /MHCdb
    >IMGTHLA: Location of the IMGTHLA database to download
    >IEDB: Location of the IEDB database to download

    Methods:

    '''
    ## ATTRIBUTES ##
    allele=None
    templates=[]
    sequence=""
    model=None
    location=""
    species="human" # Variable that can be set in future iterations
                    # to specify species-specific behavior.
    IMGTHLA="https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/hla_prot.fasta"
    IEDB="https://www.iedb.org/downloader.php?file_name=doc/iedb_3d_full.zip"         

    ## METHODS ##
    def __init__(self, location, 
            build=False,
            allele="",
            toprint="all",
            species="human"
            ):
        '''
        Description:
        Constructor for the MHCdatabase class.

        Arguments:
        >location: [str] Location of the MHC database
        >build: [bool][optional] Indication that database should be built
        >allele: [str][optional] Assign an allele for query purposes
        >toprint: [bool][optional] Message level (see utilities.Notify module)
        >species: [str][optional] Specify the allele species 
        '''
        # Set attribute values
        self.notifications=Notify(toprint=toprint)
        self.species=species
        self.location=location
        # Set the database location and ensure it is valid
        locate_success=self.locate()
        database_built=False
        # If no database found at self.location, then exit
        if not locate_success[0]:
            sys.exit(1)
        # Otherwise, if a location was not found but was
        # created and user requests build, execute self.build() 
        elif locate_success[1] and build:
            self.notifications.notice(
                f"Now building MHC database at {self.location}"
                )
            self.build()
            database_built=True
        # Prevent self.build() from executing if called from
        # the constructor IF the database already exists
        # (force user to override existing data)
        elif not locate_success[1] and build:
            self.notifications.warning(
                f"Call to build database at {self.location}, but database already exists!\n"+
                "Call MHCdatabase.build() to force an update."
                )

        # Only assign an allele if database exists or was sucessfully built
        if (locate_success[0] and allele) or (database_built and allele): 
            # If an allele is set, then assign sequence and template data, and
            # ensure we  have a valid model.
            self.set_allele(allele)
    
    def locate(self):
        '''
        Description:
        The purpose of this method is to return a valid location for the
        MHC database. Returns False if it fails to find or create a valid
        directory. Returns True if the directory is succesfully created.

        Arguments:
        > None. Acts on the current value set on the "location" attribute

        Positional returns:
        > [0]: [bool] Value indicating that location exists
        > [1]: [bool] Value indicating that location was created
        '''
        # First ensure the location is a path or path-like object
        try:
            self.location=os.path.abspath(self.location)
        except:
            return False,False
        
        # Next, check if it exists, and if not, attempt to create it
        if not Path(self.location).is_dir():
            self.notifications.warning(
                    f"No database found at {self.location}."
                    )
            self.notifications.message(
                    f"Creating specified directory {self.location}"
                    )
            try:
                Path(self.location).mkdir(parents=True,exist_ok=True)
            except:
                self.notifications.error(
                        f"Unable to create directory. Cannot access database contents!"
                        )
                return False, False
            else:
                return True, True
        else:
            return True, False

    def build(self):
        '''
        Description:
        A method for building the database.
        
        Accesses a comprehensive list of HLA alleles and their protein
        sequences through the IMGT HLA. Uses this data to match allele
        designations with protein sequences and to generate alpha1/alpha2
        domain sequences for alignments and template designations.

        Also accesses list of PDB structure templates for assigning 
        initial peptide conformations in output models through the IEDB. 

        Arguments:
        >None

        Positional returns:
        > [0]: [bool] Indicates if build was successful.
        >Generates 3 files: , build/allele_seq.info,
        and build/allele_list.info
        '''
        # Ensure the database directory exists or was created (in case self.build()
        # is called outside of the constructor)
        locate_success=self.locate()
        if not locate_success[0]:
            self.notifications.error(f"Unable to access database at {self.location}")
            sys.exit(1)

        # Create the /<MHCdatabase>/build and templates directories to store files in for the database build
        build_dir=os.path.join(self.location,"build")
        templates_dir=os.path.join(self.location,"templates")
        models_dir=os.path.join(self.location,"models")
        Path(build_dir).mkdir(parents=True, exist_ok=True)
        Path(templates_dir).mkdir(parents=True,exist_ok=True)
        Path(models_dir).mkdir(parents=True,exist_ok=True)

        # Download allele list and return the data as a pandas DataFrame object:
        fetch_allele_data=self.__fetch_mhc_allele_info()
        allele_data=fetch_allele_data[0]
        allele_name_unique=allele_data["Allele_Name"].unique()

        # Now we build the structural data for the database. Download the necessary data files from IEDB:
        IEDB_STRUCT_DATA=self.IEDB

        # Download the MHC structural database file if not in the build folder already:
        struct_data_tcr_fn=os.path.join(build_dir,"tcr_3d_assays.csv")
        struct_data_mhc_fn=os.path.join(build_dir, "mhc_3d_assays.csv")

        success=False
        if Path(struct_data_tcr_fn).is_file() and Path(struct_data_mhc_fn).is_file():
            success=True
            self.notifications.message(
                f"IEDB structure data files {struct_data_tcr_fn} and {struct_data_mhc_fn} located"
                )
        else:
            self.notifications.message("Downloading IEDB data from https://www.iedb.org...")
            success=True
            try:
                struct_data=requests.get(IEDB_STRUCT_DATA)
                if struct_data.status_code != 200:
                    success=False
            except requests.exceptions.RequestException:
                success=False

            if not success:
                self.notifications.error(
                    f"{IEDB_STRUCT_DATA} is not reachable\nUnable to build database. Terminating execution"
                    )
                return success

            # Retrieve the data in a zip folder and extract it:
            struct_data_zip_fn=os.path.join(build_dir,os.path.basename(IEDB_STRUCT_DATA))
            struct_data_zip=requests.get(IEDB_STRUCT_DATA, stream=True)
            with open(struct_data_zip_fn, mode="wb") as sfile:
                for chunk in struct_data_zip.iter_content(chunk_size=1024):
                    if chunk:
                        sfile.write(chunk)
            self.notifications.message(f"Structure data retrieved and written to {struct_data_zip_fn}")
            with ZipFile(struct_data_zip_fn,"r") as zfile:
                zfile.extractall(path=build_dir)

            # Remove the zip folder and antibody/b-cell receptor epitope structure data
            os.remove(struct_data_zip_fn)
            if Path(os.path.join(build_dir,"bcr_3d_assays.csv")).is_file():
                os.remove(os.path.join(build_dir,"bcr_3d_assays.csv"))

        # Now read the data into pandas dataframes 
        # Filter out duplicate epitopes, keeping only the highest resolution structures
        struct_data_tcr=pd.read_csv(struct_data_tcr_fn,header=1)
        # Clean up the fieldnames
        struct_data_tcr.columns=struct_data_tcr.columns.str.replace(" ","_")
        struct_data_tcr.columns=struct_data_tcr.columns.str.replace("/","-")

        struct_data_mhc=pd.read_csv(struct_data_mhc_fn,header=1)
        # Clean up the fieldnames
        struct_data_mhc.columns=struct_data_mhc.columns.str.replace(" ","_")
        struct_data_mhc.columns=struct_data_mhc.columns.str.replace("/","-")

        # Concatenate structure data from TCRs and MHC only PDB files into a single dataframe for ease of use:
        struct_data_mhc.rename(
            columns={"MHC_PDB_Chain_1":"MHC_PDB_Chain1"},
            inplace=True
            )
        cols=["Epitope_IRI",
              "Epitope_Description",
              "MHC_Allele",
              "MHC_PDB_Chain1",
              "Antigen_PDB_Chain(s)",
              "PDB_ID",
              "Resolution_(Angstrom)"
              ]

        struct_data=pd.concat([
            struct_data_tcr[cols],
            struct_data_mhc[cols]],
            axis=0,
            ignore_index=True,
            )

        # Generate column with epitope IDs
        struct_data["Epitope_ID"]=struct_data["Epitope_IRI"].str.split("/").apply(lambda x:x[-1])
        struct_data.drop(columns=["Epitope_IRI"])

        # Save structures with mouse alleles H2-Db and H2-Kb
        mouse_struct_data=struct_data.loc[
            (struct_data["MHC_Allele"]=="H2-Db") | (struct_data["MHC_Allele"]=="H2-Kb")
            ]

        # Remove structures without a corresponding class 1 HLA allele in the pre-generated list of HLA alleles 
        # (i.e., CD1, MHC class 2, etc)
        struct_data=struct_data[
            struct_data["MHC_Allele"].isin(["{}{}".format("HLA-",i) for i in allele_name_unique])
            ]
        struct_data["MHC_Allele"]=struct_data["MHC_Allele"].str.split("-").apply(lambda x:x[-1])

        # Recombine mouse and human data
        struct_data=pd.concat([
            struct_data,
            mouse_struct_data],
            ignore_index=True
            )
        stats_total=len(struct_data)

        # Remove structures with duplicate epitopes and NCAA that might not be recognized by Rosetta when threading for flexpepdock
        stats_with_NCAA=len(struct_data[struct_data["Epitope_Description"].str.contains("+",regex=False)])
        struct_data=struct_data[~struct_data["Epitope_Description"].str.contains("+",regex=False)]
        allele_list_struct=struct_data["MHC_Allele"].unique().tolist()

        # Remove structures with duplicate epitopes, keeping only the highest resolution structure
        stats_duplicates=0
        for a in allele_list_struct:
            tmp_epitopes_list=struct_data["Epitope_ID"].loc[struct_data["MHC_Allele"]==a].unique().tolist()
            for e in tmp_epitopes_list:
                e_struct=struct_data.loc[(struct_data["MHC_Allele"]==a) & (struct_data["Epitope_ID"]==e)]
                if len(e_struct)>1:
                    stats_duplicates+=len(e_struct.index[e_struct["Resolution_(Angstrom)"]>e_struct["Resolution_(Angstrom)"].min()].tolist())
                    struct_data.drop(e_struct.index[e_struct["Resolution_(Angstrom)"]>e_struct["Resolution_(Angstrom)"].min()].tolist(),inplace=True)
        struct_data.reset_index(inplace=True,drop=True)

        # Print query statistics
        title1="MHC-I structures total"
        title2="Structures with NCAA"
        title3="Structures duplicate epitopes"
        title4="Total MHC-I structures with unique epitopes"
        title5="Unique alleles with PDB structures"
        self.notifications.message(f"{title1:<45} : {stats_total}")
        self.notifications.message(f"{title2:<45} : {stats_with_NCAA}")
        self.notifications.message(f"{title3:<45} : {stats_duplicates}")
        self.notifications.message(f"{title4:<45} : {len(struct_data)}")
        self.notifications.message(f"{title5:<45} : {len(allele_list_struct)}")

        # Iterate through all HLAs in the database. For each allele, search the structure data. 
        # If the allele exists, download it from the PDB to <allele_name>. 
        # Keep only the MHC and peptide chains.

        # This is for updating the database only. Check if the info file is found.
        allele_struct_info_fn=os.path.join(self.location,"allele_struct.info")
        if Path(allele_struct_info_fn).is_file():
            database_info=pd.read_csv(allele_struct_info_fn)
            self.notifications.message(
                f"{allele_struct_info_fn} located. Updating database with new structures..."
                )
        else:
            database_info=pd.DataFrame(columns=struct_data.columns)
            self.notifications.message(
                f"{allele_struct_info_fn} not found. Building database from the ground up..."
                )

        title1="Allele"
        title2="#PDBs"
        title3="Now Downloading..."
        self.notifications.message(
            f"\n{title1:^15}|{title2:^5}|{title3:^15}\n===============|=====|==============="
            )
        download_count=0

        pdbs_already_downloaded=[
            x.replace(".pdb","") for x in os.listdir(templates_dir)
            ]
        for a in allele_name_unique:
            if a in allele_list_struct:
                # If an allele has at least one structure in the PDB, retrieve them
                tmp_pdb=struct_data["PDB_ID"].loc[(struct_data["MHC_Allele"]==a)].tolist()
                for (i,p) in enumerate(tmp_pdb):
                    # Download each PDB file that does not already exist in the database
                    self.notifications.message(
                        f"{a:^15}|{i+1:^5}|{p:^15}| Downloaded:{download_count:<5}Remaining:{len(struct_data)-download_count:<5}", 
                        endline="\r"
                        )
                    if (p.upper() not in database_info["PDB_ID"].tolist()) and (p.upper() not in pdbs_already_downloaded):
                        tmp_fn=PDBList(verbose=False).retrieve_pdb_file(
                            pdb_code=p,
                            file_format="pdb",
                            pdir=templates_dir
                            )
                        # Clean the PDB file: remove HETATMs, 
                        # keep only alpha chain of MHC and peptide
                        tmp_struct=PDBParser(PERMISSIVE=1, QUIET=1).get_structure(p,tmp_fn)[0]
                        chain_MHC=struct_data["MHC_PDB_Chain1"].loc[
                            (struct_data["PDB_ID"]==p)
                            ].tolist()[0]
                        chain_peptide=struct_data["Antigen_PDB_Chain(s)"].loc[
                            (struct_data["PDB_ID"]==p)
                            ].tolist()[0]
                        io=PDBIO()
                        io.set_structure(tmp_struct)
                        tmp_remove=tmp_fn
                        tmp_fn=os.path.join(
                            os.path.dirname(tmp_fn),
                            os.path.basename(tmp_fn)[-8:-4].upper()+".pdb"
                            )
                        # Save the file with the custom overloaded class that will select only for 
                        # the MHC and peptide chains, and remove HETATMs
                        io.save(
                            tmp_fn,
                            chainSelect(chain_MHC,chain_peptide)
                            )
                        os.remove(tmp_remove) # Remove the downloaded file
                    download_count+=1
                    self.notifications.message(f"{a:^15}|{i+1:^5}|{p:^15}|                                ")
        self.notifications.message(
            f"Database built at {self.location}. Error checking PDB files..."
            )

        # Next we error check the PDBs for bad files
        struct_data=remove_bad_pdbs(templates_dir,struct_data,self.notifications)

        # Now save our  file for future reference
        struct_data.to_csv(os.path.join(self.location,"allele_struct.info"),index=False)

        self.notifications.message(
            f"Database info saved to {os.path.join(self.location,'')}"
            )
        self.notifications.notice(
            "Database build complete"
        )
        # Print statistics about the database for user reference
        self.query()

        return True

    def query(self, action="get_alleles", template=""):
        '''
        Description:
        Utility to query information about structures in MHC/peptide database 
        by accessing and returning requested data from .
        Each behavior output is formatted as a tuple, with the second value 
        a boolean indicating success or failure of the requested operation.

        Arguments:
        >action: [str] The type of query. Options include:
            'get_pdbs':     Prints and returns a list of PDBs with templates for 
                            self.allele.
                            Positional returns:
                            [0] [list] PDB IDs for HLAs matching self.allele
                            [1] [bool] Indicates PDB found for self.allele
            'get_epitopes': Returns epitope sequences for self.allele.
                            Positional returns:
                            [0] [dict] <PDB ID>:<epitope sequence> for all PDBs
                            of self.allele
                            [1] [bool] Indicates whether or not at least one
                            epitope was found
            'get_templates':Searches for the user-requested PDB ID and returns
                            database info for that PDB. Set search template
                            with optional 'template' argument.
                            Positional returns:
                            [0] [dict] Database keys:values for <PDB ID>
                            [1] [bool] Indicates if <PDB ID> found 
            'get_alleles':  Prints a formatted table and returns a list of all
                            available alleles with PDB structures in database.
                            Positional returns:
                            [0] [list] HLAs with PDB structures in database
                            [1] [bool] Indicates if self.allele in [0]
        >template: [str] PDB ID of the template to retrieve. Use when action=
        "get_templates"
        '''
        # Define the query action(s):
        available_actions={
            "get_pdbs":False,
            "get_epitopes":False,
            "get_templates":False,
            "get_alleles":False
        }

        # Check to make sure the user has not muted feedback:
        if self.notifications.toprint < self.notifications._message_level_keys["all"]:
            self.notifications.warning(
                "MHC database query initiated, but notification level set to a value other than 'all'. "+
                "Information regarding this query will not be printed."
            )

        # Quick check to make sure the  file exists
        allele_struct_info_fn=os.path.join(self.location,"allele_struct.info")
        if not Path(allele_struct_info_fn).is_file():
            self.notifications.error(
                f"MHCdatabase().query() called, but {allele_struct_info_fn} not located."
                )
            allele_list=[]
            allele_present=False
            return allele_list, allele_present

        # If it exists, get the database information contained within
        allele_struct_info=pd.read_csv(allele_struct_info_fn)

        # Select the requested actions and execute:
        if action not in available_actions.keys():
            self.notifications.warning(f"Selection {action} invalid for MHCdatabase.query()!"+
                                       f" Available actions include {', '.join(available_actions.keys())}"
                                       )
            action="get_allele"
        available_actions[action]=True

        ## GET PDB ##
        # Return a list of PDBs for self.allele if set. 
        # By default, returns PDBs for the allele set in the database instance. 
        # This behavior can be overwritten if MHCdatabase().allele=None
        if available_actions["get_pdbs"]:
            if self.allele==None:
                pdbs=allele_struct_info
                self.notifications.warning(
                    f"No allele value set for database instance."+ 
                    "Retrieving PDBs for all available"
                    )
            else:
                pdbs=allele_struct_info.loc[
                    allele_struct_info["MHC_Allele"]==self.allele
                    ].copy()
            pdbs.sort_values(
                by="Resolution_(Angstrom)",
                ignore_index=True,
                inplace=True,
                ascending=True
                )
            pdb_ids=pdbs["PDB_ID"].tolist()
            if len(pdb_ids)==0:
                self.notifications.message(
                    f"No PDBs found in database for {self.allele}"
                    )
            else:
                self.notifications.message(
                    "{:^15}{:^20}{:^10}{:^20}".format(
                        "PDBs",
                        "Epitope",
                        "Length",
                        "Resolution (A)"
                        )
                    )
                for pdb in pdb_ids:
                    tmp_epitope=pdbs['Epitope_Description'].loc[
                        pdbs['PDB_ID']==pdb
                        ].values[0]
                    tmp_length=len(tmp_epitope)
                    tmp_res=pdbs['Resolution_(Angstrom)'].loc[
                        pdbs['PDB_ID']==pdb
                        ].values[0]
                    self.notifications.message(
                        f"{pdb:^15}|{tmp_epitope:<20}|{tmp_length:^10}|{tmp_res:^20}"
                        )

            if len(pdb_ids)==0:
                pdb_found=False
            else:
                pdb_found=True

            return pdb_ids, pdb_found

        ## GET EPITOPES ##
        # Return a dictionary with key:value pairs corresponding to PDB:epitope 
        # sequences for the requested allele. If database instance has no allele 
        # set, will return a dict with all PDB and all epitope sequences.
        if available_actions["get_epitopes"]:
            if self.allele==None:
                pdbs=allele_struct_info
                epitopes=allele_struct_info["Epitope_Description"].tolist()
                pdb_ids=allele_struct_info["PDB_ID"].tolist()
                self.notifications.warning(
                    f"No allele value set for database instance. "+
                    "Retrieving epitope sequences for all available."
                    )
            else:
                pdbs=allele_struct_info.loc[
                    allele_struct_info["MHC_Allele"]==self.allele
                    ].copy()
                epitopes=allele_struct_info["Epitope_Description"].loc[
                    allele_struct_info["MHC_Allele"]==self.allele
                    ].tolist()
                pdb_ids=allele_struct_info["PDB_ID"].loc[
                    allele_struct_info["MHC_Allele"]==self.allele
                    ].tolist()

            if len(pdb_ids)==0:
                self.notifications.message(
                    f"No epitope sequences found in database for {self.allele}"
                    )
            else:
                self.notifications.message(
                    "{:^15}{:^20}{:^10}".format(
                        str(self.allele)+" PDBs",
                        "Epitope",
                        "Length"
                        )
                    )
                for pdb in pdb_ids:
                    tmp_epitope=pdbs['Epitope_Description'].loc[
                        pdbs['PDB_ID']==pdb
                        ].values[0]
                    tmp_length=len(tmp_epitope)
                    self.notifications.message(
                        f"{pdb:^15}|{tmp_epitope:^20}|{tmp_length:^10}"
                        )

            if len(pdb_ids)==0:
                epitopes_found=False
            else:
                epitopes_found=True

            return dict(zip(pdb_ids,epitopes)), epitopes_found

        ## GET TEMPLATES ##
        # Searches for a user-requested template based on the PDB ID
        # provided in the 'template' argument. Returns a dictionary
        # of values correpsonding to template data (including the HLA
        # allele, resolution, and epitope sequence)
        if available_actions["get_templates"]:
            if template in allele_struct_info["PDB_ID"].tolist():
                return dict(zip(
                    allele_struct_info.columns.tolist(),
                    allele_struct_info.loc[
                        allele_struct_info["PDB_ID"]==template
                        ].values.tolist()[0]
                        )), True
            else:
                self.notifications.message(
                    f"Template {template} not found in database\n"+
                    "Available templates include:"
                    )
                allele_struct_info.sort_values(
                    by=["MHC_Allele"],
                    inplace=True,
                    ignore_index=True
                    )
                pdb_ids=allele_struct_info["PDB_ID"].tolist()
                self.notifications.message("{:^10}|{:^10}".format("Allele","PDB"))
                self.notifications.message("---------------------")
                for pdb in pdb_ids:
                    mhc_match=allele_struct_info['MHC_Allele'].loc[
                        allele_struct_info['PDB_ID']==pdb
                    ].values[0]
                    self.notifications.message(
                        f"{mhc_match:^10} {pdb}"
                        )
                return dict(zip(
                    allele_struct_info.columns.tolist(),
                    [None]*len(allele_struct_info.columns.tolist())
                    )), False

        ## GET ALLELES ##
        # Returns a list of alleles available in the database, and a boolean
        # value indicating if self.allele has available structural templates.
        if available_actions["get_alleles"]:
            # Print basic statistics about the database:
            allele_list=allele_struct_info["MHC_Allele"].sort_values().unique().tolist()

            self.notifications.message(
                "{:^50}".format("==== MHC-1 Database @"+self.location+"====\n")
                )
            # Generate statistics for the various HLAs
            A_alleles=sum("A" in s for s in allele_list)
            A_count=sum(
                    allele_struct_info["MHC_Allele"].str.contains("A").tolist()
                    )
            B_alleles=sum("B" in s for s in allele_list)
            B_count=sum(
                    allele_struct_info["MHC_Allele"].str.contains("B").tolist()
                    )
            C_alleles=sum("C" in s for s in allele_list)
            C_count=sum(
                    allele_struct_info["MHC_Allele"].str.contains("C").tolist()
                    )
            total_alleles=len(allele_list)
            total_PDBs=len(allele_struct_info)
            self.notifications.message(
                "{:^50}".format("-----------Summary Statistics-----------\n")
                )
            self.notifications.message(
                "{:<10}{:^20}{:^25}".format(
                    "Gene",
                    "Unique Alleles",
                    "Structures with Unique Epitopes"
                    )
                )
            self.notifications.message("{:<10}{:^20}{:^25}".format("All",total_alleles,total_PDBs))
            self.notifications.message("{:<10}{:^20}{:^25}".format("HLA-A",A_alleles,A_count))
            self.notifications.message("{:<10}{:^20}{:^25}".format("HLA-B",B_alleles,B_count))
            self.notifications.message("{:<10}{:^20}{:^25}\n".format("HLA-C",C_alleles,C_count))
            self.notifications.message(
                "{:^50}".format(
                    "----------All Available Alleles----------\n"
                    )
                )
            self.notifications.message(
                "{:^10}{:^20}".format(
                    "Allele",
                    "# of Structures"
                    )
                )
            for a in allele_list:
                num_pdbs=len(
                    allele_struct_info["PDB_ID"].loc[
                        allele_struct_info["MHC_Allele"]==a
                        ]
                    )
                self.notifications.message(f"{a:^10}{num_pdbs:^15}")

            if self.allele in allele_list:
                allele_present=True
            else:
                allele_present=False

            return allele_list, allele_present

    def set_allele(self, allele):
        '''
        Description:
        Method to assign an allele to the database instance. At the time of
        allele assignment, also sets PDB templates in priority (i.e. sequence
        similarity) order, and determines an alpha1/alpha2 domain sequence.

        Arguments:
        >allele: [str] The HLA allele in the format <gene>*<allele_num>:<prot>

        Positional returns:
        [0] [bool] Indicates whether or not  was successful.
        '''
        self.allele=allele
        # Retrieved the alpha 1/alpha 2 domain sequences and find the one for this allele
        allele_sequences=self.__fetch_mhc_allele_info()[1]
        allele_list=allele_sequences["Allele_Names"].apply(
                lambda x: x.split()
                ).apply(pd.Series).stack().tolist()
        if not allele or allele not in allele_list:
            # Exit without setting the allele if not found in the list of available HLA sequences
            self.notifications.error(
                f"Failed to set allele to {self.allele}. "+
                "The requested allele has no available sequence!\n"+
                "Make sure allele is formatted <gene>*<allele_number>:<protein>.\n"+
                "Example: A*02:01"
                )
            self.allele=None
            return False
        else:
            # Otherwise, find the sequence (already trimmed)
            self.sequence=allele_sequences["Sequence"].loc[
                allele_sequences["Allele_Names"].str.contains(self.allele,regex=False)
                ].values[0]
            self.notifications.message(
                f"{self.allele} alpha 1 and alpha 2 domain amino acid sequence retrieved:"+
                f"\n{self.sequence}"
            )
            # And assign structural templates based on sequence similarity
            self.templates=self.__set_mhc_templates()
            if self.templates:
                self.notifications.message(
                    "Structural templates identified in sequence identity order:\n{}".format(
                        ", ".join(self.templates)
                        )
                    )
            self.notifications.notice(
                f"Database instance successfully assigned to allele {self.allele}"
            )
            return True

    def __fetch_mhc_allele_info(self):
        '''
        Description:
        Private method to download and return data from the IMGT HLA database,
        including HLA sequences.

        Arguments:
        >none

        Positional returns:
        [0] [DataFrame] allele_list.info
        [1] [DataFrame] allele_seq.info
        '''
        allele_list=self.IMGTHLA
        allele_list_fn=os.path.join(self.location,"build","allele_list.fasta")

        # Check if the default allele list exists and return requested values without update (or proceed with new build)
        seq_data_default_fn=os.path.abspath(os.path.join(self.location,"build","allele_seq.info"))
        allele_data_default_fn=os.path.abspath(os.path.join(self.location,'build','allele_list.info'))
        if Path(allele_data_default_fn).is_file() and Path(seq_data_default_fn).is_file():
            allele_data=pd.read_csv(allele_data_default_fn)
            seq_data=pd.read_csv(seq_data_default_fn)
            return allele_data, seq_data

        ## Download alleles and sequences from IMGT ##
        else:
            # A quick check to make sure the default URL is reachable if it is needed
            success=True
            try:
                get=requests.get(allele_list)
                if get.status_code != 200:
                    success=False
            except requests.exceptions.RequestException:
                success=False
            if not success:
                self.notifications.error(f"Unable to reach default allele list location at\n{allele_list}\nTerminating execution")
                sys.exit(1)

            # Otherwise proceed with the download
            tmp_data=requests.get(allele_list, stream=True)
            with open(allele_list_fn, "wb") as afile:
                for chunk in tmp_data.iter_content(chunk_size=1024):
                    if chunk:
                        afile.write(chunk)
            self.notifications.message(f"Allele fasta downloaded to {allele_list_fn}")

        ## Create allele_list.info ##
        # Now load the data into a pandas dataframe. Limit to A, B, and C genes
        self.notifications.message("Now processing allele sequence list...",endline="")
        self.notifications.message("removing non-class-I MHCs...",endline="")
        with open(allele_list_fn) as allele_fasta:
            allele_name=[]
            allele_seq=[]
            for record in SeqIO.parse(allele_fasta, "fasta"):
                tmp_des=record.description.split(" ")[1]
                tmp_gene=tmp_des.split("*")[0]
                if (tmp_gene=="A" or tmp_gene=="B" or tmp_gene=="C"):
                    allele_name.append(tmp_des)
                    allele_seq.append(str(record.seq))
        allele_data=pd.DataFrame(columns=["Allele","Sequence"])
        allele_data["Allele"]=allele_name
        allele_data["Sequence"]=allele_seq

        # Remove sequence duplicates
        self.notifications.message("removing duplicate MHCs...",endline="")
        allele_data.drop_duplicates(subset="Sequence",inplace=True,ignore_index=True)

        # Remove alleles with N (null), S (secreted), C (cytoplasmic), A (aberrant), and Q (questionable) suffixes
        self.notifications.message(
            "removing alleles with N (null), S (secreted), C (cytoplasmic), A (aberrant), and Q (questionable) suffixes...",
            endline=""
            )
        discard_filter="NSCAQ"
        allele_data["tmp_alleleName"]=allele_data["Allele"].str.split("*").apply(lambda x:x[1])
        allele_data=allele_data[~allele_data.tmp_alleleName.str.contains("|".join(discard_filter))]
        allele_data.drop(columns=["tmp_alleleName"])
        allele_data.reset_index(inplace=True,drop=True)

        # We want only specific proteins. Create a new column that describes the specific HLA protein (i.e., A*02:01). 
        # For all synonymous DNA substitution variants, select the longest sequence length. 
        # We will save the sequence data for creating models to be used in future peptide binding simulations.
        self.notifications.message("removing synonymous DNA substitution variants...")
        allele_data["Allele_Name"]=allele_data["Allele"].str.split(":").apply(lambda x:x[0:2]).apply(lambda x:":".join(x))
        allele_data["tmp_seqLen"]=allele_data["Sequence"].apply(lambda x:abs(len(x)-181))
        allele_name_unique=allele_data.Allele_Name.unique().tolist()

        # If a given allele has synonymous variants, keep only the variant with the sequence closest in length to 181 
        # (length of alpha1 and alpha2 domains)
        for a in allele_name_unique:
            tmp=allele_data.loc[(allele_data["Allele_Name"]==a)].copy()
            if len(tmp)>1: 
                allele_data.drop(tmp.index[tmp["tmp_seqLen"]>tmp["tmp_seqLen"].min()].tolist(),inplace=True)
        allele_data.drop(columns=["tmp_alleleName","tmp_seqLen","Allele"],inplace=True) # Discard unused columns

        # Define full sequences for H2-Kb and H2-Db for inclusion of mouse data
        H2Db_seq="MGAMAPRTLLLLLAAALAPTQTRAGPHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAGGSHTLQQMSGCDLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRTDSPKAHVTHHPRSKGEVTLRCWALGFYPADITLTWQLNGEELTQDMELVETRPAGDGTFQKWASVVVPLGKEQNYTCRVYHEGLPEPLTLRWEPPPSTDSYMVIVAVLGVLGAMAIIGAVVAFVMKRRRNTGGKGGDYALAPGSQSSEMSLRDCKA"
        H2Kb_seq="MVPCTLLLLLAAALAPTQTRAGPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRTDSPKAHVTHHSRPEDKVTLRCWALGFYPADITLTWQLNGEELIQDMELVETRPAGDGTFQKWASVVVPLGKEQYYTCHVYHQGLPEPLTLRWEPPPSTVSNMATVAVLVVLGAAIVTGAVVAFVMKMRRRNTGGKGGDYALAPGSQTSDLSLPDCKVMVHDPHSLA"
        
        # Include mouse data here
        mouse_data=pd.DataFrame(data={"Sequence":[H2Db_seq,H2Kb_seq],"Allele_Name":["H2-Db","H2-Kb"]})
        # Recombine data and save to allele_list.info
        allele_data=pd.concat([allele_data,mouse_data],ignore_index=True)
        allele_data.reset_index(inplace=True,drop=True)
        allele_data.reindex(columns=allele_data.columns.tolist()+["Templates"],copy=False)
        allele_data.to_csv(allele_data_default_fn,index=False)
        self.notifications.message(
            f"Allele list and corresponding sequences generated and saved to {allele_data_default_fn}."+
            f" {len(allele_data)} sequenced alleles retrieved."
            )

        ## Create allele_seq.info ##
        # Here we build a file with non-redudant sequences of the alpha 1 and alpha 2 domains. 
        # The "Allele_Names" column will be a list of alleles that have these sequences.
        # We will also store template (sequence alignment) data in these files (generated ad hoc)
        seq_trim=[]
        self.notifications.message(f"Trimming sequences and saving only non-redundant values...")
        for a in allele_data["Allele_Name"].tolist():
            if a not in ["H2-Db","H2-Kb"]:
                tmpseq=self.__trim_seq(
                    allele_data["Sequence"].loc[allele_data["Allele_Name"]==a].values[0],
                    a[0]
                    )
                seq_trim.append(tmpseq)
            else:
                tmpseq=self.__trim_seq(
                    allele_data["Sequence"].loc[allele_data["Allele_Name"]==a].values[0],
                    a
                    )
                seq_trim.append(tmpseq)
            if len(seq_trim)%10==0 or len(seq_trim)==len(allele_data):
                self.notifications.message(
                            f"Trimmed {len(seq_trim):<5} of {len(allele_data)} sequences",
                            endline="\r"
                            )
        self.notifications.message("")
        allele_data["Sequence_trimmed"]=seq_trim
        seq_trim=allele_data.Sequence_trimmed.unique().tolist()

        # Save unique sequences to a new dataframe, and then find their matching HLAs
        seq_trim=pd.DataFrame(seq_trim,columns=["Sequence"])
        hla_list=[]
        for s in seq_trim["Sequence"].tolist():
            tmp=allele_data["Allele_Name"].loc[allele_data["Sequence_trimmed"]==s].values
            tmp_hla_list=""
            for v in tmp:
                if tmp_hla_list=="":
                    tmp_hla_list=v
                else:
                    tmp_hla_list=tmp_hla_list+" "+v
            hla_list.append(tmp_hla_list)
        seq_trim["Allele_Names"]=hla_list

        seq_trim.to_csv(os.path.join(self.location,"build","allele_seq.info"),index=False)
        self.notifications.notice("Non-redundant alpha1 and alpha2 domain sequences saved to"+
                                  f"{os.path.join(self.location,'build','allele_seq.info')}"
                                  )

        return allele_data, seq_trim

    def __set_mhc_templates(self):
        '''
        Description:
        Private method that lists closest-matching HLA alleles by
        sequence identity to self.allele.

        Arguments:
        >None

        Positional returns:
        [0] [list] List of templates (empty if failed to set) 
        '''
        toMatch=self.allele
        if not self.allele:
            self.notifications.error(
                "No allele set to database instance. Unable to assign templates")
            return []

        # Add a line to handle mouse alleles
        if toMatch in ["H2-Kb","H2-Db"]:
            self.notifications.warning(
                f"No similar alleles to mouse allele {toMatch}"
                )
            return []

        # Retrieve the sequence of the allele we want to find templates for 
        # (set it and trim it with self.() if not done already)
        if self.sequence=="":
            self.notifications.error(
                f"No sequence found for {self.allele}."
                )
            return []

        # Make sure that the  file exists
        allele_struct_info_fn=os.path.join(self.location,"allele_struct.info")
        if not Path(allele_struct_info_fn).is_file():
            self.notifications.error(
                f"{allele_struct_info_fn} not found. Does it exist? If not, call MHCdatabase.build()"
                )
            return []
        else:
            allele_struct_info=pd.read_csv(allele_struct_info_fn)

        # Check and make sure the allele does not already have associated templates
        all_allele_data=self.__fetch_mhc_allele_info()
        allele_list_info=all_allele_data[0]
        allele_seq_info=all_allele_data[1]
        if "Templates" in allele_list_info.columns.tolist():
            hla_templates=allele_list_info["Templates"].loc[
                allele_list_info["Allele_Name"]==self.allele
                ].values[0]
            if not pd.isna(hla_templates):
                self.notifications.message(
                    f"Ranked templates generated previously for {self.allele}"
                    )
                return hla_templates.split(" ")
        # This check is to ensure we have a "Templates" column in our .csv file
        # This will only ever be executed the first time we attempt to create
        # templates from a brand new database build.
        else:
            allele_list_info.reindex(
                columns=allele_list_info.columns.tolist()+["Templates"],
                copy=False
                )

        self.notifications.message(
            f"No ranked templates identified for {self.allele}"
            )

        # Store a list of all available unique HLAs with PDBs in the structure database
        allele_list_struct=allele_struct_info["MHC_Allele"].unique().tolist()

        # Create a ranked list of structures based on sequence similarity. 
        # Do a pairwise alignment for each allele with a structure, and rank 
        # those alignments by score.
        alignments=pd.DataFrame(columns=["Aligned_to","score"])
        a_type=toMatch[0]
        toMatch_seq=allele_seq_info["Sequence"].loc[
            allele_seq_info["Allele_Names"].str.contains(toMatch)
            ].values[0]
        # Create our Biopython aligner object, heavily penalizing indels
        aligner=Align.PairwiseAligner(
            mode="global",
            gap_score=-10,
        )
        # Assign the BLOSUM62 scoring matrix to the alignment object
        aligner.substitution_matrix=substitution_matrices.load("BLOSUM62")
        # Compare it to every allele with at least 1 available structure and score the alignment
        for a2 in allele_list_struct: 
            # A conditional to ensure we only compare HLA-A to A, HLA-B to B, etc
            if a_type in a2:
                self.notifications.message(
                    f"Generating ranked alignments for: {toMatch:<10} vs {a2:<10}",
                    endline="\r"
                    )
                # Obtain the trimmed alpha 1/alpha 2 domain sequence for the allele in the 
                # structure database
                a2_seq=allele_seq_info["Sequence"].loc[
                    allele_seq_info["Allele_Names"].str.contains(a2,regex=False)
                    ].values[0]
                # Now do a global pairwise alignment and extract the best scoring result
                tmp_alignments=aligner.align(toMatch_seq,a2_seq)
                tmp_alignments=list(tmp_alignments)
                ava2_best=tmp_alignments[0]
                a_scores_best=ava2_best.score
                if len(tmp_alignments)>1:
                    for each_alignment in tmp_alignments:
                        if each_alignment.score<a_scores_best:
                            ava2_best=each_alignment
                            a_scores_best=ava2_best.score
                alignments=pd.concat(
                    [alignments,
                     pd.DataFrame(
                        data={"Aligned_to":[a2],"score":[a_scores_best]},
                        columns=["Aligned_to","score"])
                    ],
                    ignore_index=True
                )
        self.notifications.message("")
        # Sort templates by alignment score, with the highest score at index=0
        alignments.sort_values(
            by="score",
            ascending=False,
            inplace=True
            ) 
        # Write the data to allele_list.info in the /build folder so we don't have to do this again. Should speed up runtime for future applications, and with this approach we don't have to build templates for 13,000 HLAs (which would take...about 4.5 days of runtime, and then would have to be repeated if the database got updated...)
        allele_list_info.at[
            allele_list_info.loc[allele_list_info["Allele_Name"]==self.allele].index[0],
            "Templates"
            ]=" ".join(alignments["Aligned_to"].tolist())
        # re-write the data file if updated
        allele_list_info.to_csv(
            os.path.join(self.location,"build","allele_list.info"),
            index=False
            ) 

        # Return the list of templates
        return alignments["Aligned_to"].tolist()

    def __trim_seq(self,seq,HLA_allele):
        '''
        Description:
        Private method that takes in an HLA sequence and returns 
        the alpha 1 and alpha 2 domains only, based on a simple 
        sliding window alignment of the canonical domain sequences 
        for each allele and the allele in question.

        Arguments:
        >seq: [str] The full HLA sequence to trim
        >HLA_allele: [str] The HLA allele to trim (A, B, C, H2-Db, or H2-Kb)
        '''
        canon_domains={"A":"GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQIMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAAHEAEQLRAYLDGTCVEWLRRYLENGKETLQRT",
                       "B":"GSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDRNTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQSMYGCDVGPDGRLLRGHDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGECVEWLRRYLENGKDKLERA",
                       "C":"CSHSMRYFDTAVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQADRVSLRNLRGYYNQSEDGSHTLQRMSGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKLEAARAAEQLRAYLEGTCVEWLRRYLENGKETLQRA",
                       "H2-Db":"GPHSMRYFETAVSRPGLEEPRYISVGYVDNKEFVRFDSDAENPRYEPRAPWMEQEGPEYWERETQKAKGQEQWFRVSLRNLLGYYNQSAGGSHTLQQMSGCDLGSDWRLLRGYLQFAYEGRDYIALNEDLKTWTAADMAAQITRRKWEQSGAAEHYKAYLEGECVEWLHRYLKNGNATLLRT",
                       "H2-Kb":"GPHSLRYFVTAVSRPGLGEPRYMEVGYVDDTEFVRFDSDAENPRYEPRARWMEQEGPEYWERETQKAKGNEQSFRVDLRTLLGYYNQSKGGSHTIQVISGCEVGSDGRLLRGYQQYAYDGCDYIALNEDLKTWTAADMAALITKHKWEQAGEAERLRAYLEGTCVEWLRRYLKNGNATLLRT"
                       }
        
        canon_domain=canon_domains[HLA_allele]
        # Return hardcoded mouse alpha 1 and alpha 2 domains
        if HLA_allele in ["H2-Kb","H2-Db"]:
            return canon_domain

        # For now we will assume that if the sequence passed is less than the length of the canonical sequence, 
        # it contains both alpha 1 and alpha 2 domains (this appears to be true based on a review of the curated AA sequences)
        len_diff=len(seq)-len(canon_domain)
        if len_diff<=0:
            return seq

        # Otherwise do a simple sliding window of the canonical alpha1/alpha2 sequence vs the allele in question:
        # First create two lists of the same length with no overlapping characters:
        canon_domain_win=list(canon_domain)+[0]*(len(seq))
        seq_win=[0]*len(canon_domain)+list(seq)
        scores=[]
        alignments=[]
        # Now iterate the canonical sequence window over the entire length of the sequence in question, scoring at each position. The maximum score should correspond to the best possible alignment of the sequence with the canonical alpha 1/alpha 2 domains
        for ioffset in range(len(canon_domain_win)):
            tmpscore=0
            for (i,aa) in enumerate(seq_win):
                if aa!=0 and aa==canon_domain_win[i]:
                    tmpscore+=1
            scores.append(tmpscore)
            alignments.append([canon_domain_win,seq_win])
            canon_domain_win.pop()
            canon_domain_win=[0]+canon_domain_win

        # Now trim the sequence of the best alignment to fit the length of the canonical domains
        best_alignment=alignments[scores.index(max(scores))]
        newseq=""
        for i in range(len(best_alignment[0])):
            if best_alignment[0][i]!=0 and best_alignment[1][i]!=0:
                newseq=newseq+best_alignment[1][i]

        return newseq

    def get_peptide_template(self, query_sequence, method="best", cutoff=0, omit=[]):
        '''
        Description:
        Method to retrieve PDB from structural database with a peptide backbone
        sequence most closely matching a query sequence.

        Arguments:
        >query_sequence: [str] AA sequence of sequence to find matching template
        >method: [str][optional] Scoring method to use for match
            "best": [default] Finds the best matching template by alignment score
            "worst": As above, but returns the template with the worst score
            "bestabove": Finds the best matching template by % identity above <cutoff> 
            "worstbelow": Finds the worst matching template by % identity below <cutoff>    
        >cutoff: [int][optional] Score cutoff threshold to accept/decline templates [0-1]
        >omit: [str][optional] List of PDB codes to omit. A value of "self" indicates
            any template with a 100% sequence identify match will be ignored.

        Positional returns:
        [0] [str] 4-character PDB ID of selected template. Empty string if not found.
        '''
        # List of allowed values for the method variable
        methods=[
            "best",
            "worst",
            "bestabove",
            "worstbelow"
            ] 
        # Do some basic error checking of our inputs:
        # Ensure the database instance has an assigned allele:
        if not self.allele:
            self.notifications.error(
                "No allele assigned to database instance!"+
                "Unable to retrieve PDB template."
            )
            return ""

        # Make sure the instance has MHC templates set
        self.__set_mhc_templates()
        # Set methods to default value if user makes an invalid selection
        if method not in methods:
            self.notifications.warning(
                f"{method} is an invalid method."+
                f" Options include: {', '.join(methods)}."+
                " Setting to default 'best'"
                )
            method="best"
        elif (method=="bestabove" or method=="worstbelow") and not 0<=cutoff<=1:
            bad_cutoff=cutoff
            if method=="bestabove":
                cutoff=0
            else:
                cutoff=1
            self.notifications.warning(
                f"{bad_cutoff} is an invalid threshold (must be between 0 and 1)."+
                f" Setting to default of {cutoff}..."
                )
        # Read in our database data file
        db_records=pd.read_csv(os.path.join(self.location,"allele_struct.info")) 
        found_templates=False
        # Ensure we don't accidentally count placeholder 3-letter codes for sequences with a NCAA
        for allele in self.templates:
            db_records=db_records.loc[
                db_records["MHC_Allele"]==allele
                ].copy()
            db_records["Seq_Length"]=db_records["Epitope_Description"].apply(lambda x:len(x))
            #Define our PDBs to omit if the user requests omitting identical peptides
            if "self" in omit: 
                omit_self=db_records["PDB_ID"].loc[
                    db_records["Epitope_Description"]==query_sequence
                    ].tolist()
            else:
                omit_self=["NONE"] # Arbitrary value that will never equal "equal_seq_pdbs"
            equal_seq_lengths=db_records["Seq_Length"].loc[
                db_records["Seq_Length"]==len(query_sequence)
                ].any()
            equal_seq_pdbs=db_records["PDB_ID"].loc[
                db_records["Epitope_Description"]==query_sequence
                ].tolist()
            # Proceed with template selection from this allele only if:
            #    (1) There exists at least one template with the same length as the query sequence, and
            #    (2) There exists at least one equal length template that is not in omit_self
            if omit_self != equal_seq_pdbs and equal_seq_lengths:
                self.notifications.message(
                    f"Found acceptable peptide templates using allele {allele}"
                    )
                found_templates=True
                break
        if not found_templates:
            self.notifications.error(
                f"No suitable peptide backbone templates found for {query_sequence}"+
                " in the MHC database. Unable to thread peptide and generate starting model."+
                " Consider using FlexPepDock ab-initio.")
            return ""
        
        # Append our list of PDBs to omit with the one containing the query sequence if requested by user
        template_seqs=db_records["Epitope_Description"].tolist()
        if "self" in omit:
            omit.extend(
                db_records["PDB_ID"].loc[
                    db_records["Epitope_Description"]==query_sequence
                    ].tolist()
                )
            omit.remove("self")
            if len(omit)!=0:
                self.notifications.message(
                    f"Omitting PDB templates {' '.join(omit)}"
                    )
        omit=[x.lower() for x in omit]

        # Create our Biopython aligner object, heavily penalizing indels
        aligner=Align.PairwiseAligner(
            mode="global",
            gap_score=-10,
        )
        # Assign the BLOSUM62 scoring matrix to the alignment object
        aligner.substitution_matrix=substitution_matrices.load("BLOSUM62")
        # Set starting thresholds for methods that score alignments
        if method=="best":
            thresh=-10000
        elif method=="worst":
            thresh=10000
        else:
            thresh=cutoff
        no_match_found=True
        scores=[]
        self.notifications.message(
            "{:^20}|{:^20}|{:^10}|{:^10}".format(
                "Query sequence",
                "Database match",
                "PDB",
                "Score"
                )
            )
        for record in template_seqs:
            # Do an alignment with the next record 
            # (exclude matches where the template sequence is shorter than the input sequence)
            tmp_pdb=db_records["PDB_ID"].loc[
                db_records["Epitope_Description"]==record
                ].iloc[0]
            # First check: does the query sequence equal the template sequence length, 
            # and is the template not excluded from comparison?
            if len(record)==len(query_sequence) and tmp_pdb.lower() not in omit:
                # If not then move on to the next step on a per method basis
                if method=="best" or method=="worst":
                    tmp_alignments = aligner.align(query_sequence, record)
                    tmp_alignments=list(tmp_alignments) 
                    for a in tmp_alignments:
                        # For each alignment, find the lowest score 
                        # and save it if less than current best score 
                        # (or reverse if method set to "best")
                        scores.append(a.score)
                        if (method=="worst" and a.score<thresh) or (method=="best" and a.score>thresh):
                            best_sequence=record
                            best_alignment=a
                            thresh=a.score
                            best_pdb=tmp_pdb
                            self.notifications.message(
                                f"{query_sequence:^20}|"+
                                f"{best_sequence:^20}|"+
                                f"{best_pdb:^10}|"+
                                f"{thresh:^10}"
                                )
                            if no_match_found:
                                no_match_found=False
                else:
                    matching_res_count=0
                    for i,aa in enumerate(query_sequence):
                        if aa==record[i]:
                            matching_res_count+=1
                    percent_match=matching_res_count/len(query_sequence)
                    scores.append(percent_match)
                    if (method=="bestabove" and percent_match>thresh) or (method=="worstbelow" and percent_match<thresh):
                        best_sequence=record
                        thresh=percent_match
                        best_pdb=tmp_pdb
                        if no_match_found:
                            no_match_found=False
        if no_match_found:
            if len(scores)==0:
                bestscore="None"
            elif method=="best" or method=="bestabove":
                bestscore=max(scores)
            else:
                bestscore=min(scores)
            self.notifications.warning(
                f"No satisfactory alignment found for input sequence {query_sequence}."+
                f"Adjust thresh to a value less than an acceptable minimum score"+
                f"(best alignment score found is {bestscore}, method chosen is {method})"
                )
            return False
        else:
            self.notifications.message(
                f"\nTarget sequence:{query_sequence:>15}"
                )
            self.notifications.message(
                f"Database match: {best_sequence:>15}"
                )
            self.notifications.message(
                f"\nAlignment:\n{best_alignment}"
            )
        best_pdb_fn = os.path.abspath(
            os.path.join(
                self.location,
                "templates",
                best_pdb+".pdb"
                )
            )

        # Make sure we have the matching PDB for threading the sequence, otherwise abort
        if not os.path.isfile(best_pdb_fn):
            self.notifications.error(
                f"Matching PDB file {best_pdb} not found in database"
                )
            return False
        else:
            self.notifications.notice(
                f"Matching PDB found in database: {best_pdb}"
                )
        return best_pdb

    def help(self):
        '''
        Description:
        Method for printing the module docstring and viewing
        module documentation.
        '''
        # Print the docs and exit
        print(__doc__)
        return
