#!usr/bin/python3

from Proteogenome3 import pg_indexes as pg_i
from Proteogenome3 import pg_input
from Proteogenome3 import pg_utils as u
from Proteogenome3 import pg_data_cleaning as dc

from absl import app
from absl import flags

import subprocess
from os import chdir

import shutil

from Test_PoGo_absl_subprocess import printinput as pi

# ****************************************************************************** #
#                               PATHS DEFINITION                                 #
home_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/'
pogo_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/Proteogenome-master/Proteogenome-master/PoGo/Windows'
# ****************************************************************************** #


FLAGS = flags.FLAGS

flags.DEFINE_string('pogo_windows_exe', None, 'The executable PoGo file for running the software on Windows OS.')
flags.DEFINE_string('protein_FASTA_seq', None, 'The protein sequences related to the reference organism genome.')
flags.DEFINE_string('protein_GFF3_annots', None, 'The protein annotations related to the reference organism genome in GFF3.')
flags.DEFINE_string('protein_GTF_annots', None, 'The protein annotations related to the reference organism genome in GTF.')
flags.DEFINE_string('peptides_table', None, 'The table containing the peptides sequences with the format:\n'
                                            '- Protein ID\n'
                                            '- Peptide sequence\n'
                                            '- Peptide modification - PTM\n'
                                            '- Peptide PSM\n'
                                            '- Peptide intensity\n')

def main(argv):

    # ***** FLAGS DEFINITION ***** #
    print("-----------------------------------------------------".center(shutil.get_terminal_size().columns))
    print("START PROTEOGENOM".center(shutil.get_terminal_size().columns))
    print("-----------------------------------------------------".center(shutil.get_terminal_size().columns))
    pogo_windows_exe_path = FLAGS.pogo_windows_exe
    protein_FASTA_seq_path = FLAGS.protein_FASTA_seq
    protein_GFF3_annots_path = FLAGS.protein_GFF3_annots
    protein_GTF_annots_path = FLAGS.protein_GTF_annots
    peptides_table_path = FLAGS.peptides_table

    print(protein_GFF3_annots_path)
    print(protein_GTF_annots_path)
    # Set the annotation format flag
    if (protein_GFF3_annots_path != None) and (protein_GTF_annots_path == None): annotations_format = 'gff3'
    elif (protein_GTF_annots_path != None) and (protein_GFF3_annots_path == None): annotations_format = 'gtf'
    else: raise Exception("Apparently have been provided two genome annotations. Please input one annotation only.")

    # Upload input files
    print('\n')
    print("UPLOAD INPUT FILES".center(shutil.get_terminal_size().columns))
    peptides_input_table = pg_input.load_input_table(peptides_table_path)       # Upload peptides table
    FASTA_seq_lst = pg_input.file_to_lst(protein_FASTA_seq_path)                # Upload FASTA file
    FASTA_chr_to_remove, compact_FASTA = dc.check_FASTA_format(FASTA_seq_lst)   # Find undesired characters in the FASTA
    if FASTA_chr_to_remove:
        print('The input FASTA format is not compliant with PoGo requirements.\nStarting cleaning FASTA')
        # u.print_lst(FASTA_seq_lst)
        char_to_remove = [(ch, '') for ch in FASTA_chr_to_remove]
        FASTA_seq_lst = dc.rectify_rows(FASTA_seq_lst, target_sub_str=char_to_remove) # Remove undesired chars from FASTA
        print('Undesired characters removed')
    if compact_FASTA:
        print('----------- BEFORE')
        # u.print_lst(FASTA_seq_lst, limit=100)
        FASTA_seq_lst = dc.FASTA_cpt_seq(FASTA_seq_lst)
        print('----------- AFTER')
        # u.print_lst(FASTA_seq_lst,limit=100)
        a=input()
    else: print('FASTA content is PoGo compliant')



# ********************* REMAINING ROWS FROM PROTEOGENOME2
    prot_CDS_index = {}  # Initialise dictionary for protein ---> CDS index
    prot_PSMint_index = {}  # Initialise dictionary for protein ---> PSM - intensity - RGB intensity index
# *********************

    # Generate indexes
    CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format=annotations_format)

    #TESTING COMMANDS
    pi.print_input(pogo_windows_exe_path, protein_FASTA_seq_path, protein_GTF_annots_path, peptides_table_path)
    chdir(pogo_dir)
    subprocess.run(['.\PoGo.exe', 'capture_output=True'])
    subprocess.run(['dir'], shell=True)




if __name__ == '__main__':
    """
    Proteogenome
    
    Version: 3.0.0
    
    Author: Giammarco Ferrari
    """
    app.run(main)