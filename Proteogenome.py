#!usr/bin/python3

from Proteogenome3 import pg_indexes as pg_i
from Proteogenome3 import pg_input
from Proteogenome3 import pg_output
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
flags.DEFINE_string('out_dir', None, 'The output directory')

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
    output_dir_path = FLAGS.out_dir

    print(pogo_windows_exe_path)
    print(protein_FASTA_seq_path)
    print(protein_GFF3_annots_path)
    print(protein_GTF_annots_path)
    print(peptides_table_path)
    print(output_dir_path)
    # Set the annotation format flag
    if (protein_GFF3_annots_path != None) and (protein_GTF_annots_path == None): annotations_format = 'gff3'
    elif (protein_GTF_annots_path != None) and (protein_GFF3_annots_path == None): annotations_format = 'gtf'
    else: raise Exception("Apparently have been provided two genome annotations. Please input one annotation only.")

    # Upload input files
    print('\n')
    print("UPLOAD INPUT FILES".center(shutil.get_terminal_size().columns))
    peptides_input_table = pg_input.load_input_table(peptides_table_path)       # Upload peptides table
    FASTA_seq_lst = pg_input.load_protein_FASTA_seq(protein_FASTA_seq_path)     # Upload protein FASTA sequences
    annot_lst = pg_input.file_to_lst(protein_GTF_annots_path)                   # Upload reference genome annotations

    # FASTA rectification
    FASTA_seq_lst = dc.rectify_rows(FASTA_seq_lst, target_sub_str=[('lcl|NC_006273.2_prot_', '')])
    FASTA_seq_lst = dc.rectify_rows(FASTA_seq_lst, target_patterns=[('gene=.*?\s', '')])
    FASTA_seq_lst = dc.locus_tag_substitution(FASTA_seq_lst)

    # GTF rectification
    annot_lst = dc.rectify_rows(annot_lst, target_sub_str=[('gene-', ''), ('rna-', '')])
    annot_lst = dc.rectify_rows(annot_lst, open_patterns=[('.*?gene_id\s\"(.*?)\"', '.*?locus_tag\s\"(.*?)\"')])

    # SAVE cleaned files
    protein_FASTA_seq_path = output_dir_path + 'PoGo_prot_seq.fasta'
    pg_output.make_sep_file(protein_FASTA_seq_path, FASTA_seq_lst, sep=' ')
    protein_GTF_annots_path = output_dir_path + 'PoGo_annot.gtf'
    pg_output.make_sep_file(protein_GTF_annots_path, annot_lst, sep=' ')

# ********************* REMAINING ROWS FROM PROTEOGENOME2
    prot_CDS_index = {}  # Initialise dictionary for protein ---> CDS index
    prot_PSMint_index = {}  # Initialise dictionary for protein ---> PSM - intensity - RGB intensity index
# *********************

    # Generate indexes
    CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format=annotations_format)

    #TESTING COMMANDS
    PoGo_input_table_path = output_dir_path + 'PoGo_Input_Table.txt'
    print(f'Pogo input table file - {PoGo_input_table_path}')
    PoGo_input_table = pg_input.gen_PoGo_input_table(peptides_input_table, out_file_name=PoGo_input_table_path)
    pi.print_input(pogo_windows_exe_path, protein_FASTA_seq_path, protein_GTF_annots_path, PoGo_input_table_path)
    chdir(pogo_windows_exe_path)

    PoGo_command =[pogo_windows_exe_path+'PoGo.exe', '-fasta', protein_FASTA_seq_path, '-gtf', protein_GTF_annots_path,
                   '-in', PoGo_input_table_path]
    print(f'################\n{PoGo_command}\n################')
    # PoGo_command = '.\PoGo.exe'
    subprocess.run(PoGo_command, capture_output=True)
    # subprocess.run(['dir'], shell=True)
    # subprocess.run(['.\PoGo.exe', 'capture_output=True'])


    #u.print_lst(FASTA_seq_lst)
    # u.print_lst(annot_lst, limit=20)



if __name__ == '__main__':
    """
    Proteogenome
    
    Version: 3.0.0
    
    Author: Giammarco Ferrari
    """
    app.run(main)