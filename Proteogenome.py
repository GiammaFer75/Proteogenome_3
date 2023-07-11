#!usr/bin/python3
import os

from Proteogenome3 import pg_indexes as pg_i
from Proteogenome3 import pg_input
from Proteogenome3 import pg_output
from Proteogenome3 import pg_utils
from Proteogenome3 import pg_data_cleaning as dc
from Proteogenome3 import pg_data_processing as dp


from absl import app
from absl import flags

import subprocess
from os import chdir, makedirs

import shutil
import pathlib

from Test_PoGo_absl_subprocess import printinput as pi

# ****************************************************************************** #
#                               PATHS DEFINITION                                 #
home_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/'
pogo_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/Proteogenome-master/Proteogenome-master/PoGo/Windows'
# ****************************************************************************** #


FLAGS = flags.FLAGS

flags.DEFINE_string('project_dir', None, 'The output directory')
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

    if FLAGS.project_dir:
        project_dir_path = pathlib.Path(FLAGS.project_dir)                  # Project Directory
    else: project_dir_path = pathlib.Path.cwd()
    pogo_windows_exe_path = pathlib.Path(FLAGS.pogo_windows_exe)            # PoGo Executable
    protein_FASTA_seq_path = pathlib.Path(FLAGS.protein_FASTA_seq)          # FASTA protein sequences
    if FLAGS.protein_GFF3_annots:
        protein_GFF3_annots_path = pathlib.Path(FLAGS.protein_GFF3_annots)  # GFF3 Reference Genome Annotations
    else: protein_GFF3_annots_path = None
    if FLAGS.protein_GTF_annots:
        protein_GTF_annots_path = pathlib.Path(FLAGS.protein_GTF_annots)    # GTF Reference Genome Annotations
    else: protein_GTF_annots_path = None
    peptides_table_path = pathlib.Path(FLAGS.peptides_table)                # Peptides Input Table

    print(project_dir_path)
    print(pogo_windows_exe_path)
    print(protein_FASTA_seq_path)
    print(protein_GFF3_annots_path)
    print(protein_GTF_annots_path)
    print(peptides_table_path)


    # GENERATE WORKING DIRECTORIES
    project_dir_path, pg_output_path, PoGo_input_path, PoGo_output_path = pg_utils.gen_dir_structure(project_dir_path)

    # Set the annotation format flag
    if (protein_GFF3_annots_path != None) and (protein_GTF_annots_path == None): annotations_format = 'gff3'
    elif (protein_GTF_annots_path != None) and (protein_GFF3_annots_path == None): annotations_format = 'gtf'
    else: raise Exception("Apparently have been provided two genome annotations. Please input one annotation only.")

    # UPLOAD INPUT FILES
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
    protein_FASTA_seq_path = PoGo_input_path.joinpath('PoGo_prot_seq.fasta')
    pg_output.make_sep_file(protein_FASTA_seq_path, FASTA_seq_lst, sep=' ')
    protein_GTF_annots_path = PoGo_input_path.joinpath('PoGo_annot.gtf')
    pg_output.make_sep_file(protein_GTF_annots_path, annot_lst, sep=' ')

    # GENERATE PoGo INPUT PEPTIDES FILE and TABLE
    PoGo_input_table_path = PoGo_input_path.joinpath('PoGo_Input_Table.txt')
    PoGo_input_table = pg_input.gen_PoGo_input_table(peptides_input_table, out_file_name=PoGo_input_table_path)

# ********************* REMAINING ROWS FROM PROTEOGENOME2
    prot_CDS_index = {}  # Initialise dictionary for protein ---> CDS index
    prot_PSMint_index = {}  # Initialise dictionary for protein ---> PSM - intensity - RGB intensity index
# *********************

    # GENERATE INDEXES
    CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format=annotations_format)

    # TESTING COMMANDS
    pi.print_input(pogo_windows_exe_path, protein_FASTA_seq_path, protein_GTF_annots_path, PoGo_input_table_path)
    chdir(pogo_windows_exe_path)              # Go in the PoGo executable folder

    PoGo_command =[pogo_windows_exe_path.joinpath('PoGo.exe'),
                   '-fasta', protein_FASTA_seq_path,
                   '-gtf', protein_GTF_annots_path,
                   '-in', PoGo_input_table_path]
    print(f'\n################\n{PoGo_command}\n################\n')
    # PoGo_command = '.\PoGo.exe'
    subprocess.run(PoGo_command, capture_output=True)
    pg_utils.move_files(PoGo_input_path, PoGo_output_path,filenames_patterns=['*.bed','*.gct',
                                                                              '*_unmapped.txt', '*_out.gtf'])

    # subprocess.run(['dir'], shell=True)
    # subprocess.run(['.\PoGo.exe', 'capture_output=True'])


    # GENERATE MAPS

    # Proteins MAP
    proteogenome_peptide_MAP_path = pg_output_path.joinpath('Proteins_MAP.bed')
    pg_output.gen_protein_track(protein_pep_index, prot_CDS_index, bed_fn=proteogenome_peptide_MAP_path)

    # Peptides MAP
    PoGo_peptide_map_path = PoGo_output_path.joinpath('PoGo_Input_Table.bed')
    PoGo_peptide_map_table = pg_input.load_generic_table(PoGo_peptide_map_path)
    proteogenome_peptide_MAP_path = pg_output_path.joinpath('Peptides_MAP.bed')
    print(f'Path for the peptide MAP ---- > {proteogenome_peptide_MAP_path}')
    dp.filter_peptides(PoGo_peptide_map_table, pep_protein_index, prot_CDS_index,out_file_name=proteogenome_peptide_MAP_path)

    # PTMs MAP
    PoGo_PTM_map_path = PoGo_output_path.joinpath('PoGo_Input_Table_ptm.bed')
    PoGo_PTM_map_table = pg_input.load_generic_table(PoGo_PTM_map_path)
    proteogenome_PTM_MAP_path = pg_output_path.joinpath('PTMs_MAP.bed')
    print(f'Path for the PTM MAP ---- > {PoGo_PTM_map_path}')
    dp.filter_peptides(PoGo_PTM_map_table, pep_protein_index, prot_CDS_index, out_file_name=proteogenome_PTM_MAP_path)




if __name__ == '__main__':
    """
    Proteogenome
    
    Version: 3.0.0
    
    Author: Giammarco Ferrari
    """
    app.run(main)