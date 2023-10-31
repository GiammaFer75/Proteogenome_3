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
import json

from Test_PoGo_absl_subprocess import printinput as pi

# ****************************************************************************** #
#                               PATHS DEFINITION                                 #
home_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/'
pogo_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/Proteogenome-master/Proteogenome-master/PoGo/Windows'
install_folder = 'C:/Bioinformatics/Proteogenome_v3/Proteogenome_Install'
# ****************************************************************************** #



FLAGS = flags.FLAGS

flags.DEFINE_string('install_folder', None, 'The installation folder')

flags.DEFINE_string('species', None, 'The species to map')

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
    print("START PROTEOGENOME".center(shutil.get_terminal_size().columns))
    print("-----------------------------------------------------".center(shutil.get_terminal_size().columns))

    if FLAGS.project_dir:
        project_dir_path = pathlib.Path(FLAGS.project_dir)                  # Project Directory
    else:
        print(f'{FLAGS.project_dir} - NOT FOUND')
        project_dir_path = pathlib.Path.cwd()

    # INPUT WITH JSON FILE
    if pg_utils.check_input_json(project_dir_path):
        pg_install_folder_path = pathlib.Path(install_folder) # Now the installation folder is static
        print('Input JSON file detected')
        print(f'Project Directory : {project_dir_path}')
        pogo_windows_exe_path, protein_FASTA_seq_path, protein_GFF3_annots_path, protein_GTF_annots_path, \
        peptides_table_path, species = pg_utils.input_json(project_dir_path)
        #print('PATH - ',pogo_windows_exe_path, protein_FASTA_seq_path, protein_GFF3_annots_path, protein_GTF_annots_path, peptides_table_path, species )
    else:
        # INPUT WITH FLAGS
        if not FLAGS.install_folder:
            pg_install_folder_path = pathlib.Path(install_folder)               # Now the installation folder is static

        pogo_windows_exe_path = pathlib.Path(FLAGS.pogo_windows_exe)            # PoGo Executable
        protein_FASTA_seq_path = pathlib.Path(FLAGS.protein_FASTA_seq)          # FASTA protein sequences
        if FLAGS.protein_GFF3_annots:
            protein_GFF3_annots_path = pathlib.Path(FLAGS.protein_GFF3_annots)  # GFF3 Reference Genome Annotations
        else: protein_GFF3_annots_path = None
        if FLAGS.protein_GTF_annots:
            protein_GTF_annots_path = pathlib.Path(FLAGS.protein_GTF_annots)    # GTF Reference Genome Annotations
        else: protein_GTF_annots_path = None
        peptides_table_path = pathlib.Path(FLAGS.peptides_table)                # Peptides Input Table

        species = FLAGS.species                                                 # Species to process

    print(f'Install folder - {pg_install_folder_path}')

    print('species                  - ',species)
    print('project_dir_path         - ',project_dir_path)
    print('pogo_windows_exe_path    - ',pogo_windows_exe_path)
    print('protein_FASTA_seq_path   - ',protein_FASTA_seq_path)
    print('protein_GFF3_annots_path - ',protein_GFF3_annots_path)
    print('protein_GTF_annots_path  - ',protein_GTF_annots_path)
    print('peptides_table_path      - ',peptides_table_path)

    # GENERATE WORKING DIRECTORIES
    project_dir_path, pg_output_path, pg_data_structure_path, PoGo_input_path, PoGo_output_path = pg_utils.gen_dir_structure(project_dir_path)

    # Set the annotation format flag
    if (protein_GFF3_annots_path != None) and (protein_GTF_annots_path == None): annotations_format = 'gff3'
    elif (protein_GTF_annots_path != None) and (protein_GFF3_annots_path == None): annotations_format = 'gtf'
    else: raise Exception("Apparently have been provided two genome annotations. Please input one annotation only.")

    # UPLOAD INPUT FILES
    print('\n')
    print("UPLOAD INPUT FILES".center(shutil.get_terminal_size().columns))

    peptides_input_table = pg_input.load_input_table(peptides_table_path)       # Upload peptides table
    FASTA_seq_lst = pg_input.load_protein_FASTA_seq(protein_FASTA_seq_path)     # Upload protein FASTA sequences
    if protein_GTF_annots_path and not protein_GFF3_annots_path:                # Upload reference genome annotations
        annot_lst = pg_input.file_to_lst(protein_GTF_annots_path)
    elif protein_GFF3_annots_path and not protein_GTF_annots_path:
        annot_lst = pg_input.file_to_lst(protein_GFF3_annots_path)
    else:
        print('Please check your annotation file path')
        exit

    uniprot_to_ensembl_acc=[]
    # CONVERT UniProt accession in Ensemble accession in protein_pep_index
    if species == 'homo sapiens':
        input_protein_list = peptides_input_table[:,0]
        # Conversion table DO NOT EXISTS
        if not pg_utils.check_input_json(pg_install_folder_path, proteogenome_input_json = 'uniprot_to_ensembl_acc.json'):
            print('Protein conversion JSON file - NOT FOUND\nSTART the conversion of UniProt accession codes to Ensembl accession codes\n')
            uniprot_to_ensembl_acc = pg_utils.UniProt2Ensembl(accession_lst=input_protein_list,
                                                              speces='homo sapiens')  # Dictionary of UniProt -> Ensemble acc
            # Save the conversion table from uniprot to Ensembl in a JSON file
            f = open(pg_install_folder_path.joinpath('uniprot_to_ensembl_acc.json'), 'w')
            json.dump(uniprot_to_ensembl_acc, f)
            f.close()
        # Conversion Table already EXIST
        else:
            print(f'Protein conversion JSON file - FOUND\nUpload converted codes')
            f = open(pg_install_folder_path.joinpath('uniprot_to_ensembl_acc.json'), 'r')
            uniprot_to_ensembl_acc = json.load(f)
            uniprot_to_ensembl_acc_keys = uniprot_to_ensembl_acc.keys()  # Extract UniProt codes from the conversion table for comparison
            print(f'{len(uniprot_to_ensembl_acc_keys)} already converted\n')
            f.close()
            uniprot_not_matched = [uni_pro_acc for uni_pro_acc in input_protein_list if uni_pro_acc not in uniprot_to_ensembl_acc_keys]
            if uniprot_not_matched:
                print(f'Detected new UniProt accession codes not yet converted\n-----------------------\n{uniprot_not_matched}\nSTART conversion of unmatched codes')
                uniprot_not_matched_converted = pg_utils.UniProt2Ensembl(accession_lst=uniprot_not_matched,
                                                                         speces='homo sapiens')  # Dictionary of new UniProt -> Ensemble acc
                uniprot_to_ensembl_acc.update(uniprot_not_matched_converted)
                print('Update conversion table')
                pg_output.save_json_file(pg_install_folder_path.joinpath('uniprot_to_ensembl_acc.json'), uniprot_to_ensembl_acc)


        # Update the UniProt accession with Ensemble accession (in GENCODE they use sequences reference to Ensembl)
        for uniprot_acc, ensembl_acc in uniprot_to_ensembl_acc.items():
            if ensembl_acc: peptides_input_table[peptides_input_table==uniprot_acc] = ensembl_acc
        print(f'\npeptides_input_table converted in Ensembl\n--------------------')
        for i in peptides_input_table: print(i)

        # GENERATE INDEXES - for homo sapiens
        CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = \
        pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format='gtf_compliant')

        # PoGo version for humans
        pogo_version = 'PoGo_windows_v1.0.0.exe'


    elif species == 'virus':

        # GENERATE INDEXES - for virus
        print('Generate indexes for virus')
        CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = \
        pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format='not_gtf_compliant')

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

        pg_output.make_sep_file(pg_data_structure_path.joinpath('CDS_matrix.txt'), CDS_matrix, sep='') # Save the CDS annotation rows in a separate matrix (ONLY for annotation not in GTF)

        # PoGo version for virus
        pogo_version = 'PoGo_windows_v1.2.3.exe'

    # GENERATE PoGo INPUT PEPTIDES FILE and TABLE
    PoGo_input_table_path = PoGo_input_path.joinpath('PoGo_Input_Table.txt')
    PoGo_input_table = pg_input.gen_PoGo_input_table(peptides_input_table, out_file_name=PoGo_input_table_path)

    # prot_CDS_index = {}     # Initialise dictionary for protein ---> CDS index
    prot_PSMint_index = {}  # Initialise dictionary for protein ---> PSM - intensity - RGB intensity index

    # GENERATE INDEXES
    # CDS_matrix, prot_CDS_index, protein_pep_index, pep_protein_index = \
    # pg_i.initialise_indexes(protein_GTF_annots_path, peptides_input_table, annot_format=annotations_format)

    # SAVE INDEXES
    print('Save Data Structures')
    if species != 'homo sapiens': pg_output.make_sep_file(pg_data_structure_path.joinpath('CDS_matrix.txt'),
                                                          CDS_matrix, sep='')
    pg_output.save_dict_list(prot_CDS_index, pg_data_structure_path.joinpath('prot_CDS_index.txt'), 'd')
    pg_output.save_dict_list(protein_pep_index, pg_data_structure_path.joinpath('protein_pep_index.txt'), 'd')
    pg_output.save_dict_list(pep_protein_index, pg_data_structure_path.joinpath('pep_protein_index.txt'), 'd')
    print('Data Structure SAVED\n')

    PoGo_command =[pathlib.Path(pogo_windows_exe_path, pogo_version),
                   '-fasta', protein_FASTA_seq_path,
                   '-gtf', protein_GTF_annots_path,
                   '-in', PoGo_input_table_path]
    print('################'.center(shutil.get_terminal_size().columns))
    print('RUN PoGo'.center(shutil.get_terminal_size().columns))
    print('################'.center(shutil.get_terminal_size().columns))
    print(f'PoGo Command --- \n{PoGo_command}')
    # PoGo_command = '.\PoGo.exe'
    subprocess.run(PoGo_command, capture_output=True)
    pg_utils.move_files(PoGo_input_path, PoGo_output_path,filenames_patterns=['*.bed','*.gct',
                                                                              '*_unmapped.txt', '*_out.gtf'])

    # GENERATE MAPS

    # Proteins MAP
    proteogenome_peptide_MAP_path = pg_output_path.joinpath('Proteins_MAP.bed')
    proteins_not_found = pg_output.gen_protein_track(protein_pep_index, prot_CDS_index,
                                                     species=species, bed_fn=proteogenome_peptide_MAP_path)
    print(f'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nPROTEINS NOT FOUND\n{proteins_not_found}\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')

    # Peptides MAP
    PoGo_peptide_map_path = PoGo_output_path.joinpath('PoGo_Input_Table.bed')
    PoGo_peptide_map_table = pg_input.load_generic_table(PoGo_peptide_map_path)
    proteogenome_peptide_MAP_path = pg_output_path.joinpath('Peptides_MAP.bed')
    print(f'Path for the peptide MAP ---- > {proteogenome_peptide_MAP_path}')
    dp.filter_peptides(PoGo_peptide_map_table, pep_protein_index, prot_CDS_index,out_file_name=proteogenome_peptide_MAP_path)

    # PTMs MAP
    PoGo_PTM_map_path = PoGo_output_path.joinpath('PoGo_Input_Table_ptm.bed')
    if PoGo_PTM_map_path.stat().st_size != 0:                     # Check if there is a PTM map
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