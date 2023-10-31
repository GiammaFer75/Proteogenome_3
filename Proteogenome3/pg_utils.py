import os
import sys
import json
import shutil
import pathlib
import requests


def check_input_json(proteogenome_project_dir, proteogenome_input_json = 'proteogenome_input.json'):
    input_file = pathlib.Path(pathlib.PurePath(proteogenome_project_dir, f'./{proteogenome_input_json}'))
    print('Looking for - ', input_file)
    return input_file.is_file()
def input_json(project_directory):
    """
    Version : 1.0
    Name History : input_json
    The function check for the json input file in the project_directory.
    Returns the parameters found in the json file. They are empty if json is not presents.
    :param  project_directory   Str path to the project directory:
    :return pogo_windows_exe_path:
    :return protein_FASTA_seq_path:
    :return protein_GFF3_annots_path:
    :return protein_GTF_annots_path:
    :return peptides_table_path:
    """

    input_json_path = pathlib.Path(project_directory)

    input_json_path = list(input_json_path.glob("proteogenome_input.json"))
    if not input_json_path:
        print(f'The proteogenome_input.json NOT FOUND in the current project folder')
        pogo_windows_exe_path, protein_FASTA_seq_path, \
            protein_GFF3_annots_path, protein_GTF_annots_path, peptides_table_path, species = '', '', '', '', ''
    else:
        # str(pathlib.PureWindowsPath(list(input_json_path.glob("proteogenome_input.json"))[0]))
        input_json_path = str(pathlib.PureWindowsPath(input_json_path[0]))
        input_json = open(input_json_path, 'r')
        input_json = json.load(input_json)
        input_json_keys = input_json.keys()
        # print(input_json_keys)
        # Fill the field not input by the user
        if not 'species' in input_json_keys: input_json['species'] = None
        if not 'pogo_windows_exe' in input_json_keys: input_json['pogo_windows_exe'] = None
        if not 'protein_FASTA_seq' in input_json_keys: input_json['protein_FASTA_seq'] = None
        if not 'protein_GFF3_annots' in input_json_keys: input_json['protein_GFF3_annots'] = None
        if not 'protein_GTF_annots' in input_json_keys: input_json['protein_GTF_annots'] = None
        if not 'peptides_table' in input_json_keys: input_json['peptides_table'] = None

        species = input_json['species']

        # Upload the input path
        pogo_windows_exe_path = input_json['pogo_windows_exe']
        protein_FASTA_seq_path = input_json['protein_FASTA_seq']
        protein_GFF3_annots_path = input_json['protein_GFF3_annots']
        protein_GTF_annots_path = input_json['protein_GTF_annots']
        peptides_table_path = input_json['peptides_table']

    print('From JSON - {}\n{}\n{}\n{}\n{}\n' \
          .format(pogo_windows_exe_path, protein_FASTA_seq_path, protein_GFF3_annots_path, \
          protein_GTF_annots_path, peptides_table_path))
    return pogo_windows_exe_path, protein_FASTA_seq_path, protein_GFF3_annots_path, protein_GTF_annots_path, \
           peptides_table_path, species


def print_lst(input_list, limit=10, en_sep=True, sep_type='-'):
    """
    Version : 2.0
    Name History : print_lst (from Proteogenome 1.0)

    Print the list content limitely to the element indicated by the parameter 'limit'.

    :param  input_list  List    Items to print.                                                  :
    :param  limit       Int     Number of items to print.                                        :
    :param  en_sept     Bool    Enable the separation of each element of the list with a string. :
    :param  sep_type    Char    The character that will make up the separator string.            :

    :return     Print the data structure content:
    """

    for ind, i in enumerate(input_list):
        if ind < limit:
            print(i)
            if en_sep: print(sep_type * 100)

def create_path(input_filename, add_prefix='', add_suffix=''):
    """
    Version : 1.0
    Name History : create_path

    This function take the absolute path of a file and add a prefix or a suffix to the original filename, based on the
    user input.

    :param      input_filename  String  Absolute path to modify           :
    :param      add_suffix      String  String to add after the filename  :
    :param      add_prefix      String  String to add before the filename :
    :return     new_filename    String  The original filename modified    :
    """
    filename_parent = pathlib.PurePath(input_filename).parents[0]
    filename = pathlib.PurePath(input_filename).name
    if add_suffix and not add_prefix:
        new_filename = filename.rsplit('.', 1)[0] + add_suffix + filename.rsplit('.', 1)[1]              # Add suffix
    elif add_prefix and not add_suffix:
        new_filename = add_prefix + filename                                                             # Add prefix
    elif add_prefix and add_suffix:
        new_filename = add_prefix + filename.rsplit('.', 1)[0] + add_suffix + filename.rsplit('.', 1)[1] # Add both
    else:
        print('Input filename it is not modified. Please provide a prefix or suffix. ')
        new_filename = input_filename
    new_filename = filename_parent.joinpath('PoGo_Output', new_filename)
    print(new_filename)
    return new_filename

def gen_dir_structure(proj_dir_path):
    """
    Version : 1.0
    Name History : gen_dir_structure

    This function generates the directories structure tu run Proteogenome 3
    It receives a directory path as input. Inside this directory generates the 'Proteogenome_Project'. Therefore, inside
    this parent directory generates two subdirs for running PoGo

    :param      proj_dir_path       String      The path where create the Proteogenome project folder:
    :return:
    """
    # Create the Proteogenome_Project folder
    os.chdir(proj_dir_path)
    proj_dir_path = proj_dir_path.joinpath('Proteogenome_Project')
    os.mkdir(proj_dir_path)
    # Create the subfolder for the Proteogenome OUTPUT
    proteogenome_output_path = proj_dir_path.joinpath('Proteogenome_Output')
    os.mkdir(proteogenome_output_path)
    print(proteogenome_output_path)

    # Create the subfolder for the Proteogenome data structure OUTPUT
    proteogenome_data_structure = proj_dir_path.joinpath('Proteogenome_Data_Structure')
    os.mkdir(proteogenome_data_structure)

    # Create the subfolder for the PoGo INPUT files
    pogo_input_path = proj_dir_path.joinpath('PoGo_Input')
    os.mkdir(pogo_input_path)
    # Create the subfolder for the PoGo OUTPUT files
    pogo_output_path = proj_dir_path.joinpath('PoGo_Output')
    os.mkdir(pogo_output_path)
    return proj_dir_path, proteogenome_output_path, proteogenome_data_structure, pogo_input_path, pogo_output_path

def move_files(origin_dir, destination_dir, filenames_patterns=[]):
    """

    :param origin_dir:
    :param destination_dir:
    :param filenames_patterns:
    :return:
    """
    print(f'Moving these files from {origin_dir} to {origin_dir}\nFile list:\n')
    for file_name_pattern in filenames_patterns:
        for origin_abs_file_name in origin_dir.glob(file_name_pattern):
            file_name = origin_abs_file_name.name  # Extract only the filename with its extension from the absolute path
            dest_abs_file_name = destination_dir.joinpath(file_name)
            print(f'{origin_abs_file_name} --------> {dest_abs_file_name}')
            origin_abs_file_name.rename(dest_abs_file_name)

def Ensembl_API_symbol(accession, speces):
    """
    Version : 1.0

    The API entry point for the Ensemble database (adapted from the https://rest.ensembl.org/documentation/info/xref_external)

    The conversion will be performed for a single accession at once (only GET method was available)
    Quote Ensemble documentation for this entry point:
    "Looks up an external symbol and returns all Ensembl objects linked to it. This can be a display name for a
    gene/transcript/translation, a synonym or an externally linked reference. If a gene's transcript is linked to the
    supplied symbol the service will return both gene and transcript (it supports transient links)."
    """
    #import requests, sys
    decoded = 0
    server = "https://rest.ensembl.org"
    entry_point = f'/xrefs/symbol/{speces}/{accession}?'

    r = requests.get(server + entry_point, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    # print(accession)
    # print(decoded)
    # a=input()
    return decoded


def UniProt2Ensembl(accession_lst, speces, isoforms=False, reqs_per_sec=15):
    """
    Version : 1.0
    Convert a list of UniProt codes in Ensemble protein accession codes.
    The conversion could be performed to obtain only the canonical accession or including isoforms.
    :param  reqs_per_sec        Int     The limit of request per secons allowed for the Ensembl server:
    :return conversion_dict     Dict    The dictionary with k: UniProt accession
                                                            v: List of conversion in Ensembl
    """
    import time

    start = time.time()

    req_count = 0
    last_req = 0

    print(f'Protein to convert - {len(accession_lst)}')
    conversion_dict = {}
    for accession_n, accession in enumerate(accession_lst):
        # print(accession)
        print(f'{accession_n + 1} - {accession}                  ', end='\r')
        # check if we need to rate limit ourselves
        if req_count >= reqs_per_sec:
            delta = time.time() - last_req
            if delta < 1:
                time.sleep(1 - delta)
            last_req = time.time()
            req_count = 0

        conversion_lst = Ensembl_API_symbol(accession, speces)
        req_count += 1
        conversion_lst = [id_type['id'] for id_type in conversion_lst if
                          id_type['type'] == 'translation']  # Extract canonical and isoforms accessions

        # Only if you want the first Ensemble translation code of the list (should be the canonical form)
        if not isoforms and conversion_lst:
            conversion_dict[accession] = conversion_lst[0]
        else:
            conversion_dict[accession] = conversion_lst
    end = time.time()
    print(f'\n{end - start} sec')
    return conversion_dict


