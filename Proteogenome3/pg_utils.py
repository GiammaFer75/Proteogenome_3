import os
import shutil

def print_lst(input_list, limit=10, en_sep=True, sep_type='-'):
    """
    Version : 2.0
    Name History : print_lst (from Proteogenome 1.0)

    Print the list content limitely to the element indicated by the parameter 'limit'.

    :param      input_list  List    Items to print.                                                  :
    :param      limit       Int     Number of items to print.                                        :
    :param      en_sept     Bool    Enable the separation of each element of the list with a string. :
    :param      sep_type    Char    The character that will make up the separator string.            :

    :return     Print the data structure content:
    """

    for ind, i in enumerate(input_list):
        if ind < limit:
            print(i)
            if en_sep: print(sep_type * 100)


def create_path(input_filename, add_suffix='', add_prefix=''):
    """
    Version : 1.0
    Name History : create_path

    This function take the absolute path of a file and add a prefix or a suffix to the original filename, based on the
    user input.

    :param      input_filename:
    :param      add_suffix:
    :param      add_prefix:
    :return     new_filename:
    """
    cwd = os.getcwd()
    filename = os.path.basename(input_filename)
    if add_suffix and not add_prefix:
        new_filename = filename.rsplit('.', 1)[0] + add_suffix + filename.rsplit('.', 1)[1]              # Add suffix
    elif add_prefix and not add_suffix:
        new_filename = add_prefix + filename                                                             # Add prefix
    elif add_prefix and add_suffix:
        new_filename = add_prefix + filename.rsplit('.', 1)[0] + add_suffix + filename.rsplit('.', 1)[1] # Add both
    else:
        print('Input filename it is not modified. Please provide a prefix or suffix. ')
        new_filename = input_filename
    new_filename = os.path.abspath(new_filename)
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
    os.chdir(proj_dir_path)
    proj_dir_path = proj_dir_path.joinpath('Proteogenome_Project')
    os.mkdir(proj_dir_path)  # Create the Proteogenome_Project folder
    pogo_input_path = proj_dir_path.joinpath('PoGo_Input')
    os.mkdir(pogo_input_path)
    pogo_output_path = proj_dir_path.joinpath('PoGo_Output')
    os.mkdir(pogo_output_path)
    return proj_dir_path, pogo_input_path, pogo_output_path

def move_files(origin_dir, destination_dir, filenames_patterns=[]):
    print(f'Moving these files from {origin_dir} to {origin_dir}\nFile list:\n')
    for file_name_pattern in filenames_patterns:
        for origin_abs_file_name in origin_dir.glob(file_name_pattern):
            file_name = origin_abs_file_name.name  # Extract only the filename with its extension from the absolute path
            dest_abs_file_name = destination_dir.joinpath(file_name)
            print(f'{origin_abs_file_name} --------> {dest_abs_file_name}')
            origin_abs_file_name.rename(dest_abs_file_name)

