import os
def print_lst(input_list, limit=10, en_sep=True, sep_type='-'):
    """
    Version : 2.0

    Name History : print_lst (from Proteogenome 1.0)

    Print the list content limitely to the element indicated by the parameter 'limit'.

    :param  input_list  List    Items to print.
            limit       Int     Number of items to print.
            en_sept     Bool    Enable the separation of each element of the list with a string.
            sep_type    Char    The character that will make up the separator string.
    :
    :return
    :
    """

    for ind, i in enumerate(input_list):
        if ind < limit:
            print(i)
            if en_sep: print(sep_type * 100)


def create_path(input_filename, add_suffix, add_prefix):
    """
    Version : 1.0
    Name History : create_path

    This function take the absolute path of a file and add a prefix or a suffix to the original filename, based on the
    user input.

    :param input_filename:
    :param add_suffix:
    :param add_prefix:
    :return:
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