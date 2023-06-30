def print_lst(input_list, limit=0, en_sep=True, sep_type='-'):
    """
    Version : 2.0

    Name History : print_lst (from Proteogenome 1.0)

    Print the list content limitely to the element indicated by the parameter 'limit'.

    INPUT:  input_list  List    Items to print.
            limit       Int     Number of items to print.
            en_sept     Bool    Enable the separation of each element of the list with a string.
            sep_type    Char    The character that will make up the separator string.


    OUTPUT :
    """

    for ind, i in enumerate(input_list):
        if ind < limit:
            print(i)
            if en_sep: print(sep_type * 100)
