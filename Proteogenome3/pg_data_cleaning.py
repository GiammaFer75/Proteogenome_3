def FASTA_cpt_seq(list_rows):
    """
    Version: 1.0

    Name History: FASTA_cpt_seq

    This function compact the sequence after each FASTA header.
    It is possible that FASTA files downloaded from online sources could have
    the sequences splitted in multiple rows by '\n'.
    The purpose is to arrange each sequence in a unique string.

    INPUT : list_rows   List[Str]   List of rows of a FASTA file.
    OUTPUT: seq_compa   List[Str]   List of rows of a FASTA file with the sequences
                                    compacted in unique rows.
    """
    seq_compa = []
    sequence_row = ''

    for row in list_rows:
        if row[0] == '>':
            if sequence_row != '':  # If the sequence line is empty but there is a FASTA header than this is the first header.
                seq_compa.append(
                    sequence_row)  # Otherwise, it is a sequence that belongs to the current header and then it can be appended.
            seq_compa.append(row)  # If the current line is a FASTA header then append to the list.
            sequence_row = ''  # This means that you are expecting for a new sequence into the next rows.
        else:
            sequence_row += row  # Append the current part of sequence to the whole sequence.
    return seq_compa


def check_format(list_rows):
    """
    Version: 1.0

    Name History: check_format

    This function controls if the first row in a datafile contains the expected data.
    The default is the GFF3 control where '##' chars represent comment lines and
    the first column in a .split() row contains the strain identifier.

    INPUT :
    OUTPUT:
    """
    format_control = True  # I suppose that the file has the expected format
    strain_ID = 'NA'

    for row in list_rows:
        # print(row[:2])             # Find the strain ID as future refernce
        if row[:2] != '##':  # The '##' represent comment lines in the GFF3 format
            strain_ID = row.split()[0]  # In a normal row the strain ID is in the first column
            break
    print(f'Row header - {strain_ID}')
    for row in list_rows:  # Restart from the beginning of the list

        if row[:2] != '##':
            if row:
                current_strain_ID = row.split()[0]
                if current_strain_ID != strain_ID:
                    print('INVALID ROW')
                    print(row)
                    format_control = False
                    break
    if format_control:
        print('FILE FORMAT - OK')
    else:
        print('XXXXXX  WARNING - WRONG FILE FORMAT  XXXXXX')

    print('++ FILE FORMAT CONTROL FINISHED ++')



def locus_tag_substitution(FASTA_lst):
        """
        Version: 1.0

        Name History: locus_tag_substitution

        This function substitutes the tag 'locus_tag='' with 'gene:gene-' and 'transcript:rna-'
        into the FASTA rows.
        This function is specific to the single operation of fitting the FASTA header into
        a PoGo readable format.

        INPUT :     FASTA_lst       List    List of FASTA rows including also FASTA sequences.
        OUTPUT:     FASTA_lst4PoGo  List    List of FASTA rows with the rectified headers and
                                            the FASTA sequences.
        """
        FASTA_lst4PoGo=[]

        locus_tag_patt=re.compile(r'.*?locus_tag=(.*?)\s')
        for row in FASTA_lst:
            if row[0] =='>':
                target_tag=locus_tag_patt.search(row).group(1)
                row=row.replace('locus_tag='+target_tag,
                                'gene:'+target_tag+' transcript:'+target_tag)
                                #gene:gene-   The first suggestion was to add gene:gene- but later was to remove them.

            FASTA_lst4PoGo.append(row)
                #print(row)
                #a=input()
        return FASTA_lst4PoGo


def rectify_rows(list_rows, target_sub_str=[], target_patterns=[],
                 open_patterns=[], rows_not_modif=False):
    """
    Version: 1.0

    Name History: rectify_rows

    This function receives a list of rows. Therefore cleans every single row replacing the
    target substrings with the proper substitution.
    Target substrings can be substitute in two ways:

            - target_sub_str:   are substrings known before the substituting operation
                                and are applied through the command ".replace()".
                                These patterns must be passed through the list
                                "target_sub_string"

            - target_patterns:  substring that must be replaced in the file rows occur in
                                the same pattern. This pattern must be passed as a regerx
                                Python patterns enclosed in ' '. This option search for the
                                pattern in the row. If it is found, the pattern is replaced
                                by the substring passed in the tuple with the pattern itself.

            - open patterns:    In this case also the substring that will replace the target
                                pattern is not specified explicitely as second argument in the
                                tuple. Instead, it must be extracted by a generic pattern as well.


    INPUT : list_rows           [List]  List of rows to clean
            target_sub_str      [List]  List of [tuples] of two elements. Replacement through
                                        .replace() command.
                                        Tuple position reference:
                                        [0] target substring that must be substitute.
                                        [1] replacement string used in the substitution.
            target_patterns     [List]  List of [tuples] of two elements. Replacement through
                                        .sub() command in a regex search.
                                        Tuple position reference:
                                        [0] target pattern that must be substitute.
                                        [1] replacement string used in the substitution.
            open_patterns       [List]  List of [tuples] of two elements. Replacement through
                                        .sub() command in a regex search.
                                        Tuple position reference:
                                        [0] target pattern that must be substitute.
                                        [1] replacement pattern used to find the substring
                                            that will substitute the value of the pattern [0].

    OUTPUT: list_rows           [List]  The original list of rows cleaned from the
                                        substrings provided as input.
    """
    rows_not_modif_lst = []
    rows_not_modif_flag = False

    for ind, row in enumerate(list_rows):
        # print('*'*30,' NEW ROW ','*'*30)
        # print('original row - ', row)
        if target_sub_str:
            for substitution in target_sub_str:
                row = row.replace(substitution[0], substitution[1])
                # print(row)

        if target_patterns:
            for substitution in target_patterns:
                #        |target pattern|    |replacement sub|
                row = re.sub(r'' + substitution[0] + '', substitution[1], row)
                # print(row)

        if open_patterns != []:
            for substitution_tuple in open_patterns:
                try:
                    to_be_replaced_pat = re.compile(r'' + substitution_tuple[0] + '')
                    replace_pat = re.compile(r'' + substitution_tuple[1] + '')

                    to_be_replaced_value = to_be_replaced_pat.match(row).group(1)
                    replace_value = replace_pat.match(row).group(1)

                    row = re.sub(to_be_replaced_value, replace_value, row)
                    # print(row)
                except:
                    rows_not_modif_lst.append(row)
                    rows_not_modif_flag = True
                    # print(row)
                    # print('+'*100)
                    # a=input()

        # print(row)
        # a=input()
        list_rows[ind] = row
        # print('\nlist_rows\n',list_rows)

    if rows_not_modif_flag == True:
        print('§§§§§§§§§§§ Some rows are not affected by the substitution §§§§§§§§§§§')

    if (rows_not_modif == False):
        return list_rows
    else:
        return list_rows, rows_not_modif_lst
