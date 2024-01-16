import re

def FASTA_cpt_seq(list_rows):
    """
    Version: 1.0
    Name History: FASTA_cpt_seq

    This function compact the sequence after each FASTA header.
    It is possible that FASTA files downloaded from online sources could have
    the sequences splitted in multiple rows by '\n'.
    The purpose is to arrange each sequence in a unique string.

    :param  list_rows   List[Str]   List of rows of a FASTA file.:
    :return seq_compa   List[Str]   List of rows of a FASTA file with the sequences compacted in unique rows.:
    """
    seq_compa = []
    sequence_row = ''

    for row in list_rows:
        if (len(row) > 0) and (row[0] == '>'):
            if sequence_row != '':  # If the sequence line is empty but there is a FASTA header than this is the first header.
                seq_compa.append(
                    sequence_row)  # Otherwise, it is a sequence that belongs to the current header and then it can be appended.
            seq_compa.append(row)  # If the current line is a FASTA header then append to the list.
            sequence_row = ''  # This means that you are expecting for a new sequence into the next rows.
        else:
            sequence_row += row  # Append the current part of sequence to the whole sequence.
    return seq_compa


def check_GFF3_format(list_rows):
    """
    Version: 1.0
    Name History: check_format (from Proteogenome 2), check_GFF3_format

    This function controls if the first row in a datafile contains the expected data.
    The default is the GFF3 control where '##' chars represent comment lines and
    the first column in a .split() row contains the strain identifier.

    :param:
    :return:
    """
    format_control = True  # I suppose that the file has the expected format
    strain_ID = 'NA'

    for row in list_rows:
        # print(row[:2])                # Find the strain ID as future refernce
        if row[:2] != '##':             # The '##' represent comment lines in the GFF3 format
            strain_ID = row.split()[0]  # In a normal row the strain ID is in the first column
            break
    print(f'Row header - {strain_ID}')
    for row in list_rows:                # Restart from the beginning of the list

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

    print('++ GFF3 - FILE FORMAT CONTROL FINISHED ++')
    return format_control

def check_FASTA_format(list_rows):
    """
    Version 1.0
    Name History: check_FASTA_format

    This function check if a list of FASTA rows contains undesired characters.
    If yes, then returns the list of undesired characters found.
    Moreover, check if the sequences in the FASTA are represented in singe lines (NO MULTILINES SEQUENCES)

    :param      list_rows        :

    :return     unique_char_found:
    :return     compact_flag     :
    """
    compact_flag = True
    upper_alfa = {'A','B','C','D','E','F','G','H','I','J','K','L','M',
                  'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'}
    compact_FASTA = False
    fasta_header_flag = False
    FASTA_seq_flag = False
    FASTA_seq_lines_count = 0

    unique_char_found = []                          # Unique list of undesired characters found in the input list of rows
    undesired_chars_set = {'[', ']', '{', '}','\n'} # Undesired characters set, not compliant with PoGo
    format_control = True                           # I suppose that the file has the expected format
    for row in list_rows:
        row_set = set(row)   # Make a set of the current row for further comparison
        chars_found = []
        chars_found = list(set(row).intersection(undesired_chars_set))           # Find undesired chars
        if chars_found:
            print('Undesired characters found in the FASTA file.')
            update_unique_char_found =[]
            # Check if the character found are already found before
            update_unique_char_found = set(chars_found) - set(unique_char_found)
            # If NOT update the final list of undesired characters found
            if update_unique_char_found: unique_char_found.extend(list(update_unique_char_found))

        if not compact_flag:      # Check for NO multiline FASTA sequences
            if (row[0] == '>') and (FASTA_seq_flag == False) and (FASTA_seq_lines_count == 0): # FASTA header condition
                fasta_header_flag = True
                FASTA_seq_flag = False
                FASTA_seq_lines_count = 0
            elif row.issubset(upper_alfa):  # This is a FASTA sequence
                fasta_header_flag = False
                FASTA_seq_flag = True
                FASTA_seq_lines_count += 1

            if FASTA_seq_lines_count > 1:
                print('+++++ Multiline FASTA sequence detected. The FASTA sequences will be compacted +++++')
                compact_flag = True
            elif (fasta_header_flag == True) and (FASTA_seq_flag == False) and (FASTA_seq_lines_count == 0) and (row[0] == '>'):
                print('--- WARNING --- \nApparently one FASTA header has no sequence.\nPlease check your FASTA file')
                break

    return unique_char_found, compact_flag


def locus_tag_substitution(FASTA_lst):
        """
        Version: 1.0
        Name History: locus_tag_substitution

        This function substitutes the tag 'locus_tag='' with 'gene:gene-' and 'transcript:rna-'
        into the FASTA rows.
        This function is specific to the single operation of fitting the FASTA header into
        a PoGo readable format.

        :param      FASTA_lst       List    List of FASTA rows including also FASTA sequences.                    :
        :return     FASTA_lst4PoGo  List    List of FASTA rows with the rectified headers and the FASTA sequences.:
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

    This function receives a list of rows. Therefore, cleans every single row replacing the
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


    :param  list_rows           List    List of rows to clean:
    :param  target_sub_str      List    List of [tuples] of two elements. Replacement through .replace() command.
                                        Tuple position reference:
                                        [0] target substring that must be substitute.
                                        [1] replacement string used in the substitution.                             :
    :param  target_patterns     List    List of [tuples] of two elements. Replacement through
                                        .sub() command in a regex search.
                                        Tuple position reference
                                        [0] target pattern that must be substitute.
                                        [1] replacement string used in the substitution.                             :
    :param  open_patterns       List    List of [tuples] of two elements. Replacement through .sub() command in
                                        a regex search.
                                        Tuple position reference:
                                        [0] target pattern that must be substitute.
                                        [1] replacement pattern used to find the substring that will substitute the
                                        value of the pattern [0].                                                    :

    :return list_rows           List    The original list of rows cleaned from the substrings provided as input.     :
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

def FASTA_cleaning(FASTA_path, out_file_path, chr_to_remove='[\[,\],{,}]', remove_original_file=False):
    """
    Version: 1.0

    Name History: FASTA_cleaning

    This function cleans the FASTA file from undesired characters.
    In particular remove the '\n' character from the FASTA sequences. This is necessary to further sequencing process.

    :param      FASTA_path     :
    :param      out_file_path  :
    :param      chr_to_remove  :
    :return                    :
    """
    with open(out_file_path, 'w') as outfile:
        for record in SeqIO.parse(FASTA_path, 'fasta'):
            # Clean the header
            header = record.description
            header = re.sub(chr_to_remove, '', header)
            record.description = header

            # Clean the sequence
            sequence = str(record.seq)
            if '\n' in sequence:
                sequence = sequence.replace('\n', '')

            # Update the record's sequence
            record.seq = sequence

            # Write the modified record to the output file
            SeqIO.write(record, outfile, 'fasta')
