from Proteogenome3 import pg_indexes as pg_i


def make_sep_file(out_file_name, input_array, sep=''):
    """
    Version: 1.0
    Name History: make_tab_file - make_sep_file

    This function creates a file from an input list or an np.array made of strings.
    Each string has separated fields.

    The function verify the data type of the input array (list or np.array). Three cases are managed:
    1 - List of strings
    2 - List of lists
    3 - np.array

    The function retrieves the number of columns looking at the number of columns stored
    in the first row of the input list/array.
    The column's content will be written into the file separated by the value of the sep variable.

    :param      out_file_name   String      The file path where the new file will be created
                input_array     np.array    Array of data to write in the output file
                                List        List of data to write in the output file

    :return The file with the input array content
    """

    out_file_hand = open(out_file_name, 'w')
    if type(input_array) == list and type(input_array[0]) == list:  # LIST of LISTS
        number_of_columns = len(input_array[0])
        for row in input_array:
            out_row = ''
            for col_ind, col_data in enumerate(row):
                out_row += col_data
                if (col_ind < number_of_columns): out_row += sep
            out_row += '\n'
            out_file_hand.write(out_row)
        print('1 - LIST of LISTS')

    if type(input_array) == list and type(input_array[0]) == str:   # LIST of STRINGS
        number_of_columns = 1
        for row in input_array:
            out_file_hand.write(row + '\n')
        print('2 - LIST of STRINGS')

    if 'numpy.ndarray' in str(type(input_array)):                   # NUMPY np.array
        number_of_columns = input_array.shape[1]
        for row in input_array:
            out_row = ''
            for col_ind, col_data in enumerate(row):
                out_row += col_data
                if (col_ind < number_of_columns): out_row += sep
            out_row += '\n'
            out_file_hand.write(out_row)
        print('3 - NUMPY np.array')

    out_file_hand.close()

def gen_protein_track(prot_pep_index, prot_CDS_index, prot_list=[], strand='NC_006273.2', bed_fn='test1.bed'):
    """
    Version: 1.0
    Name History: protein_set_bed - protein_track (from Proteogenome 3) - gen_protein_track

    Generates a .bed track with the map showing the protein location. The map is created as an heat-map based on the
    protein expression quantification.

    Example of multi exon protein in input:
      UL37 - [['52573', '53060', '-'], ['51302', '51344', '-'], ['50262', '51197', '-']]

    :param      prot_pep_index    :
    :param      prot_CDS_index    :
    :param      prot_list         List      Specifying the list of proteins, it is possible to visualise only a subset
                                            of the entire set of proteins founded in the sample :
    :return                       :
    """

    prot_PSMint_index = pg_i.protein_PSM_int_index(prot_pep_index)  # Generates the color code for protein intensities

    if prot_list == []:
        prot_list = prot_pep_index.keys()

    print(f'Start processing {len(prot_list)} proteins')
    prot_track_fh = open(bed_fn, 'w')

    chr_name = ''
    if strand == 'NC_006273.2':
        chr_name = 'NC_006273.2'  # Find the chromosome name
    else:
        print('Option under construction. Check how to manage chromosome name in Human annotations\FASTA')

    BED_rows = []
    prot_row = ''

    for protein in prot_list:
        prot_row = ''
        CDS_block = prot_CDS_index[protein]  # Extract the block of CDS coordinates

        Strand = CDS_block[0][2]  # The strand is in the third position of the first CDS features record

        if Strand == '+':  # If strand positive
            chromStart = CDS_block[0][0]  # Chromosome Start is in the first position of the first record
            chromEnd = CDS_block[-1][1]  # Chromosome End is in the second position of the last record
        else:
            chromStart = CDS_block[-1][1]  # In the negative case the bonduaries coordinates
            chromEnd = CDS_block[0][0]  # inverted their orders in the respective records

        print(prot_CDS_index[protein])
        print(f'{chromStart} - {chromEnd}')

        # Score='1'
        tickStart = chromStart
        tickEnd = chromEnd
        blockCount = str(len(CDS_block))

        itemRGB = prot_PSMint_index[protein][2]  # Fetch the protein intensity into the proper index
        Score = str(prot_PSMint_index[protein][1])  # Set the Score column in the bed protein file
        # to the total intensity of the current protein.

        blockSizes_str = ''
        blockStarts_str = ''
        for CDS in CDS_block:
            print('CDSb', CDS_block)
            if Strand == '+':
                CDS_start = int(CDS[0])
                CDS_end = int(CDS[1])
            else:
                CDS_start = int(CDS[1])
                CDS_end = int(CDS[0])

            blockSizes_str += str(abs(CDS_end - CDS_start)) + ','
            blockStarts_str += str(abs(int(chromStart) - CDS_start)) + ','

        prot_row += (chr_name + '\t' + chromStart + '\t' + chromEnd + '\t' + protein + '\t' + Score + '\t' +
                     Strand + '\t' + tickStart + '\t' + tickEnd + '\t' + itemRGB + '\t' + blockCount + '\t' +
                     blockSizes_str + '\t' + blockStarts_str + '\n')

        print(prot_row)
        print('-' * 100)

        BED_rows.append(prot_row)

    for row in BED_rows:
        prot_track_fh.write(row)
    prot_track_fh.close()

