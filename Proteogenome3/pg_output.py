from Proteogenome3 import pg_indexes as pg_i

import json


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

    print(f'Saving procedure of variable type - {type(input_array)}')

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

def save_dict_list(input_data, out_file_name, input_data_type='d'):
    """
    This function save the content of a list or a dictionary
    :param input_data:
    :param out_file_name:
    :param input_data_type:
    :return: A file containing the input data structure
    """
    out_file_hand = open(out_file_name, 'w')
    if input_data_type == 'd':
        for k,v in input_data.items():        # Dictionary
            out_file_hand.write(f'{k} - {v}\n')
    elif input_data_type == 'l':              # List
        for row in input_data:
            out_file_hand.write(f'{row}\n')
    out_file_hand.close()

def save_json_file(input_file_path, json_object):
    # Save the conversion table from uniprot to Ensembl in a JSON file
    f = open(input_file_path, 'w')
    json.dump(json_object, f)
    f.close()
def gen_protein_track(prot_pep_index, prot_CDS_index, prot_list=[], strand='', bed_fn='test1.bed'): # strand for HCMV NC_006273.2
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

    proteins_not_found = []
    prot_PSMint_index = pg_i.protein_PSM_int_index(prot_pep_index)  # Generates the color code for protein intensities

    if prot_list == []:
        prot_list = prot_pep_index.keys()

    print(f'Start processing {len(prot_list)} proteins')
    prot_track_fh = open(bed_fn, 'w')

    strand_position = 0
    start_p, start_n, end_p, end_n = 0, 0, 0, 0
    organism = ''

    # Define the array parser based on the organism
    if 'NC_006273' in strand:                              ## HCMV herpes virus 5 ##
        organism = 'virus'
        chr_name = 'NC_006273.2'  # The chromosome name (only one in viruses) is always the same in the HCMV
        start_p, end_p = 0, 1
        start_n, end_n = 1, 0
        strand_position = 2

    elif 'GRCh38' in strand:                                ## HOMO SAPIENS GRCh38 ##
        # print('Option under construction. Check how to manage chromosome name in Human annotations\FASTA')
        organism = 'homo'
        start_p, end_p, start_n, end_n = 1, 2, 1, 2
        # start_p, end_p = 1, 2
        # start_n, end_n = 2, 1
        strand_position = 3

    BED_rows = []
    prot_row = ''

    for protein in prot_list:
        prot_row = ''
        try:
            CDS_block = prot_CDS_index[protein]  # Extract the block of CDS coordinates
        except:                                  # If the protein code is not in the index exit the loop
            proteins_not_found.append(protein)
            continue

        # Strand = CDS_block[0][2]  # The strand is in the third position of the first CDS features record
        #
        # if Strand == '+':  # If strand positive
        #     chromStart = CDS_block[0][0]  # Chromosome Start is in the first position of the first record
        #     chromEnd = CDS_block[-1][1]  # Chromosome End is in the second position of the last record
        # else:
        #     chromStart = CDS_block[-1][1]  # In the negative case the bonduaries coordinates
        #     chromEnd = CDS_block[0][0]  # inverted their orders in the respective records

        Strand = CDS_block[0][strand_position]  #

        if Strand == '+':  # If strand positive
            chromStart = CDS_block[0][start_p]  # Chromosome Start is in the first position of the first record
            chromEnd = CDS_block[-1][end_p]  # Chromosome End is in the second position of the last record
        else:
            chromStart = CDS_block[-1][start_n]  # In the negative case the bonduaries coordinates
            chromEnd = CDS_block[0][end_n]  # inverted their orders in the respective records

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
        print('CDSb', CDS_block)
        for CDS in CDS_block:
            if organism == 'homo': chr_name = CDS[0] # For human genome the chromosome name is in the first position

            if Strand == '+':
                CDS_start = CDS[start_p]
                CDS_end = CDS[end_p]
            else:
                CDS_start = CDS[start_n]
                CDS_end = CDS[end_n]

            blockSizes_str += str(int(abs(CDS_end - CDS_start))) + ','
            blockStarts_str += str(int(abs(chromStart - CDS_start))) + ','
        if organism == 'homo':
            blockStarts_str = blockStarts_str.split(',')
            blockStarts_str = [x for x in blockStarts_str if x != '']
            blockStarts_str.reverse()
            blockStarts_str = ','.join(blockStarts_str)

            blockSizes_str = blockSizes_str.split(',')
            blockSizes_str = [x for x in blockSizes_str if x != '']
            blockSizes_str.reverse()
            blockSizes_str = ','.join(blockSizes_str)

        if type(CDS_start) == float:        # For the genomes manipulated with pd.DataFrame the coordinates are floats NOT str
            tickStart = str(int(tickStart))
            tickEnd = str(int(tickEnd))
            chromStart = str(int(chromStart))
            chromEnd = str(int(chromEnd))
        #
        # print(f"""
        # chr_name        - {type(chr_name)}\n
        # chromStart      - {type(chromStart)}\n
        # chromEnd        - {type(chromEnd)}\n
        # protein         - {type(protein)}\n
        # Score           - {type(Score)}\n
        # Strand          - {type(Strand)}\n
        # tickStart       - {type(tickStart)}\n
        # tickEnd         - {type(tickEnd)}\n
        # itemRGB         - {type(itemRGB)}\n
        # blockCount      - {type(blockCount)}\n
        # blockSizes_str  - {type(blockSizes_str)}\n
        # blockStarts_str - {type(blockStarts_str)}\n""")
        prot_row += (chr_name + '\t' + chromStart + '\t' + chromEnd + '\t' + protein + '\t' + Score + '\t' +
                     Strand + '\t' + tickStart + '\t' + tickEnd + '\t' + itemRGB + '\t' + blockCount + '\t' +
                     blockSizes_str + '\t' + blockStarts_str + '\n')

        print(prot_row)
        print('-' * 100)

        BED_rows.append(prot_row)

    for row in BED_rows:
        prot_track_fh.write(row)
    prot_track_fh.close()

    return proteins_not_found

