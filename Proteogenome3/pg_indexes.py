from Proteogenome3 import pg_data_processing as pg_dp
from Proteogenome3 import pg_input

import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import re

def initialise_indexes(prot_annots_path, input_table, annot_format='gff3'):
    """
    Version: 1.1
    Name History: initialise_indexes

    This function create the index for
        - CDS     ---> annotation - in the genome
        - Protein ---> CDS        - that cotribute to the protein production
        - Protein ---> Peptides   - that confirmed the protein existence in the proteomics data
        - Peptide ---> Proteins   - where the peptide sequence has been found

    """
    print("""
           **********************
           INDEXES INITIALISATION
           **********************
           """)

    genome_annot_lst = pg_input.file_to_lst(prot_annots_path)                             # Upload the genome annotation file into a list
    CDS_matrix = pg_dp.CDS_annot_matrix(genome_annot_lst)                                 # Matrix with CDS annotations info only
    prot_CDS_index = gen_protein_CDS_index(CDS_matrix, annot_format=annot_format)         # Dictionary of proteins and their group of CDS
    prot_pep_index = protein_peptide_index(input_table, )
    pep_prot_index = peptide_protein_index(input_table)

    # try:
    #     prot_pep_index = protein_peptide_index(input_table)
    #     pep_prot_index = peptide_protein_index(input_table)
    # except:
    #     print('Proteomics data not found')

    return CDS_matrix, prot_CDS_index, prot_pep_index, pep_prot_index

def protein_peptide_index(pep_input_table):
    """
    Version: 1.0

    History Name: protein_peptide_index

    This function parse the input table and group all the peptides by protein ID

    INPUT : pep_input_table         np.array

    OUTPUT: prot_pep_index     [Dict]  The protein index with the CDS coordinates.
                                        (Key)   [Str]               Protein ID
                                        (Value) [List][List][Str]   peptides =
                                                                    pepsequence, intensisty
    """
    prot_pep_index = {}  # Initialise dictionary for protein ---> peptide index
    prot_ID_array = np.unique(pep_input_table[:, 0])  # Extract the protein IDs
    print('Protein ID array\n----------------\n', prot_ID_array)
    print('\n-------------------------\nGROUP PEPTIDES BY PROTEIN\n')

    for protein in prot_ID_array:
        peptides_block = pg_dp.groupby_protein(pep_input_table,
                                              protein)  # Group the peptides that belong to the current protein.
        prot_pep_index[protein] = peptides_block
        print(protein, '--\n', peptides_block)
    return prot_pep_index


def peptide_protein_index(pep_input_table):
        """
        Version: 1.0

        History Name: peptide_protein_index

        This function parse the input table and group all the protein ID by peptide sequence.

        INPUT : pep_input_table         np.array

        OUTPUT: pep_prot_index     [Dict]  The peptide sequence.
                                                (Key)   [Str]       peptide sequence
                                                (Value) [List][Str] Protein IDs where the peptide come from.

        """
        pep_prot_index = {}  # Initialise dictionary for peptide ---> protein index
        peptide_seq_array = np.unique(pep_input_table[:, 1])  # Extract the protein IDs
        print('\nUnique peptide sequences array\n-----------------------\n', peptide_seq_array)
        print('\nNumber of unique peptide sequences - ', len(peptide_seq_array),
              '\n----------------------------------------\n')

        for peptide in peptide_seq_array:
            protein_block = pg_dp.groupby_peptide(pep_input_table,
                                                 peptide)  # Group the peptides that belong to the current protein.
            pep_prot_index[peptide] = protein_block
        return pep_prot_index


def gen_protein_CDS_index(annotations, annot_format='gff3'):
    """
    Version: 3.0

    Name History: prot_index, protein_CDS_index, gen_protein_CDS_index

    Generate a dictionary with the protein ID ('gene' field in GFF3) as key and
    its CDS list as value.

        Example:
         {'RL1'   : [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
          'RL5A'  : [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
          'RNA2.7': [ [chr,start,end,strand], [chr,start,end,strand], ....... ]
          .....
          }

    Two alternative patterns could be used to parse the protein identifier
    into the annotation file:
        - GFF3 ----> .*?;gene=(.*?);
        - GTF  ----> .*?;\sgene\s\"(.*?)\";


    INPUT : annotations     [List]  The annotation file stored in a list
            annot_format    [Str]   Tag that indicates which kind of format the annotation file is.

    OUTPUT: prot_CDS_index     [Dict]  The protein index with the CDS coordinates.
                                    (Key)   [Str]               Protein ID
                                    (Value) [List][List][Str]   CDS coordinates =
                                                                chr,start,end,strand
    """
    prot_CDS_index = {}

    # **************** Parsing annotations in GFF3 format
    if annot_format == 'gff3':  # specific patterns for the
        # gene_pat = re.compile(r'.*?;gene=(.*?);')  # annotation in GFF3 format
        gene_pat = re.compile(r'.*?;protein_id=(.*?);')  # annotation in GFF3 format


    # **************** Parsing annotations in GTF format
    elif annot_format == 'gtf':  # specific patterns for the
        # gene_pat = re.compile(r'.*?;\sgene\s\"(.*?)\";')  # annotation in GTF format
        gene_pat = re.compile(r'.*?;\sprotein_id\s\"(.*?)\";')  # annotation in GTF format

    for row in annotations:

        strand = row[6]

        if strand == '+':
            coord_1 = row[3]
            coord_2 = row[4]
        else:  # If the coding region is in reverse starnd
            coord_1 = row[4]
            coord_2 = row[3]

        print(row[-1])
        protein_ID = gene_pat.match(row[-1]).group(1)  # Find the current protein identifier.
        # print(row)
        # print(protein_ID)
        # print('-----------------------')
        # a=input()

        CDS_feat = [coord_1, coord_2, strand]

        if (protein_ID in prot_CDS_index):
            # Append the features of the current CDS on the list of CDS that belongs to the current protein.
            prot_CDS_index[protein_ID].append(CDS_feat)

        else:
            prot_CDS_index[protein_ID] = [CDS_feat]

    return prot_CDS_index


def protein_PSM_int_index(prot_pep_index,
                          color_gradient=['black', 'blue', 'cyan', 'green', 'greenyellow', 'yellow', 'orange', 'red']):
    # def protein_PSM_int_index(color_gradient=['blue','cyan','lime','yellow']):

    """
    Version: 1.0
    Name History: protein_PSM_int_index

    This function creates the protein index for the PSM and intensities.
    Moreover, it generates the RGB code for each protein intensity.
    The RGB codes will be used for the creation of the protein map.

    :param      prot_pep_index    :
    :param      prot_CDS_index    :
    :return     prot_PSMint_index :
    """
    import math
    def generate_color_gradient(color_lst, reverse_gradient=False):
        """
        Version: 1.0

        Name History: generate_color_gradient

        This function generates a list of tuples of RGB codes converted as three floating numbers.
        For each color in the input color_lst will be generated a colour fader block that goes from the previous colour
        in the list toward the next.
        All these colour fader blocks will be combined all together.
        The output is a list of RGB codes that represent each colour point in the color_lst evenly faded.

        ** MAIN FUNCTION ** ---> protein_PSM_int_index

        INPUT  :  color_lst     List    List of strings with the colour names that will be compose the colour gradient
                                        Example: ['black','gray','blue','green','yellow','red']
        OUTPUT :
        """

        def colorFader(c1, c2, mix=0):  # fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
            """
            Version: 1.0

            Name History: colorFader

            ** MAIN FUNCTION ** ---> generate_color_gradient

            This function generates an array that contains an RGB gradient from color1 to the color2.
            The gradient is expressed in RGB codes.
            INPUT :
            OUTPUT:
            """

            c1 = np.array(colors.to_rgb(c1))
            c2 = np.array(colors.to_rgb(c2))
            return colors.to_rgb((1 - mix) * c1 + mix * c2)

        def rgb_block(color1, color2, n_data_point):
            col = []
            for data_point in range(n_data_point):
                col.append(colorFader(color1, color2, data_point / n))
            return col

        n = 30
        color_blocks_lst = []
        for ind, color in enumerate(color_lst):  # For each color in the C list generate a gradient
            color_blocks_lst.append([color, color_lst[ind + 1]])
            if ind == (len(color_lst) - 2): break
        gradient = []
        for block in color_blocks_lst:
            gradient += rgb_block(block[0], block[1],
                                  n)  # Combine together all the color fader gradient for each color in the color_tuples
        if reverse_gradient: gradient.reverse()

        return gradient

    def exprlev_resc_RGB(RGB_scale):
        """
        Version: 2.0

        Name History: exprlev_resc_RGB

        For each protein, this function find the RGB codes for each intensity value
        in the range of the intensities.
        Returns two vctors one for the protein IDs and une for the RGB codes that represents
        the relative intensities.
        The two vectors are ordered by RGB codes from the brightest to the darkest color.

        INPUT : RGB_scale   [List]  List of tuples with the RGB codes converted in real numbers
                                    and ordered in ascending order.

        OUTPUT: proteins    [List]  Collection of protein IDs where the index of each ID is the index of
                                    its RGB code in RGB_vec.
                RGB_vec     [List]  RGB codes in tuples.
        """
        PSM_inten_values = prot_PSMint_index.values()
        inten_values = [int(x[1]) for x in PSM_inten_values]
        print('\nExtraction of protein intensities\n', inten_values)

        rescaled_values = []  # Intensities are converted to the index for the corresponding color code

        max_int_vec = [0]  # Intensities sorted in descending order
        proteins = []  # Protein IDs
        int_scaled = []  # Indexes referring to the RGB codes vctor
        RGB_vec = []  # RGB codes sorted in ascending order

        insert_ind = 0
        print('\n')
        print('Number of proteins                             - ',
              len(prot_pep_index.keys()))  # Number of proteins
        print('Number of intensities to convert in RGB tuples - ',
              len(prot_pep_index.keys()))  # Number of intensities
        print('\n')

        # Print the proteins with their intensities and RGB codes sorted descending
        print(
            'Protein Number\t-\tProtein Code\t-\tProtein intensity\t-\tPosition of RGB value\t-\tRGB Tuple')  # Must be the same
        print('--------------\t \t------------\t \t-----------------\t \t---------------------\t \t---------')

        # ************************************************************************* #
        # Based on the level of intensity, find the relative index in the RGB_scale #
        # The resulting RGB code will represent the relative protein intensity      #
        # translated in the proper color code.                                      #

        old_max = max(inten_values)
        old_min = min(inten_values)
        new_max = len(RGB_scale)
        new_min = 0

        for prot, PSM_intensity in prot_PSMint_index.items():
            old_value = int(PSM_intensity[1])

            NewValue = int((((old_value - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min)
            if NewValue == 0: NewValue = 1  # 0 values for the less intense proteins use to assign the higer intensity RGB code
            # For this reason, force to 1
            grather = False
            for ind, val in enumerate(max_int_vec):  # Find the position of the current intensity
                if old_value > val:  # in the intensities vector.
                    insert_ind = ind
                    grather = True
                    break

            if not grather:  # If the current intensity is the SMALLEST intensty in the vector:
                proteins.append(prot)  # - put this intensity in the last position of the vector.
                max_int_vec.append(old_value)  # - update all the correnspondent vectors.
                int_scaled.append(NewValue)
                RGB_vec.append(RGB_tup[NewValue - 1])
            else:
                proteins.insert(insert_ind,
                                prot)  # If the current intensity is AT LEAST greater than the smallest intensity in the vector:
                max_int_vec.insert(insert_ind,
                                   old_value)  # - put this intensity in the right position in the vector.
                int_scaled.insert(insert_ind,
                                  NewValue)  # - update all the correnspondent vectors in the same position.
                RGB_vec.insert(insert_ind, RGB_tup[NewValue - 1])
            rescaled_values.append(NewValue)

        # ************************************************************************* #

        # Print the content of the vectors used for the RGB assignment
        prot_number = 0
        for i in range(0, len(max_int_vec[:-1])):
            RGB = RGB_tup[int_scaled[i] - 1]
            print(
                f'    {prot_number}\t\t-\t     {proteins[i]}\t-\t      {max_int_vec[i]}       \t-\t         {int_scaled[i] - 1}       \t-\t{RGB}   ')
            # print(f'\t\t\t\t\t\t\t\t\t{RGB_vec[i]}')
            prot_number += 1
        return proteins, RGB_vec

        #  ------- MAIN --------
        #  protein_PSM_int_index

    intensities = {}

    max_intensity = 0
    min_intensity = 0

    prot_PSMint_index = {}
    for protein, pep_array in prot_pep_index.items():
        PSM_sum = 0
        inten_sum = 0
        for pep_row in pep_array:
            PSM_sum += int(pep_row[3])
            inten_sum += int(pep_row[4])
        prot_PSMint_index[protein] = [PSM_sum, inten_sum]

        if max_intensity < inten_sum: max_intensity = inten_sum
        if min_intensity > inten_sum: min_intensity = inten_sum

    #print('Protein-Peptide Index - Protein-CDS Index')
    #print(f'               {len(prot_pep_index)}     -     {len(prot_CDS_index)}')

    RGB_tup = generate_color_gradient(color_lst=color_gradient, reverse_gradient=False)  # 'gray',

    # ----------------------- #
    # DRAW THE COLOR GRADIENT #
    # ----------------------- #
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.yaxis.set_visible(False)
    # x_tick_locations=range(0,max_intensity-min_intensity,10)
    # x_tick_labels=range(min_intensity,max_intensity,1000)
    # plt.xticks(x_tick_locations, x_tick_labels)
    for x, colorRGB in enumerate(RGB_tup):
        ax.axvline(x, color=colorRGB, linewidth=4)
    plt.show()
    # ----------------------- #

    # PRINT THE RGB TUPLES FOR THE COLOR GRADIENT
    print('RGB TUPLES\n----------')
    rgbind = 0
    for rgb in RGB_tup:
        print(rgbind, '  -  ', round(rgb[0], 3), '-', round(rgb[1], 3), '-', round(rgb[2], 3))
        rgbind += 1

    # +++++++++++++++++++++++++++++++++++++++++++ #
    #      Convert intensities in RGB codes       #
    prot_vec, RGB_vector = exprlev_resc_RGB(RGB_tup)
    # +++++++++++++++++++++++++++++++++++++++++++ #

    # print(f'{len(prot_vec)} - {len(prot_CDS_index)}')
    # print(prot_expressions_RGB)
    # print(prot_CDS_index)

    # Update the dictionary of prot_CDS_index with the RGB intensities
    print('\nUPDATING PROTEIN INDEX WITH RGB INTENSITIES\n-------------------------------------------')
    ind = 0
    for ind, prot in enumerate(prot_vec):
        RGB_tup = RGB_vector[ind]
        RGB_code = str(RGB_tup[0]) + ',' + str(RGB_tup[1]) + ',' + str(RGB_tup[2])
        prot_PSMint_index[prot].append(RGB_code)  # prot_expressions_RGB[ind]
        print(prot, '-', prot_PSMint_index[prot][-1])  # prot_expressions_RGB[ind]
        # ind+=1
    return prot_PSMint_index
