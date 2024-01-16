from Proteogenome3 import pg_data_processing as pg_dp
from Proteogenome3 import pg_utils
from Proteogenome3 import pg_input

import pandas as pd
import numpy as np
import re

def initialise_indexes(prot_annots_path, input_table, annot_format='gff_compliant'):
    """
    Version: 1.1
    Name History: initialise_indexes

    This function create the index for (the '--->' means 'to')
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

    prot_pep_index = protein_peptide_index(input_table)
    print('prot_pep_index - DONE')
    pep_prot_index = peptide_protein_index(input_table)
    print('pep_prot_index - DONE')
    unique_protein_list = prot_pep_index.keys()

    if annot_format == 'not_gtf_compliant':
        genome_annot_lst = pg_input.file_to_lst(prot_annots_path)                             # Upload the genome annotation file into a list
        CDS_matrix = pg_dp.CDS_annot_matrix(genome_annot_lst)                                 # Matrix with CDS annotations info only
        prot_CDS_index = gen_protein_CDS_index(CDS_matrix, annot_format='not_gtf_compliant')  # Dictionary of proteins and their group of CDS
        return CDS_matrix, prot_CDS_index, prot_pep_index, pep_prot_index
    else:
        ensembl_prot_CDS_index = {}
        CDS_matrix = []
        prot_CDS_index = gen_protein_CDS_index(prot_annots_path, annot_format=annot_format)
        for ensembl_acc in unique_protein_list:
            prot_block = prot_CDS_index[prot_CDS_index['attribute'].str.contains(ensembl_acc)] \
                                                                                 [['seqname', 'start', 'end', 'strand']]
            if prot_block.empty == False: # Consider only the the protein with CDS. UniProt accession not converted will be skipped
                ensembl_prot_CDS_index[ensembl_acc] = prot_block.values.tolist()
        print(ensembl_prot_CDS_index)
        print(f'\nCDS_matrix\n----------\n{CDS_matrix}')
        # a=input('Look at the prot_CDS ----- STOP')
        return CDS_matrix, ensembl_prot_CDS_index, prot_pep_index, pep_prot_index
                # try:
    #     prot_pep_index = protein_peptide_index(input_table)
    #     pep_prot_index = peptide_protein_index(input_table)
    # except:
    #     print('Proteomics data not found')



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


def gen_protein_CDS_index(annotations, annot_format=''):
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

    # **************** Parsing annotations in GFF format
    if annot_format == 'not_gtf_compliant':  # specific patterns for the
        # gene_pat = re.compile(r'.*?\sgene\s\"(.*?)\";')          # annotation in GTF format
        gene_pat = re.compile(r'.*?\sprotein_id\s\"(.*?)\";')       # annotation in GTF format
        # gene_pat = re.compile(r'.*?;gene=(.*?);')        # annotation in GFF3 format
        # gene_pat = re.compile(r'.*?;protein_id=(.*?);')  # annotation in GFF3 format
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
            print(row)
            print(protein_ID)
            print('-----------------------')

            CDS_feat = [coord_1, coord_2, strand]   ---- O
            # CDS_feat = [int(coord_1), int(coord_2), strand]

            if (protein_ID in prot_CDS_index):
                # Append the features of the current CDS on the list of CDS that belongs to the current protein.
                prot_CDS_index[protein_ID].append(CDS_feat)
            else:
                prot_CDS_index[protein_ID] = [CDS_feat]

    # **************** The GTF format is converted straight to the CDS index
    elif annot_format == 'gtf_compliant':
        col_names = ["seqname","source","feature","start","end","score","strand","frame","attribute"]
        annot_to_df = pd.read_csv(annotations, sep="\t", names=col_names)
        prot_CDS_index = annot_to_df[annot_to_df['feature'] == "CDS"]                # Extract CDS

    print(f'\n***********\nprot_CDS_index\n***********\n{prot_CDS_index}')

    return prot_CDS_index


def protein_PSM_int_index(prot_pep_index, peptides_input_table, color_gradient,
                          quantitation_method = 'log2Fold-1', coloring_method = 'volcano',
                          ):
#=['black', 'blue', 'cyan', 'green', 'greenyellow', 'yellow', 'orange', 'red']

    """
    Version: 2.0

    Name History: protein_PSM_int_index

    This function creates the protein index for the PSM and intensities.
    Moreover, it generates the RGB code for each protein intensity.
    The RGB codes will be used for the creation of the protein map.

    :param      prot_pep_index    :
    :param      prot_CDS_index    :
    :return     prot_PSMint_index :
    """
        #  ------- MAIN --------
        #  protein_PSM_int_index

    treshold_l2f = 0
    intensities = {}

    max_intensity = np.min(peptides_input_table[:,5])
    min_intensity = np.min(peptides_input_table[:,5])

    prot_PSMint_index = {}
    if ('log2Fold' in quantitation_method) or ('pag' in quantitation_method):
        if ('log2Fold' in quantitation_method):
            treshold_l2f = float(quantitation_method.split('-')[1]) # Only for log2fold extract the threshold value
        sub_peptides_input_table = peptides_input_table[:, np.r_[0, 5]]
        # sub_peptides_input_table = sub_peptides_input_table[abs(sub_peptides_input_table[:,1]) >= threshold_l2f]
        prot_tup_set = [tuple(row) for row in sub_peptides_input_table]
        prot_tup_set = list(set(prot_tup_set))
        prot_tup_set = np.array(prot_tup_set, dtype = 'object')

        prot_abun_float =  list(map(float, prot_tup_set[:, 1]))
        prot_abundance = np.vstack([prot_tup_set[:, 0], prot_abun_float]).transpose()
        for protid_abundance in prot_abundance:
            prot_PSMint_index[protid_abundance[0]] = [1, protid_abundance[1]]
        prot_PSMint_index = pg_utils.dict_numlst_2_dict_strlst(prot_PSMint_index)

    elif quantitation_method == 'pas': # Protein Abundance = Peptides Abundance Sum
        for protein, pep_array in prot_pep_index.items():
            PSM_sum = 0
            inten_sum = 0
            for pep_row in pep_array:
                PSM_sum += int(pep_row[3])
                inten_sum += float(pep_row[4])
            prot_PSMint_index[protein] = [PSM_sum, inten_sum]

            if max_intensity < inten_sum: max_intensity = inten_sum
            if min_intensity > inten_sum: min_intensity = inten_sum

    if coloring_method == 'continuous_color_fade':
        # print('Protein-Peptide Index - Protein-CDS Index')
        # print(f'               {len(prot_pep_index)}     -     {len(prot_CDS_index)}')

        RGB_tup = pg_utils.generate_color_gradient(color_lst=color_gradient, reverse_gradient=False)  # 'gray',

        print('Color gradient done\n+++++++++++++++++++\n ')

        # PRINT THE RGB TUPLES FOR THE COLOR GRADIENT
        print('RGB TUPLES\n----------')
        rgbind = 0
        for rgb in RGB_tup:
            print(rgbind, '  -  ', round(rgb[0], 3), '-', round(rgb[1], 3), '-', round(rgb[2], 3))
            rgbind += 1

        # +++++++++++++++++++++++++++++++++++++++++++ #
        #      Convert intensities in RGB codes       #
        prot_vec, RGB_vector = pg_utils.exprlev_resc_RGB(RGB_tup, prot_PSMint_index)
        # +++++++++++++++++++++++++++++++++++++++++++ #

        # Update the dictionary of prot_CDS_index with the RGB intensities
        print('\nUPDATING PROTEIN INDEX WITH RGB INTENSITIES\n-------------------------------------------')
        ind = 0
        for ind, prot in enumerate(prot_vec):
            RGB_tup = list(
                map(pg_utils.roundfloat, RGB_vector[ind]))  # round the RGB values (3 of them in the RGB_vector)
            RGB_code = str(RGB_tup[0]) + ',' + str(RGB_tup[1]) + ',' + str(RGB_tup[2])
            prot_PSMint_index[prot].append(RGB_code)
            print(prot, '-', prot_PSMint_index[prot][-1])  # prot_expressions_RGB[ind]
            # ind+=1

    elif coloring_method == 'volcano':
        for prot, psm_int in prot_PSMint_index.items():
            prot_PSMint_index[prot].append(pg_utils.volcano_coloring(psm_int[1], treshold=treshold_l2f))

    return prot_PSMint_index
