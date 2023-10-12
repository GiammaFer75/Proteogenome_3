from Proteogenome3 import pg_output
from Proteogenome3 import pg_utils

import pandas as pd
import numpy as np
import re

def CDS_annot_matrix(annotations_in_lst):
    """
    Version: 1.1
    Name History: annot_matrix, CDS_annot_matrix

    This function receives a list of annotations rows and extracts only the block of
    rows regarding the protein CDS
    INPUT : annotations_in_lst  List        Contains the rows of the genome annotation file
    OUTPUT: annotation_matrix   np.array    Matrix with the annotation splitted and fitered to return only the CDS
                                            information.
    """

    first = True
    for row in annotations_in_lst:
        if row[0] != '#':

            if first:
                row_array = np.array(row.split('\t'), dtype=object)  # Split the current annotation rows in columns
                annotation_matrix = np.array([row_array])  # increase the annotation matrix dimension
                first = False
            else:
                row_array = np.array(row.split('\t'), dtype=object)

                annotation_matrix = np.concatenate((annotation_matrix, [row_array]))
                # print('*'*50)
                # print(self.annotation_matrix)
                # a=input()

    return annotation_matrix[annotation_matrix[:, 2] == 'CDS']  # Filter only for coding regions


def annot_to_df(annotations_path, annot_format ='gff3'):
    """
    Version: 2.0
    Name History: annot_to_df

    This function upload the genome annotation in a data frame

    :param      annotations_path    Str         The file name of the GFF/GTF genome annotation
    :return     annotation_df       DataFrame   The data frame with the annotation file content
    """
    file_type = annotations_path[-3:]
    gtf_gff_column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    bed_column_names = ["chrom", "start", "end", "name", "score", "strand", "thickStart",    "thickEnd", "itemRrgb",
                        "blockCount", "blockSizes", "blockStarts"]
    if (file_type == "gtf") or (file_type == "GTF") or (file_type == "gff") or (file_type == "GFF"):
        column_names = gtf_gff_column_names
    elif (file_type == "bed") or (file_type == "BED"):
        column_names = bed_column_names
    annotation_df = pd.read_csv(annotations_path, sep = "\t", names = column_names)
    return annotation_df

def apply_PTMs_to_pep(peptide_tab, PTMs_to_remove=[]):
    """
    Version: 2.0
    Name History: apply_PTMs_to_pep

    This function updates the peptide sequences with the PTMs encoding (see PoGo software instructions).

    INPUT  :  peptide_tab              np.array   The table with proteomics data.

              PTMs_to_remove           List       Strings with the type of PTMs that must not be applied to the peptides.
    OUTPUT :  peptides_PTMs_updated    List       he peptides sequence updated with the PTMs encoding (format in PoGo instruction)

    References: PoGo - https://github.com/cschlaffner/PoGo
    """

    def PTMs_type_set(modifications, PTMs_remove):
        """
        Version: 2.0
        Name History: PTMs_type_set
        """

        # Apply the unique function but the result might be something like that [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...] because the PTM |amino acid(position)| part might affect the uniqueness of the PTM type (Oxidation in this example)
        PTM_types = []

        if 'None' in modifications:  # If 'None' value is included into the PTM values.
            modifications = list(dict.fromkeys(modifications))  # Clean from 'None' duplications.
            None_position = modifications.index('None')  # Find the position of the last 'None' value.
            del modifications[None_position]  # Remove the last 'None' value.

        for PTM in modifications:  # Iterate over the PTMs list. Because the PTMs could appear like that Example: [Oxidation M(6), Oxidation M(10), Oxidation M(1),  ...]
            # we are interested only in the PTM type (in the example: 'Oxidation').
            PTM_types.append(PTM.split()[
                                 0])  # Splitting the original PTM and taking [0] we extract always the type (  Example of split result ['Oxidation', 'M(10)']  )
        PTM_types = list(
            set(PTM_types))  # Apply the unique (set) function for the list in order to have only one element for each PTM type.

        for mod_type in PTMs_remove:  # Loop over the list of PTMs to exclude
            if (mod_type in PTM_types): PTM_types.remove(mod_type)  # Remove the PTM not desired

        return PTM_types

    # print(' --------------- PTMs_to_remove ---------------\n',PTMs_to_remove)
    # print(' ***************** peptide_tab *****************\n',peptide_tab)

    PTMs_code = {'phosphorylation': '(phospho)', 'phospho': '(phospho)','acetylation': '(acetyl)', 'acetyl': '(acetyl)',
                 'amidation': '(amidated)', 'oxidation': '(oxidation)', 'methylation': '(methyl)',
                 'methyl': '(methyl)', 'ubiquitinylation': '(glygly; gg)', 'sulfation': '(sulfo)',
                 'palmitoylation': '(palmitoyl)', 'formylation': '(formyl)', 'deamidation': '(deamidated)',
                 'nitrosyl': '(oxidation)'
                 }  # Any other post-translational modification

    peptides_PTMs_updated = np.array([[]])
    PTM_position_pat = re.compile(r'.*?\([A-Z](.*?)\)')

    # Extract the column with the PTM encoding.
    PTM_column = peptide_tab[:, 1]

    # Create a list of valid PTMs to insert in the peptide sequences.
    PTMs_types = PTMs_type_set(PTM_column, PTMs_to_remove)

    for peptide_ind, peptide_row in enumerate(peptide_tab):  # Loop over the peptide table.
        peptide_PTM = peptide_row[0]  # From the row take only the peptide SEQUENCE.

        for apply_PTM in PTMs_types:  # Loop over the PTMs list.

            peptide_modification = peptide_row[1]  # From the row take only the peptide MODIFICATION
            if peptide_modification != 'None':
                peptide_modification = peptide_modification.split()  # Split the peptide modification string in the two main components: Modification Type - Modification Position
                modificatio_type = peptide_modification[0]  # Modification Type

                # try:
                modification_position = peptide_modification[
                    1]  # Modification Position. Only for those modification other than 'None'
                if (modificatio_type.lower() == apply_PTM.lower()):  # Find the specific PTM allowed to be applyed
                    modification_position = int(PTM_position_pat.match(modification_position).group(1))

                    # modification_position = int(modification_position.split('(')[1].replace(')',''))  # Clean the Modification position. Example: in general modification position might appear like that M(10) and we want only '10'
                    PTM_encoded = PTMs_code[modificatio_type.lower()]  # Translate the PTM name in the PoGo PTM name
                    peptide_PTM = peptide_PTM[:modification_position] + PTM_encoded + peptide_PTM[
                                                                                      modification_position:]  # Insert the PTM encoded in the peptide sequence

                # except:
                #    pass  # If you are on a peptide that reports 'None' in the modification field you have to skip the entire updating part

        peptide_tab[peptide_ind, 0] = peptide_PTM

    return peptide_tab



def filter_peptides(PoGo_peptides, pep_prot_index, prot_CDS_index, out_file_name=''):
    """
    Version: 2.0
    Name History: filter_peptides

    This function  filter the peptides mapped by PoGo.
    The filter criteria is based on the proteins identified in the proteomics data
    provided to Proteogenome.
    The peptides coordinates are compared with the CDS coordinates.
    If the peptide coordinates between the CDS coordinates, then the peptides will
    be included into the peptide map.
    If the PTMs have been applied to the peptide sequences the function detects the
    modifications, extracts the original peptides sequences and performs the peptide
    filtration. However, the peptide sequences already updated with the PTM encoding
    will not be modified.

    :param      out_file_name      String     File name for the output file :
    :param      PoGo_peptides                                               :
    :param      pep_prot_index                             :
    :param      prot_CDS_index                  :
    :return:
    """

    row, col = PoGo_peptides.shape

    # Initialise a smaller peptide table #
    coord_pep_strand = PoGo_peptides[:, 1:6]  # Slice the PoGo mapped peptide table considering only:
    # | start | end | peptide sequence | score | strand |
    coord_pep_strand = np.delete(coord_pep_strand, 3, 1)  # Delete the 'score' column

    print('Subsetting of PoGo_peptides\n', coord_pep_strand)

    protein_set = np.array([])  # Set of unique protein codes

    peptide_seq_array = coord_pep_strand[:, 2]  # Consider all the peptide sequence
    print('peptide_seq_array ----\n', peptide_seq_array)
    for pep_seq in peptide_seq_array:

        # If a '(' is in the peptide sequence, then means that the PTMs have been already applied to the peptide sequences.
        if '(' in pep_seq:
            pep_seq = re.sub('\(.*\)', '',
                             pep_seq)  # Remove the PTM encoding restoring the original peptide sequence

        protein_block = pep_prot_index[pep_seq]  # Fetch the protein codes where the current peptide has been found
        protein_set = np.concatenate((protein_set, protein_block))
    protein_set = np.unique(protein_set)  # Shrink the protein code collection to the unique codes
    print(f'protein_set - {protein_set}')
    ################### GENERATE THE ALLOWED GENOMIC SPACE ###################
    allowed_genomic_space = {}

    # Check the type of data in the first position of the CDS array
    p_code = list(prot_CDS_index.items())[0][1] # Extract the first CDS record
    print(f'p_code - {p_code}')
    first_cds = p_code[0]
    print(f'first_cds - {first_cds}')
    cpo, cpt = None, None
    if first_cds[0].replace('.','').isdigit():
        cpo, cpt = 0, 2   # For virus there is no chr number
    else:
        cpo, cpt = 1, 3   # For human the chr number is in position 0. Then take position 1 and 2 in the array


    for unique_prot_code in protein_set:
        int_CDS_coordinates = []
        try:
            for str_CDS_coord in prot_CDS_index[unique_prot_code]:  # Considering the current protein code coordinates.
                from_str_to_int = list(map(int, str_CDS_coord[cpo:cpt]))  # Convert from Str to Int all the CDS coordinates.
                from_str_to_int.append(str_CDS_coord[-1])  # Append the strand tag in a Str format
                int_CDS_coordinates.append(from_str_to_int)  # Increase the coordinate bolck
            allowed_genomic_space[unique_prot_code] = int_CDS_coordinates  # Append the new piece of allowed genomic coordinates.
        except Exception as excpt:
            print(f'ERROR - {unique_prot_code}')
            print(excpt)

    print(f'allowed_genomic_space\n---------------------\n{allowed_genomic_space}')
    ##########################################################################

    for pep_row_index, peptide_row in enumerate(coord_pep_strand):  # Iterate over the PoGo peptide map
        peptide_coord_1 = int(peptide_row[0])  # Fetch peptide genomic coordinates
        peptide_coord_2 = int(peptide_row[1])
        peptide_sequence = peptide_row[2]
        if '(' in peptide_sequence:
            peptide_sequence = re.sub('\(.*\)', '', peptide_sequence)
        peptide_strand = peptide_row[3]
        peptide_to_protein = pep_prot_index[peptide_sequence]   # Fetch the set of proteins where
                                                                # the peptide has been found.

        # ---------- COORDINATES COMPARISON ---------- #

        for protein in peptide_to_protein:  # Iterate the set of protein
            try:
                CDS_block = allowed_genomic_space[protein]  # Fetch the genomic coordinates of the CDS of the protein where the peptide has been found.
                for CDS in CDS_block:
                    CDS_coord_1 = CDS[0]
                    CDS_coord_2 = CDS[1]
                    if peptide_strand == '+':
                        # if (peptide_coord_1>=CDS_coord_1) and (peptide_coord_2<=CDS_coord_2): # VALID genomic coordinates
                        #     pass
                        # else:                                                                 # INVALID genomic coordinates
                        #     PoGo_peptides=np.delete(PoGo_peptides,pep_row_index,0)     # REMOVE THE PEPTIDE FROM THE PoGo PEPTIDE TABLE

                        if (peptide_coord_1 < CDS_coord_1) and (
                                peptide_coord_2 > CDS_coord_2):  # INVALID genomic coordinates
                            try:
                                # PoGo_peptides=np.delete(PoGo_peptides,pep_row_index,0)     # REMOVE THE PEPTIDE FROM THE PoGo PEPTIDE TABLE
                                PoGo_peptides[pep_row_index, 0] = ''
                            except:
                                print('Invalid row index -----> ', PoGo_peptides.shape, '-', pep_row_index)
            except:
                print(f'Protein {protein} NOT FOUND. Genomic coordinates for this position of the peptide in the genome skipped.')
    # print('PoGo_peptides')
    # print(PoGo_peptides)

    if out_file_name:
        # print(out_file_name)
        # filtered_PoGo_peptides_file_path = pg_utils.create_path(out_file_name, add_prefix='proteogenome3_')
        # print(filtered_PoGo_peptides_file_path)
        pg_output.make_sep_file(out_file_name, PoGo_peptides, sep='\t')




def groupby_peptide(pep_input_table, peptide_sequence):
    """
    Version: 1.0

    Name History: groupby_peptide

    INPUT : pep_input_table
            peptide_sequence
    OUTPUT: protein_peptide_block
    """
    peptide_protein_block=pep_input_table[pep_input_table[:,1]==peptide_sequence]
    peptide_protein_block=peptide_protein_block[:,0]
    return peptide_protein_block


def groupby_protein(pep_input_table, protein_ID):
    """
    Version: 1.0

    Name History: groupby_protein

    INPUT :
    OUTPUT:
    """
    protein_peptide_block=pep_input_table[pep_input_table[:,0]==protein_ID]
    return protein_peptide_block
