def CDS_annot_matrix(annotations):
    """
    Version: 1.0

    Name History: annot_matrix, CDS_annot_matrix

    This function receives a list of annotations rows and extracts only the block of
    rows regarding the protein CDS
    INPUT :
    OUTPUT:
            self.annotation_matrix
            self.CDS_matrix
    """

    first = True
    for row in annotations:
        if row[0] != '#':

            if first:
                row_array = np.array(row.split('\t'), dtype=object)  # Split the current annotation rows in columns
                self.annotation_matrix = np.array([row_array])  # increase the annotation matrix dimension
                first = False
            else:
                row_array = np.array(row.split('\t'), dtype=object)

                self.annotation_matrix = np.concatenate((self.annotation_matrix, [row_array]))
                # print('*'*50)
                # print(self.annotation_matrix)
                # a=input()

    self.CDS_matrix = self.annotation_matrix[self.annotation_matrix[:, 2] == 'CDS']  # Filter only for coding regions



def annot_to_df(annotations, annot_format ='gff3'):
    """
    Version: 1.0

    Name History: annot_to_df

    This function extracts relevant data from the annotation file and put them in a DataFrame

    INPUT : annotations     List    List list filled with the annotation rows
    OUTPUT:
    """
    #annotations_fh = open(GFF3,'r')

    protein_index = {}
    self.annotations_df = pd.DataFrame()

    new_row = {}      # In order to append a new row to the dataframe I need a dictionary

    # **************** Parsing file for GFF3 format
    if annot_format =='gff3':
        unique_pat = re.compile(r'(.*?)\t')
        ID_pat = re.compile(r'.*?\tID=(.*?);')            # specific patterns for the
        gene_pat = re.compile(r'.*?;gene=(.*?);')         # annotation in GFF3 format
        product_pat = re.compile(r'.*?;product=(.*?)$')
        col_index=[2, 3, 4, 6]

    elif annot_format == 'gtf':
        unique_pat = re.compile(r'(.*?)\t')
        ID_pat = re.compile(r'.*?\sgene_id\s\"(.*?)\";')   # specific patterns for the
        gene_pat = re.compile(r'.*?;\sgene\s\"(.*?)\";')   # annotation in GTF format
        product_pat = re.compile(r'.*?;\sproduct\s\"(.*?)\"$')
        col_index=[2, 3, 4, 6]

    for row in annotations:

        if (row[0]!='#') and ('country=United Kingdom: Cardiff;culture-collection=ATCC' not in row):
            match_lst=unique_pat.findall(row)

            # print(row)
            # a=input()
            try:
                new_row['row_type'] = match_lst[col_index[0]]      # 2 - GFF3 column 3
                new_row['coordinate1'] = match_lst[col_index[1]]   # 3 - GFF3 column 4
                new_row['coordinate2'] = match_lst[col_index[2]]   # 4 - GFF3 column 5
                new_row['strand'] = match_lst[col_index[3]]        # 6 - GFF3 column 7
            except:
                print('This line returns an error on parsing procedure!!!')
                print(f'-------------------------------------------------\n {row}')
                print('Please press a button to continue')
                #a=input()

            try:
                new_row['ID'] = ID_pat.match(row).group(1)
            except:
                new_row['ID'] = 'ID nf'
            try:
                new_row['gene'] = gene_pat.match(row).group(1)
            except:
                new_row['gene'] = 'gene nf'
            try:
                new_row['product'] = product_pat.match(row).group(1)
                new_row['product'] = new_row['product'].split(';')[0]
            except:
                new_row['product'] = 'product nf'

            #print(new_row)

            self.annotations_df = self.annotations_df.append(new_row, ignore_index=True)
    self.annotations_df = self.annotations_df[['row_type', 'coordinate1', 'coordinate2',
                                               'strand', 'ID', 'gene', 'product']]
    #annotations_fh.close()

    #return self.GFF3_to_df

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

    PTMs_code = {'phosphorylation': '(phospho)', 'acetylation': '(acetyl)', 'acetylation': '(acetyl)',
                 'amidation': '(amidated)',
                 'oxidation': '(oxidation)', 'methylation': '(methyl)', 'ubiquitinylation': '(glygly; gg)',
                 'sulfation': '(sulfo)',
                 'palmitoylation': '(palmitoyl)', 'formylation': '(formyl)', 'deamidation': '(deamidated)'
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



def filter_peptides(PoGo_peptides, out_file_name=''):
    """
    Version: 1.2

    Name History: filter_peptides

    This function  filter the peptides mapped by PoGo.
    The filter criteria is based on the proteins identified in the proteomics data
    provided to the software.
    The peptides coordinates are compared with the CDS coordinates.
    If the peptide coordinates between the CDS coordinates, then the peptides will
    be included into the peptide map.
    If the PTMs have been applied to the peptide sequences the function detects the
    modifications, extracts the original peptides sequences and performs the peptide
    filtration. However, the peptide sequences already updated with the PTM encoding
    will not be modified in all instance's attributes.

    INPUT : out_file_name   [Str]   File name for the output file
    OUTPUT:
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

        protein_block = self.pep_prot_index[
            pep_seq]  # Fetch the protein codes where the current peptide has been found
        protein_set = np.concatenate((protein_set, protein_block))
    protein_set = np.unique(protein_set)  # Shrink the protein code collection to the unique codes

    ################### GENERATE THE ALLOWED GENOMIC SPACE ###################
    allowed_genomic_space = {}
    for unique_prot_code in protein_set:
        int_CDS_coordinates = []
        for str_CDS_coord in self.prot_CDS_index[
            unique_prot_code]:  # Considering the current protein code coordinates.
            from_str_to_int = list(map(int, str_CDS_coord[0:2]))  # Convert from Str to Int all the CDS coordinates.
            from_str_to_int.append(str_CDS_coord[-1])  # Append the strand tag in a Str format
            int_CDS_coordinates.append(from_str_to_int)  # Increase the coordinate bolck
        allowed_genomic_space[
            unique_prot_code] = int_CDS_coordinates  # Append the new piece of allowed genomic coordinates.

    print(allowed_genomic_space)
    ##########################################################################

    for pep_row_index, peptide_row in enumerate(coord_pep_strand):  # Iterate over the PoGo peptide map
        peptide_coord_1 = int(peptide_row[0])  # Fetch peptide genomic coordinates
        peptide_coord_2 = int(peptide_row[1])
        peptide_sequence = peptide_row[2]
        if '(' in peptide_sequence:
            peptide_sequence = re.sub('\(.*\)', '', peptide_sequence)
        peptide_strand = peptide_row[3]
        peptide_to_protein = self.pep_prot_index[
            peptide_sequence]  # Fetch the set of proteins where the peptide has been found.

        # ---------- COORDINATES COMPARISON ---------- #

        for protein in peptide_to_protein:  # Iterate the set of protein

            CDS_block = allowed_genomic_space[
                protein]  # Fetch the genomic coordinates of the CDS of the protein where the peptide has been found.
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
    # print('PoGo_peptides')
    # print(PoGo_peptides)

    if out_file_name:
        self.make_sep_file(out_file_name, PoGo_peptides, sep='\t')
        # -------------------------------------------- #


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
