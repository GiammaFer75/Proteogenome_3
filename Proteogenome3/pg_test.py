def dummy_peptides(prot_sequence, pep_min_length, pep_max_length):
    """
    Version: 1.0

    Name History: dummy_peptides

    INPUT :

    OUTPUT:
    """

    prot_sequence = prot_sequence.replace('\n', '')
    prot_sequence = prot_sequence.replace('K', '|')  # LYSINE
    prot_sequence = prot_sequence.replace('R', '|')  # ARGININE

    new_peptides = prot_sequence.split('|')  # The protein sequence is now a list of peptides

    peptides = []
    for pep in new_peptides:  # Filter the peptides with the lenght raange
        if (len(pep) >= pep_min_length) & (len(pep) <= pep_max_length):
            peptides.append(pep)

    return peptides

# def dummy_input(peptides='', out_file_name='', dummy_exp_name = 'exp1',
#                 prot_sequences_FASTA_fn='', pep_min_length=4, pep_max_length=30,
#                 PSMs_range=1, pep_int_range=[100,10000]):
#     """
#     Version: 1.0

#     Name History: dummy_PoGo_peptides - dummy_input

#     This function generates a file .txt that contains dummy peptides in a format suitable for PoGo software.

#     INPUT : peptides            List[Str]       List of peptide sequence.
#             out_file_name       Str             Name of the file that will contains the dummy peptides.
#             dummy_exp_name      Str             Dummy experiment name.
#     OUTPUT: dummy_input_matrix  np.array[Str]   The matrix with the dummy peptides.
#     """

#     import copy

#     prot_ID_pat = re.compile(r'.*?gene=.*?(.*?)\slocus_tag=|.*?gene:(.*?)\s') # Pattern for protein ID

#     #FASTA_in_lst = self.file_to_lst(prot_sequences_FASTA_fn) # Upload the FASTA file with protein sequences in a list
#     #FASTA_in_lst = self.FASTA_cpt_seq(self.FASTA_lst)         # Compact possible multilines sequences

#     FASTA_in_lst= copy.deepcopy(self.FASTA_lst)

#     dummy_input_matrix = np.array([['','','','','']], dtype=object)
#     table_row=[]
#     current_protein=True # This boolean signals if the data belong to the current protein.
#                          # This is used for combine the protein ID from the FASTA header with his sequence.

#     for ind,prot_sequence in enumerate(FASTA_in_lst): # Iterate over the list that contains the FASTA file rows.
#         #if ind>4: break # ################REMOVE THIS ROW IN THE FINAL VERSION################

#         if prot_sequence[0] == '>':# If this row is a FASTA headers.
#             # print(prot_sequence)
#             # a=input()
#             protein_ID=prot_ID_pat.match(prot_sequence).group(1) # Extract the PROTEIN ID
#             current_protein=True       # Signals that the next rows belongs to this protein.

#         else:                   # Generate the group of PEPTIDES that belong to the current protein.
#             current_prot_pep=self.dummy_peptides(prot_sequence,pep_min_length, pep_max_length)
#             current_protein=False      # Signals that the next row will belong to a new protein.

#         if not current_protein:
#             for new_peptide in current_prot_pep: # Iterate over the peptides generated from the current protein sequence.

#                 current_prot_pep=np.array([protein_ID,new_peptide],dtype=object) # Convert the list in a np.array.


#             # ***** RANDOM PSMs ***** #
#                 if type(PSMs_range)==list:
#                                       # lower | upper
#                                       # bound | bound
#                     rand_PSM=randrange(PSMs[0],PSMs[1]) # Generate random PSMs
#                 else:
#                     rand_PSM=1
#                 rand_PSM=np.array(['',rand_PSM])
#                 #current_prot_pep=np.append(current_prot_pep,[rand_PSM],axis=0)
#                 current_prot_pep=np.append(current_prot_pep,rand_PSM,axis=0)

#             # ***** RANDOM INTENSITY ***** #
#                 rand_intensity=np.array(randrange(pep_int_range[0],pep_int_range[1]))
#                 current_prot_pep=np.append(current_prot_pep,[rand_intensity],axis=0)
#                 dummy_input_matrix=np.concatenate((dummy_input_matrix, [current_prot_pep]))

#     dummy_input_matrix=np.delete(dummy_input_matrix,0,0) # Remove the first initialisation row


#     #************* Output File creation
#     if out_file_name !='':

#         fh = open(out_file_name, 'w')
#         for ind, row in enumerate(dummy_input_matrix):
#             write_row = ''
#     #         print(pep)
#             write_row += dummy_exp_name + '\t' + row[1] + '\t' + str(row[2]) + '\t' + str(row[3]) + '\n'
#             # write_row += dummy_exp_name + '\t' + pep + '\t' + str(PSM[ind]) + '\t' + str(inten[ind]) + '\n'
#     #         print(row)
#             fh.write(write_row)
#     #         a=input()
#         fh.close()

#     return dummy_input_matrix

def generate_dummy_peptides(pep_table_fn='dummy_peptide_table.txt',
                            ID_search_pattern='.*?gene=(.*?)\s',
                            dummy_exp_name='exp1', pep_min_length=4, pep_max_length=30,
                            PSMs_range=1, pep_int_range=[100, 10000],
                            pept_percentage=0, PTM_percentage=0):
    """
    Version: 2.0

    Name History: dummy_PoGo_peptides - dummy_input - dummy_input_nparrays - generate_dummy_peptides

    This function generates a file .txt that contains dummy peptides in a format suitable for PoGo software.
    In order to create this peptides the function needs the protein sequences and annotations from a reference
    organism in FASTA and GFF3 format respectively


    INPUT : pep_table_fn        [Str]   The file name where will be saved the dummy proteomics data.
            ID_search_pattern   [Str]   The pattern in regex format that match the protein IDs in the FASTA file.
            dummy_exp_name      [Str]   Dummy experiment name related to the simulated proteomics data.
            pep_min_length      [Int]   The minimum length allowed in the random generation of the peptides.
            pep_max_length      [Int]   The maximum length allowed in the random generation of the peptides.
            PSMs_range=1  by default
            pep_int_range       [List]  Upper and lower bound respectively of the random range used to generate peptide intensities.
            pept_percentage     [Int]   How many peptides must be considerd over the total number of dummy peptides generated.
            PTM_percentage      [Int]   From the peptides in the pept_percentage, how many of them will report a random PTM.
    OUTPUT:
    """
    PTMs_types = ['phosphorylation', 'acetylation', 'acetylation', 'amidation', 'oxidation', 'methylation',
                  'ubiquitinylation', 'sulfation', 'palmitoylation', 'formylation', 'deamidation']

    prot_ID_pat = re.compile(r'' + ID_search_pattern + '')  # Pattern for protein ID

    # FASTA_in_lst = self.file_to_lst(prot_sequences_FASTA_fn) # Upload the FASTA file with protein sequences in a list
    # FASTA_in_lst = self.FASTA_cpt_seq(self.FASTA_lst)         # Compact possible multilines sequences

    # self.dummy_input_matrix = np.array([['','','','']])
    self.dummy_input_matrix = np.empty((0, 2))
    table_row = []
    current_protein = True  # This boolean signals if the data belong to the current protein.
    # This is used for combine the protein ID from the FASTA header with his sequence.

    for ind, prot_sequence in enumerate(self.FASTA_lst):  # Iterate over the list that contains the FASTA file rows.
        # if ind>4: break # ################REMOVE THIS ROW IN THE FINAL VERSION################
        if prot_sequence[0] == '>':  # This condition is for the FASTA headers.
            protein_ID = prot_ID_pat.match(prot_sequence).group(1)  # Extract the PROTEIN ID

            current_protein = True  # Signals that the next rows belongs to this protein.
        else:  # Generate the group of PEPTIDES that belong to the current protein.
            current_prot_pep = self.dummy_peptides(prot_sequence, pep_min_length, pep_max_length)
            current_protein = False  # Signals that the next row will belong to a new protein.

        if not current_protein:  # If the boolean for current protein is false means that is time to generate table rows.
            for new_peptide in current_prot_pep:  # Iterate over the peptides generated from the current protein sequence.

                current_prot_pep = np.array([protein_ID, new_peptide],
                                            dtype=object)  # Convert the list in a np.array.
                # print(current_prot_pep)
                # else:
                #     PSMs_col=np.ones(np.shape(self.dummy_input_matrix)[0]) # Otherwise the PSM value is 1 for all the peptides

                # print(current_prot_pep)
                # print(self.dummy_input_matrix)
                # print(self.dummy_input_matrix.shape,current_prot_pep.shape)
                self.dummy_input_matrix = np.append(self.dummy_input_matrix, [current_prot_pep], axis=0)

    total_peptides = np.shape(self.dummy_input_matrix)[0]

    print('INPUT MATRIX')
    print(f'Type input matrix - {type(self.dummy_input_matrix)}')
    print(f'Shape input matrix - {np.shape(self.dummy_input_matrix)}')
    print(f'input matrix -\n{self.dummy_input_matrix}')

    # ******* GENERATE PEPTIDE PSMs ******* #

    if type(PSMs_range) == list:  # lower | upper  | Generates an array of random values with the
        # bound | bound  | same number of rows present in the input table.
        PSMs_col = np.random.randint(PSMs[0], PSMs[1], size=np.shape(total_peptides))
    else:
        PSMs_col = np.ones(total_peptides, dtype=np.int8)  # Otherwise the PSM value is 1 for all the peptides
    PSMs_col = np.reshape(PSMs_col, (-1, 1))  # The unspecified value is inferred to be the number of rows
    print('\n\n---------------------------------------------------------------')
    print('PSMs column')
    print(f'Type new column - {type(PSMs_col)}')
    print(f'Shape new column - {np.shape(PSMs_col)}')
    print(f'new column - \n{PSMs_col}')

    # ******* GENERATE PEPTIDE INTENSTIES ******* #
    pep_inten_col = np.random.randint(pep_int_range[0], pep_int_range[1], total_peptides)
    pep_inten_col = np.reshape(pep_inten_col, (-1, 1))  # The unspecified value is inferred to be the number of rows
    print('\n\n---------------------------------------------------------------')
    print('pep_inten_col')
    print(f'Type new column - {type(pep_inten_col)}')
    print(f'Shape new column - {np.shape(pep_inten_col)}')
    print(f'new column - \n{pep_inten_col}')

    print('-' * 10, '\nCONCATENATE PSM - intensities\n')
    PSM_int_cols = np.concatenate((PSMs_col, pep_inten_col), axis=1)
    print(PSM_int_cols)

    # ******* GENERATE PTMs ******* #
    PTM_none = np.repeat('None', total_peptides)
    PTM_none = np.reshape(PTM_none, (-1, 1))  # The unspecified value is inferred to be the number of rows
    print('\n\n---------------------------------------------------------------')
    print('PTM_none')
    print(f'Type new column - {type(PTM_none)}')
    print(f'Shape new column - {np.shape(PTM_none)}')
    print(f'new column - \n{PTM_none}')

    print('-' * 10, '\nCONCATENATE PTM - PSM - intensities\n')
    PTM_PSM_int_cols = np.concatenate((PTM_none, PSM_int_cols), axis=1)
    print(PTM_PSM_int_cols)

    print('-' * 10, '\nCONCATENATE \n')
    self.dummy_input_matrix = np.concatenate((self.dummy_input_matrix, PTM_PSM_int_cols), axis=1)
    print(self.dummy_input_matrix)
    # PSMs_col=np.array([PSMs_col])
    # PSMs_col=PSMs_col.T
    # print('\n\n---------------------------------------------------------------')
    # print('NEW COLUMN')
    # print(f'Type new column - {type(PSMs_col)}')
    # print(f'Shape new column - {np.shape(PSMs_col)}')
    # print(f'new column - \n{PSMs_col}')
    # self.dummy_input_matrix=np.concatenate((self.dummy_input_matrix, PSM_int_cols),axis=1)
    # #self.dummy_input_matrix=np.append(self.dummy_input_matrix, PSMs_col,axis=1)

    # ++++++++++++++++++ EXTRACT THE PEPTIDE PERCENTAGE ++++++++++++++++++ #
    np.random.shuffle(self.dummy_input_matrix)
    peptide_shrinking = (total_peptides * pept_percentage) // 100
    self.dummy_input_matrix = self.dummy_input_matrix[0:peptide_shrinking, :]

    # /////////////////// APPLY THE PTM PERCENTAGE \\\\\\\\\\\\\\\\\\ #
    if PTM_percentage != 0:
        PTM_quantity = (peptide_shrinking * PTM_percentage) // 100
        for pep_ind in range(PTM_quantity):
            pep_row = self.dummy_input_matrix[pep_ind, :]  # Consider the current peptide in the shrinked subset
            PTM_pos = random.randint(0, len(pep_row[1]) - 1)
            PTM_amino_acid = pep_row[1][PTM_pos]
            PTM_type = random.choice(PTMs_types)
            PTM = PTM_type + ' (' + PTM_amino_acid + str(PTM_pos) + ')'  # generate the PTM encoding
            self.dummy_input_matrix[pep_ind, 2] = PTM  # write the PTM encoding in the peptide row

    self.make_sep_file(pep_table_fn, self.dummy_input_matrix, sep='\t')




def prot_not_represented(CDS_ind, pep_ind):
    """
    Version: 1.0

    History Name :

    This function cleans the self.prot_CDS_index class attribute of
    protein codes that are not present in the self.prot_pep_index class attribute.
    One of the reason for the discrepancy could comes from the dummy_input function.
    Here two thresholds have been set to define the length of the
    dummy peptides (max - min).
    This could generate the case that some proteins are not represented in the
    self.prot_pep_index attribute because none of its dummy peptides met the length
    requirements set for the simulation..

    INPUT :
    OUTPUT:
    """
    for protein in list(CDS_ind.keys()): # Consider only the protein IDs
        if protein not in list(list(pep_ind.keys())):
            del CDS_ind[protein]

    # for protein in list(self.prot_CDS_index.keys()): # Consider only the protein IDs
    #     if protein not in list(list(self.prot_pep_index.keys())):
    #         del self.prot_CDS_index[protein]
    return CDS_ind
