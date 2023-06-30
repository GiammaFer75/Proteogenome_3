def file_to_lst(file_name, rem_last_row=False):
    """
    Version : 1.1

    Name History : file2list (from Proteogenome 1.0)

    This function upload in a list the file rows.
    The '\n' character will be removed from each line.

    INPUT  :    file_name       String   The name of the file to upload.
                rem_last_row    Bool     Flag that allows to remove the last row

    OUTPUT :  file_in_lst List     The output list with the content of the file.
    """
    file_in_lst = []

    file_hand = open(file_name, 'r')
    file_in_lst=file_hand.read().splitlines()
    file_hand.close()

    if (rem_last_row == True): file_in_lst = file_in_lst[:-1]  # Remove the last row
    return file_in_lst


def load_generic_table(filename, sep='\t', PoGo=True):
    """
    Version: 2.0

    Name History: load_generic_table

    This function returns a np.array that contains data read from a formatted file.
    There is only one parsing mode readable from this function. This mode is a separated fields
    format. The separation element could be provided by the user (tab separation is the default)

    PoGo WARNING:
    Using PoGo software for mapping the peptides could generate A issues in the output file.
    This problem is related to the chromosome number/name. For instance in bacteria annotations,
    instead of the chromosome number is reported the strain code. Compared to the 22 human chromosomes
    (maximum only 2 chars long) the other organism strain IDs could be more than 10 characters.
    Unfortunately, PoGo cannot manage strain IDs and interprets them as chromosome IDs. As a result,
    the chromosome ID max length will became the strain ID max length (from max 2 to >10 char).

    Example:
    Considering the HCMV Strain Merlin ID = NC_006273.2
    The first three columns of the Peptide Map will be like this:

    Column Position  ----->   0                1            2
                    |     Chromosome     |   Start    |    End     | ..........
                    |       Number       | Coordinate | Coordinate | ..........
                    |--------------------|------------|------------| ..........
                    |   chrN C_006273.2  |  122714    |   122792   | ..........
                            ^
                            |
                        This additional space will affect the map view in the genome browser.

    If the flag PoGo is true the function will remove the space characters in the first table column.

    INPUT :
    OUTPUT:
    """
    fh = open(filename, 'r')
    first_line = fh.readline().rstrip()  # Divide in columns the first row
    number_of_columns = len(first_line.split(sep))  # Fetch the number of columns from the first row

    fh.close()

    input_tab = np.empty((1, number_of_columns))

    fh = open(filename, 'r')
    for row in fh:
        row = row.split(sep)
        row = np.array(row)
        row[-1] = row[-1].replace('\n', '')
        input_tab = np.concatenate([input_tab, [row]], axis=0)
    fh.close()
    input_tab = input_tab[1:, :]  # Remove initialisation row

    # If there is a space into the chromosome name field,
    # therefore PoGo has not been able to manage the name length.
    if (' ' in input_tab[0, 0]) & (PoGo == True):
        for ind, col_0_row in enumerate(input_tab[:, 0]):
            print(col_0_row)
            input_tab[ind, 0] = str(col_0_row).replace(' ', '')

    return input_tab

def load_input_table(filename, sep='\t'):
        """
        Version: 1.0

        Name History: load_input_table

        This function upload in a np.array the PROTEOMICS DATA that must be visualised
        on the genome browser.
        The format of this table is:

                    Protein ID | Peptide Sequence | PTM | PSM | Peptide Intensity

        INPUT :
        OUTPUT:
        """
        # import numpy as np
        input_tab = np.array([['', '', '', '', '']], dtype='object')

        fh = open(filename, 'r')
        for row in fh:
            row = row.split(sep)
            row = np.array(row[:-1])
            # row[-1]=row[-1].replace('\n','')
            input_tab = np.concatenate([input_tab, [row]], axis=0)
        fh.close()
        input_tab = input_tab[1:, :]  # Remove initialisation row
        return input_tab


def PoGo_input_table(experiment_tag='exp', out_file_name=''):
    """
    Version: 1.0

    Name History: PoGo_input_table

    This function creates a peptide input table suitable for the PoGo software.
    Considering the format of the input table for Proteogenome,
    this function removes the first column of the original input table.
    Then, insert a new column into the first position with an experiment tag
    that will be reported on each peptide row.

    INPUT : self.imput_table
    OUTPUT:
    """
    self.PoGo_it = np.delete(self.input_table, 0, 1)  # Remove the protein codes column. self.input_table[:,1::]

    PTMs = self.PoGo_it[:, 1]  # Save the PTMs type and location for each peptide in an array.
    print(self.PoGo_it)
    # peptides_updated=self.apply_PTMs_to_pep(peptide_tab=self.PoGo_it) # Apply the PTM to the peptide sequence #
    self.PoGo_it = self.apply_PTMs_to_pep(peptide_tab=self.PoGo_it)

    self.PoGo_it = np.delete(self.PoGo_it, 1, 1)  # Remove the PTM column
    print(self.PoGo_it)

    insert_experiment_name = np.full((len(self.PoGo_it), 1), experiment_tag)  # Generate the experiment tag column.
    self.PoGo_it = np.concatenate((insert_experiment_name, self.PoGo_it),
                                  axis=1)  # Insert the experiment tag column as a first column of the table.

    if out_file_name:
        print('making PoGo input table')
        self.make_sep_file(out_file_name, self.PoGo_it, sep='\t')  # Create the file with the the PoGo input table.
