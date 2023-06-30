#!usr/bin/python3
from absl import app
from absl import flags

import subprocess
from os import chdir

from Test_PoGo_absl_subprocess import printinput as pi

# ****************************************************************************** #
#                               PATHS DEFINITION                                 #
home_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/'
pogo_dir = 'C:/Bioinformatics/Proteogenome_v3/Test_PoGo_absl_subprocess/Proteogenome-master/Proteogenome-master/PoGo/Windows'
# ****************************************************************************** #


FLAGS = flags.FLAGS

flags.DEFINE_string('pogo_windows_exe', None, 'The executable PoGo file for running the software on Windows OS.')
flags.DEFINE_string('protein_FASTA_seq', None, 'The protein sequences related to the reference organism genome.')
flags.DEFINE_string('protein_GTF_annots', None, 'The protein annotations related to the reference organism genome.')
flags.DEFINE_string('peptides_table', None, 'The table containing the peptides sequences with the format:\n'
                                            '- Protein ID\n'
                                            '- Peptide sequence\n'
                                            '- Peptide modification - PTM\n'
                                            '- Peptide PSM\n'
                                            '- Peptide intensity\n')

def main(argv):

    # ***** FLAGS DEFINITION ***** #
    pogo_windows_exe_path = FLAGS.pogo_windows_exe
    protein_FASTA_seq_path = FLAGS.protein_FASTA_seq
    protein_GTF_annots_path = FLAGS.protein_GTF_annots
    peptides_table_path = FLAGS.peptides_table

    # Upload files and create indexes


    #TESTING COMMANDS
    pi.print_input(pogo_windows_exe_path, protein_FASTA_seq_path, protein_GTF_annots_path, peptides_table_path)
    chdir(pogo_dir)
    subprocess.run(['.\PoGo.exe', 'capture_output=True'])
    subprocess.run(['dir'], shell=True)




if __name__ == '__main__':
    """
    Proteogenome
    
    Version: 3.0.0
    
    Author: Giammarco Ferrari
    """
    app.run(main)