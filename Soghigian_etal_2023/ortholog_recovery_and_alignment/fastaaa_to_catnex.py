#!/usr/bin/python
from Bio import AlignIO
from Bio.Nexus import Nexus
import os
import glob
#Convert all alignments to nexus format.
for fn_nt in os.listdir("."):
    if fn_nt.endswith(".fasta"):
        input_handle_nt = open(os.path.join('.',fn_nt), "r")
        fn2_nt = fn_nt.replace(".fasta", "") + ".nex"
        output_handle_nt = open(os.path.join('.', fn2_nt), "w")
        alignments_nt = AlignIO.read(input_handle_nt, "fasta")
        AlignIO.convert(fn_nt, "fasta", output_handle_nt, "nexus", molecule_type="protein")
        output_handle_nt.close()
        input_handle_nt.close()

#Grab all nexus files.
file_list_nt = glob.glob(os.path.join('.','*.nex'))
nexi_nt =  [(fname_nt, Nexus.Nexus(fname_nt)) for fname_nt in file_list_nt]
combined_nt = Nexus.combine(nexi_nt)
combined_nt.export_phylip(filename='all_orthologs.phy')
combined_nt.write_nexus_data(filename='all_orthologs.nex')
