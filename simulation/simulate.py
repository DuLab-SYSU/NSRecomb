import os
import sys
from pathlib import Path
import subprocess
from Bio import SeqIO

# add package to python module path
p = Path.cwd()
root_dir = p.parent
bin_dir = root_dir.joinpath("source")
sys.path.append(str(bin_dir))

import Prepare
import BlastTools
import Verify
import Iteration

# specify dependecncy path
Verify.KAKS_PATH = root_dir.joinpath('package', 'KaKs_Calculator2.0', 'bin')
Verify.MAFFT_PATH = root_dir.joinpath('package', 'mafft-linux64')
BlastTools.BLAST_PATH = root_dir.joinpath('package', 'ncbi-blast-2.11.0+', 'bin')
Prepare.BLAST_PATH = root_dir.joinpath('package', 'ncbi-blast-2.11.0+', 'bin')


simu_dir = root_dir.joinpath("simulation", "replicates")
replicate_dirs = [dir_ for dir_ in simu_dir.iterdir() if dir_.is_dir()]
for replicate_dir_i in replicate_dirs:
    replicate_name = replicate_dir_i.name
    data_names = replicate_dir_i.joinpath("Results").glob("sequences*")
    for data_name_i in data_names:
        dataset_path = data_name_i
        data_name = dataset_path.stem
        
        # build local blast database for simulated phyloengy
        DATABASE = root_dir.joinpath('db', 'simu')
        Prepare.build_db(dataset_path, DATABASE)
        
        dataset = list(SeqIO.parse(dataset_path, "fasta"))
        # identify recombination event for each sequences
        for i, QUERY in enumerate(dataset):
            # exclude the query sequence from the database
            mask_content = QUERY.id
            mask_in_path = simu_dir.joinpath("mask.acc")
            mask_out_path = simu_dir.joinpath("mask.id")
            
            with open(mask_in_path, "w") as f:
                f.write(mask_content)
                
            Prepare.convert_to_binary(mask_in_path, mask_out_path)
            NEGATIVE_SEQIDLIST = mask_out_path
            
            BACKBONE = BlastTools.get_gene_backbone(QUERY, DATABASE, mask_out_path)
            OUTPUT = simu_dir.joinpath("simu_%s_%s_%s.result" % (replicate_name, data_name, i))

            NUM_CPUs = 12
            WINDOW_SIZE = 501
            STEP = 30
            CODON = False

            Iteration.main(OUTPUT, QUERY, BACKBONE, DATABASE, NEGATIVE_SEQIDLIST, NUM_CPUs, WINDOW_SIZE, STEP, codon_bootstrap=CODON)
