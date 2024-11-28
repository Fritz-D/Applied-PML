#
# Small sequence selector and pickle creation


from Bio import SeqIO
import numpy as np
import pandas as pd
import pickle

DNA = {"A": 0 , "C": 1, "G":2, "T": 3, 
       "a":0, "c":1, "g":2, "T": 3}

def seq_as_list(strseq):
    """
    str -> list[int]
    default is A
    """
    return [DNA.get(c, 0) for c in strseq]

seqin = "Makona_1610_genomes_2016-06-23.fasta"

s=6050
e=8080

ebola_seqs = SeqIO.parse(open(seqin),'fasta')

all_seqs = []
all_infos = []

for fasta in ebola_seqs:
    name, seq = fasta.id, str(fasta.seq)
    seq_keep = seq_as_list(seq[s:e])
    all_seqs.append(seq_keep)
    all_infos.append(name.split("|"))

infos = pd.DataFrame(all_infos, 
                    columns = ["virus", "ID", "accession", "country", "city", "date"])

seqs = np.array(all_seqs)
#selecting positions that change

infos["date"] = pd.to_datetime(infos["date"])
filtered_infos = infos["date"] < "2015-01-01"
seqs_filt = seqs[filtered_infos, :]

posequal = np.all(np.equal(seqs_filt, seqs_filt[0,:]), axis = 0)
seqs_filt_trim = seqs_filt[:, np.logical_not(posequal)]

with open('ebola.pkl', 'wb') as f:
    pickle.dump(seqs_filt_trim[: ,200:800],f)

filt_seqs_infos = infos[filtered_infos]
filt_seqs_infos.to_csv("ebola_sequences_information.csv")
