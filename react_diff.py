import sys
import re
from Bio import SeqIO
import math

mat_path = r"C:\Users\kunde\PycharmProjects\bachelor\matrix_with_pairs.tsv"
enz_path = r"C:\Users\kunde\Desktop\enzyme_names_tmp.tsv"
enz_map = {}

# function to extract nucleotide sequences from fasta files into a list
def get_sequences(msa):
    result = []
    for fasta in msa:
        sequence = str(fasta.seq)
        result.append(sequence.upper())
    return result

# map enzyme names back to their KEGG ids
with open(enz_path,"r") as file:
    for line in file:
        spli = line.split("\t")
        enz = spli[0].rstrip()
        name = spli[1].rstrip()
        enz_map[name]=enz

# for each pair, find it in the sequences in the corresponding file and look if they overlap 
with open(mat_path,"r") as file:
    kmer_map = {}
    kmer_count = {}
    remove = {}
    for line in file:
        att = line.split("\t")
        if len(att)>3:
            # all attributes from the row
            path = att[0].rstrip()
            org = att[1].rstrip()
            kmer1 = att[2].rstrip()
            kmer2 = att[3].rstrip()
            # map enzyme name back to kegg id
            enz = enz_map[path]
            count = 0
            # build a path to the corresponding file
            file_name ="C:\\Users\\kunde\\Desktop\\KOs\\"+enz+"_"+org+"_ids_nt.fa"
            try:
                msa = SeqIO.parse(open(file_name), 'fasta')
                alignment = get_sequences(msa)
                ind = []
                for seq in alignment:
                    a = -1
                    b = -1
                    for i in re.finditer(kmer1, seq):
                        a = i.start()
                    for j in re.finditer(kmer2, seq):
                        b = j.start()
                    if a != -1 and b != -1:
                        # if any kmer pairs lies in a distance smaller than 23, than we should remove it from our table
                        if abs(b - a) < 23 and b - a > 0:
                            remove[(kmer1, kmer2)] = (path, org)
                        elif abs(b - a) >= 23:
                            count+=1
                            ind.append(abs(b - a))
                # save the average distance for the current pair
                if len(ind) > 0:
                    kmer_map[(path,org,kmer1, kmer2)] = sum(ind) / len(ind)
                    kmer_count[(kmer1,kmer2)] = count
            except FileNotFoundError:
                print("no file for {}".format(enz))

# write all info into a table
with open("remain.tsv","w") as file:
    for key in kmer_map.keys():
        file.write("{}\t{}\t{}\t{}\t{}\n".format(key[0],key[1],key[2],key[3],kmer_map[key]))






