from Bio import SeqIO
import pandas as pd
import math
import sys
import os
from collections import OrderedDict
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from datetime import datetime
import xml.etree.ElementTree as ET
import pandas as pd



def main(alignment_dir, statistics_dir, pfm_dir, KO_name):

    df_ko = pd.DataFrame()
    df_ls = []
    kos_dict = {}
    alg_dirList = os.listdir(alignment_dir)
    stat_dirList = os.listdir(statistics_dir)
    pfm_dirList = os.listdir(pfm_dir)

    for i in alg_dirList:
        ko = os.path.basename(i).split("_")[1].split(".")[0]
        kos_dict[ko] = [alignment_dir + i]
    for j in stat_dirList:
        ko = os.path.basename(j).split("_")[1]
        kos_dict[ko].append(statistics_dir + j)
    for k in pfm_dirList:
        ko = os.path.basename(k).split("_")[1]
        kos_dict[ko].append(pfm_dir + k)

    for key in kos_dict:
        if len(kos_dict[key]) == 3 and "Burkhol" not in key:
#            print(key,KO_name)
            path_to_alignment, path_to_stat, path_to_pfm = kos_dict[key]
            stat_data = pd.read_csv(path_to_stat, sep='\t', index_col=0)
            pfm_data = pd.read_csv(path_to_pfm, sep='\t', index_col=0)
            number_of_seq = len([rec for rec in SeqIO.parse(path_to_alignment, "fasta")])
            seq_len = len(pfm_data.columns)
            msa = SeqIO.parse(open(path_to_alignment), 'fasta')
            alignment = get_alignment(msa)
            check = perform_gap_processing(number_of_seq, pfm_data, seq_len, path_to_alignment, 23,key,KO_name)
            if type(check)!=int:
                df_ls.append(check)
                print(key)
    if len(df_ls) > 0:
        df_ko = pd.concat(df_ls)
        return df_ko
    return -1

def run_blast():
    #with open("kmer.fasta","w") as file:
    #    file.write(">"+organism+"_1\n")
    #    file.write(seq)
    #file.close()
    os.system("blastn -query kmer.fasta -db /nfs/data/functional_k_mers/humGut/blast_db/HumGut975_library.fna -task blastn-short -num_alignments 500 -qcov_hsp_perc 100 -perc_identity 100 -out result.xml -outfmt 5")
    return process_xml()


def process_xml():
    taxonomy = pd.read_csv("HumGut.tsv", sep='\t', index_col=None)
    tree = ET.iterparse("result.xml")
    result = {}
    count_blast=[]
    count_blast_all=[]
    for event, elem in tree:
        if elem.tag == 'Iteration':
            id_full = elem.find("Iteration_query-def")
            sp = str(id_full.text).split("_")[0]
            path = str(id_full.text).split("_")[1]
            ls = elem.find("Iteration_hits")
            species_ls = []
            count = 0
            for element in ls:
                hit = element.find('Hit_def')
                id = hit.text.split("|")[1]
                species = str(taxonomy.loc[taxonomy["HumGut_tax_id"] == int(id), "gtdbtk_organism_name"])
                if sp in species:
                    count += 1
                species_ls.append(species)
            result[str(id_full.text)] = [path, sp, str(count) + " from " + str(len(species_ls))]
            count_blast.append((count,sp))
            count_blast_all.append(str(len(species_ls)))
    return count_blast,count_blast_all



def get_alignment(msa):
    result = []
    for fasta in msa:
        sequence = str(fasta.seq)
        result.append(sequence)
    return result


def perform_gap_processing(number_of_seq, pfm_data, seq_len, path_to_alignment, k_size, organism,KO):
    if number_of_seq < 15:
        return -1
    #if "Burkhol" in organism:
    #    return -1
    gaps_positions = extract_gaps_indices(number_of_seq, pfm_data, seq_len)
    msa = SeqIO.parse(open(path_to_alignment), 'fasta')
    alignment = []
    for fasta in msa:
        sequence = str(fasta.seq)
        alignment.append(sequence)
    if len(gaps_positions) == 0:
        sequences, pfm_data_new = alignment, pfm_data
    else:
        sequences, pfm_data_new = remove_gaps(gaps_positions, pfm_data, alignment)
    seq_len = len(sequences[0])
    conservation_positions = extract_sequence_indices(k_size, number_of_seq, pfm_data_new, seq_len)

    map_of_kmers = {}

    for indices in conservation_positions:
        start = indices[0]
        end = indices[1]
        count_kmers_variability = {}
        for sequence in sequences:
            kmer = sequence[start:end].upper()
            if kmer not in count_kmers_variability.keys():
                count_kmers_variability[kmer] = 1
            else:
                count_kmers_variability[kmer] += 1

        max_value = max(count_kmers_variability.values())
        kmers_candidates_atposition = [(k, v, number_of_seq) for k, v in count_kmers_variability.items() if
                                       (v == max_value) and (v >= math.ceil(number_of_seq / 2))]
        if len(kmers_candidates_atposition) > 0:
            map_of_kmers[(start, end)] = kmers_candidates_atposition

    d_descending = OrderedDict(sorted(map_of_kmers.items(), key=lambda kv: kv[1][0][1], reverse=True))

    column_kmer = []
    column_organims = []
    column_count = []
    column_path = []
    column_blast = []
    column_blast_hits = []
    count = 0
    if len(d_descending) > 0:
        with open("kmer.fasta","w") as file:
            for key, value in d_descending.items():
                for kmer in value:
                    column_path.append(KO)
                    column_kmer.append(kmer[0])
                    column_count.append(str(kmer[1]) + " from " + str(kmer[2]))
                    column_organims.append(organism)
                    file.write(">"+organism+"_"+KO+str(count)+"\n")
                    file.write(kmer[0]+"\n")
                    count+=1
#                blast_org,blast_all = run_blast(kmer[0],organism)
#                column_blast.append(blast_org)
#                column_blast_hits.append(blast_all)
        column_blast, column_blast_hits = run_blast()
        dt = {"KO_path": column_path,
                  "Kmer": column_kmer,
                  "Organism": column_organims,
                  "Count": column_count,
                  "Blast_result":column_blast,
                  "Blast_number_of_hits":column_blast_hits
            }
        df = pd.DataFrame(dt)
        return df

    return -1



def extract_sequence_indices(len_threshold, number_of_sequences, data_df, seq_length):
    i = 0
    gap = True
    counter = 0
    result = []
    index_saver = -1
    skipped_pos = 0
    while i < seq_length:
        k = data_df[str(i)]
        if (not (k > number_of_sequences * 0.6).any()) and counter == 0:
            i += 1
            gap = True
            counter = 0
            continue
        elif (not (k > number_of_sequences * 0.6).any()) and counter >= len_threshold:
            result.append((index_saver, i))
            i += 1
            gap = True
            counter = 0
            continue
        elif (not (k > number_of_sequences * 0.6).any()) and (0 < counter < len_threshold):
            i += 1
            index_saver = -1
            gap = True
            counter = 0
            continue
        elif ((k > number_of_sequences * 0.6).any()):
            if gap:
                counter = 1
                index_saver = i
                gap = False
                i += 1
            else:
                counter += 1
                i += 1

    kmer_indices = []

    for indices in result:
        start = indices[0]
        end = indices[1]
        if end - start > len_threshold:
            for i in range(start, end - len_threshold + 1):
                kmer_indices.append((i, i + 23))
        else:
            kmer_indices.append((start, end))

    return kmer_indices


def extract_gaps_indices(number_of_sequences, data_df, seq_length):
    i = 0
    gap = False
    saver = 0

    result = []
    while i < seq_length:
        if data_df[str(i)].sum() <= number_of_sequences * 0.2 and not gap:
            saver = i
            i += 1
            gap = True
        elif data_df[str(i)].sum() <= number_of_sequences * 0.2 and gap:
            i += 1
            gap = True

        elif data_df[str(i)].sum() > number_of_sequences * 0.2 and gap:
            result.append((saver, i))
            i += 1
            gap = False
        elif data_df[str(i)].sum() > number_of_sequences * 0.2 and not gap:
            i += 1

    if len(result) > 0:
        if result[0][0] == 0:
            result.pop(0)

    return result


def remove_gaps(gapspositions, pfm_data, alignment):
    sequences_with_nogaps = []
    new_conservation_data = []

    a = pfm_data.loc["A"].tolist()
    c = pfm_data.loc["C"].tolist()
    g = pfm_data.loc["G"].tolist()
    t = pfm_data.loc["T"].tolist()

    new_a = []
    new_c = []
    new_t = []
    new_g = []

    num_gaps = len(gapspositions)
    for sequence in alignment:
        result_sequence = []
        if num_gaps == 1:
            start_gap = gapspositions[0][0]
            end_gap = gapspositions[0][1]
            result_sequence.append(sequence[0:start_gap])
            result_sequence.append(sequence[end_gap:])
        else:
            for i in range(0, num_gaps):
                start_gap = gapspositions[i][0]
                end_gap = gapspositions[i][1]
                if i == 0:
                    result_sequence.append(sequence[0:start_gap])
                elif 0 < i < num_gaps - 1:
                    end_gap_prev = gapspositions[i - 1][1]
                    result_sequence.append(sequence[end_gap_prev:start_gap])
                else:
                    end_gap_prev = gapspositions[i - 1][1]
                    result_sequence.append(sequence[end_gap_prev:start_gap])
                    result_sequence.append(sequence[end_gap:])

        sequences_with_nogaps.append("".join(result_sequence))

    if num_gaps == 1:
        new_a = a[0:start_gap] + a[end_gap:]
        new_g = g[0:start_gap] + g[end_gap:]
        new_t = t[0:start_gap] + t[end_gap:]
        new_c = c[0:start_gap] + c[end_gap:]


    else:
        for i in range(0, num_gaps):
            start_gap = gapspositions[i][0]
            end_gap = gapspositions[i][1]
            if i == 0:
                new_a = new_a + a[0:start_gap]
                new_g = new_g + g[0:start_gap]
                new_t = new_t + t[0:start_gap]
                new_c = new_c + c[0:start_gap]

            elif 0 < i < num_gaps - 1:
                end_gap_prev = gapspositions[i - 1][1]
                new_a = new_a + a[end_gap_prev:start_gap]
                new_g = new_g + g[end_gap_prev:start_gap]
                new_t = new_t + t[end_gap_prev:start_gap]
                new_c = new_c + c[end_gap_prev:start_gap]


            else:
                end_gap_prev = gapspositions[i - 1][1]
                new_a = new_a + a[end_gap_prev:start_gap]
                new_g = new_g + g[end_gap_prev:start_gap]
                new_t = new_t + t[end_gap_prev:start_gap]
                new_c = new_c + c[end_gap_prev:start_gap]

                new_a = new_a + a[end_gap:]
                new_g = new_g + g[end_gap:]
                new_t = new_t + t[end_gap:]
                new_c = new_c + c[end_gap:]

    new_a.insert(0, "A")
    new_g.insert(0, "G")
    new_c.insert(0, "C")
    new_t.insert(0, "T")

    new_pfm = pd.DataFrame([new_a, new_c, new_g, new_t])
    new_pfm = new_pfm.set_index(0)

    new_index_data = [str(i) for i in range(0, len(sequences_with_nogaps[0]))]
    new_pfm.columns = new_index_data

    return sequences_with_nogaps, new_pfm


if __name__ == '__main__':

    KOs_dir = sys.argv[1]
    KOs_dirList = os.listdir(KOs_dir)
    df_ls = []
    for ko in KOs_dirList:
        align_dir = KOs_dir + ko + "/" + ko + "_msa/"
        stat_dir = KOs_dir + ko + "/" + "stat/"
        pfm_dir = KOs_dir + ko + "/" + "pfm/"
        check = main(align_dir, stat_dir, pfm_dir, ko)
        if type(check) != int:
            df_ls.append(check)
    if len(df_ls)>0:
        df = pd.concat(df_ls)
        print(df)
        df.to_csv(KOs_dir[:-1]+"_matrix.tsv",sep="\t")
