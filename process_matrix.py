import sys
import os

# directory with all matrices
mat_dir = sys.argv[1]
# file with enzyme names and KEGG ids
enzymes = sys.argv[2]
enz_dict = {}
# file name for the output matrix (tsv)
out_file = sys.argv[3]

# read enzyme file and save all info into hash map
with open(enzymes,"r") as ff:
    for line in ff:
        ko = line.split("\t")[0].rstrip()
        name = line.split("\t")[1].rstrip()
        enz_dict[ko]=name

mat_dirList = os.listdir(mat_dir)
distinct_paths = set()
distinct_kmers = set()

# write info from all matrices into one big matrix, but include only those sequences that have blast hits to one certain organism more than 90% 
# of the time. Additionaly count number of distinct kmers and distinct pathways/KOs
with open(out_file,"w") as file:
    for mat in mat_dirList:
        with open(mat_dir+"/"+mat, "r") as f:
            for line in f:
                if "KO_path" not in line:
                    columns = line.split("\t")
                    path = columns[1]
                    enz_name = enz_dict[path].rstrip()
                    kmer = columns[2]
                    organism = columns[3]
                    tupel = columns[5]
                    hit = int(tupel.split(",")[0][1:])
                    all_hits = int(columns[6])
                    if all_hits == hit and hit!=0:
                        distinct_paths.add(path)
                        distinct_kmers.add(kmer)
                        file.write(enz_name + "\t" + organism + "\t" + kmer + "\t" + str(hit) + "\t" + "100%\n")
                    elif hit >= all_hits*0.9 and hit !=0:
                        distinct_paths.add(path)
                        file.write(enz_name + "\t" + organism + "\t" + kmer + "\t" + str(hit) + "\t" + ">90%: "+str(hit)+" from "+str(all_hits)+"\n")
                        distinct_kmers.add(kmer)
    file.write("Distinct pathways: "+str(len(distinct_paths))+"\t---"+"\t"+"Number of kmers: "+str(len(distinct_kmers))+"\t---"+"\t"+"\t---"+"\t---"+"\n")
