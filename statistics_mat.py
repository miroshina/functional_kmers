
import sys
import os
import collections
from collections import defaultdict

map_dir = sys.argv[1]
mat_dir = sys.argv[2]

path_to_matrix = sys.argv[3]
path_to_names = sys.argv[4]

map_dirList = os.listdir(map_dir)

pathway_d = defaultdict(lambda:defaultdict(int))

ko_map_dict = defaultdict(str)


for map in map_dirList:
    full_path = map_dir+map+"/"
    path_dirList = os.listdir(full_path)
    for path in path_dirList:
        msa_path = map_dir+map+"/"+path+"/"+path+"_msa/"
        organism_count = len( os.listdir(msa_path) )
        pathway_d[map][path] = organism_count
        ko_map_dict[path] = map


for map,paths in pathway_d.items():
    print("{}\t{}".format(map,len(paths)))
    for path,count in paths.items():
        print("{}\t{}".format(path,count))
    print()


print()

empty = {}

pathways_all = 0
pathways_nonempty_before_counting = 0
file1 = open("map_info.tsv","w")
for map, y in pathway_d.items():
    zeros = collections.Counter(y.values())[0]
    empty[map]=zeros
    pathways_all+=len(y)
    pathways_nonempty_before_counting+=(len(y)-zeros)
    file1.write("{}\t{}\t{}\n".format(map,len(y),len(y)-zeros))
    print("Number of kos in map {} is {} and {} kos are empty".format(map,len(y),zeros))
file1.write("---\t---\t---\n")
print()
print("Insgesamt gibt es {} davon sind {} nicht-leer".format(pathways_all,pathways_nonempty_before_counting))
print()
after_blasting = defaultdict(lambda:defaultdict(lambda: defaultdict(int)))

mat_dirList = os.listdir(mat_dir)
conserved = defaultdict(lambda:defaultdict(int))
for name in mat_dirList:
    map = name.split("_matrix")[0]
    full_path = mat_dir + name
    with open(full_path,"r") as file:
        for line in file:
            if "KO_path" not in line:
                list = line.split("\t")
                path = list[1]
                organism = list[3]
                hit = int(list[5][1:].split(",")[0])
                all_hit = int(list[6])
                if hit == all_hit and hit!=0:
                    after_blasting[map][path][organism]+=1
                conserved[map][path]+=1



file2 = open("ko_info.tsv","w")
for key, value in conserved.items():
    print("Map - {}, number of KOs is {}".format(key,len(value)))
    print()
    file1.write("{}\t {}\t---\n".format(key,len(value)))
file1.write("---\t---\t---\n")





for key,value in after_blasting.items():
    print()
    file1.write("Map - {}\tKOs - {}\n".format(key,len(value)))
    print()
    file2.write(key+"\t"+"---\n")
    for key1, value1 in value.items():
        file2.write("{}\t{}\n".format(key1,len(value1)))

file1.close()
file2.close()

