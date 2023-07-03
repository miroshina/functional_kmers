#!/usr/bin/env python

'''
Original code by Kamil Slowikowski
(https://gist.github.com/slowkow/a2327b868f4e927ac4fb)

Overwritten and adjusted by Alexandra Miroshina
'''

import urllib.request, urllib.error, urllib.parse
from bs4 import BeautifulSoup as bs
import sys
import pandas as pd
import os

def main():
    args = sys.argv
    if len(args) != 3:
        print('Usage: ./pathway_genes.py pathway organisms_id')
        sys.exit(1)
    path = str(args[1])
    org_dir = str(args[2])
    os.system('mkdir '+path)
    orthology_ids = get_orthology_ids(path)
    print(orthology_ids)
    for ort in orthology_ids:
        main_ko(ort,path,org_dir)


def main_ko(ko,out_dir,organisms_dir):
    out_name = out_dir+'/'+ko
    os.system('mkdir '+out_name)
    os.system('mkdir '+out_name+"/"+ko+"_msa")
    os.system('mkdir '+out_name+'/stat')
    os.system('mkdir '+out_name+'/pfm')
    org_list = os.listdir(organisms_dir)
    for organisms_file in org_list:
        file_name = out_name+'/'+ko + "_" + os.path.basename(organisms_file).split(".")[0]
        file_msa = out_name + "/"+ko+"_msa/" + ko + "_" + os.path.basename(organisms_file).split(".")[0]+"_msa_nt.fa"
        file_cialign = out_name + "/stat/"+ko+'_'+ os.path.basename(organisms_file).split(".")[0][:-4]
        organisms_file_full_path = organisms_dir + organisms_file
        file_name_nt = file_name + "_nt" + ".fa"
        perform(ko, organisms_file_full_path, file_name_nt, "nt")
        os.system('muscle -align ' + file_name_nt + ' -output ' + file_msa)
        os.system('CIAlign --infile '+ file_msa +' --plot_stats_input --outfile_stem '+ file_cialign+' --pwm_input')
        os.system('mv '+out_name+'/stat/*_pfm_input.txt '+out_name+'/pfm/')
        os.system('rm '+out_name+'/stat/*.png')
        os.system('rm '+out_name+'/stat/*.fasta')
        os.system('rm '+out_name+'/stat/*_log.txt')
        os.system('rm '+out_name+'/stat/*_removed.txt')
        os.system('rm '+out_name+'/stat/*_ppm_input.txt')
        os.system('rm '+out_name+'/stat/*_pwm_input.txt')

def perform(ko, organisms_file, file_name, type_of_seq):
    print(organisms_file)
    if ".xlsx" not in ko:
        orthology_ids = [ko]


    df = pd.read_csv(organisms_file, header=None)
    organism_list = df[0].tolist()
    organism_ids = set(organism_list)
    check_list = []
    print(orthology_ids)

    for orthology_id in orthology_ids:
        gene_ids = get_gene_ids(orthology_id).split("\n")
        print('Writing  sequences')
        with open(file_name, 'w') as out:
            for line in gene_ids:
                if len(line.split("\t")) > 1:
                    gene_id = line.split("\t")[1]
                    organism_id = gene_id.split(":")[0]
                    if organism_id in organism_ids and organism_id not in check_list:
                        fasta = get_fasta(gene_id,type_of_seq)
                        out.write(fasta)
                        check_list.append(organism_id)



def get_ids(url):
    response = urllib.request.urlopen(url)
    html = response.read()
    b = bs(html, features="lxml")
    links = b.find_all('a')
    texts = []
    for link in links:
        href = link.get('href')
        # if len(link.text)>0 and not(href=="https://www.genome.jp/" or href=="/dbget/"):
        if href and "entry" in href:
            texts.append(link.text)
    return texts


def get_orthology_ids(pathway_id):
    URL = 'https://www.genome.jp'
    FUN = '/dbget-bin/get_linkdb?-t+orthology+pathway:'
    print("here")
    return get_ids(URL + FUN + pathway_id)


def get_gene_ids(orthology_id):
    URL = 'https://rest.kegg.jp/'
    FUN = '/link/genes/'
    response = urllib.request.urlopen(URL + FUN + orthology_id)
    html = bs(response.read(), features="lxml")
    return html.get_text()


def get_fasta(gene_id, type_of_seq):

    if type_of_seq == "aa":
        URL = 'https://rest.kegg.jp'
        FUN = '/get/'
        FUN2 = '/aaseq'
        try:
            response = urllib.request.urlopen(URL + FUN + gene_id + FUN2)
            html = bs(response.read(), features="lxml")
        except urllib.error.HTTPError:
            print("No such sequence")
            return "\n"
    else:
        URL = 'https://rest.kegg.jp'
        FUN = '/get/'
        FUN2 = '/ntseq'
        try:
            response = urllib.request.urlopen(URL + FUN + gene_id + FUN2)
            html = bs(response.read(), features="lxml")
        except urllib.error.HTTPError:
            print("No such sequence")
            return "\n"

    return html.get_text()




if __name__ == '__main__':
    main()
