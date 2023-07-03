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


def main():
    args = sys.argv
    if len(args) != 4:
        print('Usage: ./pathway_genes.py ORTHOLOGY organisms_id')
        sys.exit(1)

    first_arg = str(args[1])
    organisms_file = str(args[2])
    file_name = str(args[3])+".fa"
    print(args[1])
    print(args[2])
    if ".xlsx" not in first_arg:
        orthology_ids = [str(args[1])]
    else:
        df = pd.read_excel(first_arg, header=None)
        orthology_list = df[0].tolist()
        orthology_ids = set(orthology_list)


    df = pd.read_excel(organisms_file, header=None)
    organism_list = df[0].tolist()
    organism_ids = set(organism_list)

    for orthology_id in orthology_ids:
        gene_ids = get_gene_ids(orthology_id).split("\n")
        print('Writing {} FASTA gene sequences to "{}.fa"' \
              .format(len(gene_ids), orthology_id))

        with open(file_name, 'w') as out:
            for line in gene_ids:
                if len(line.split("\t")) > 1:
                    gene_id = line.split("\t")[1]
                    organism_id = gene_id.split(":")[0]
                    if organism_id in organism_ids:
                        fasta = get_fasta(gene_id)
                        out.write(fasta)


def get_ids(url):
    response = urllib.request.urlopen(url)
    html = response.read()
    b = bs(html, features="lxml")
    links = b.find_all('a')
    texts = []
    for link in links:
        href = link.get('href')
        #if len(link.text)>0 and not(href=="https://www.genome.jp/" or href=="/dbget/"):
        if href and "entry" in href:
            texts.append(link.text)
    return texts

def get_gene_ids(orthology_id):
    URL = 'https://rest.kegg.jp/'
    FUN = '/link/genes/'
    response = urllib.request.urlopen(URL + FUN + orthology_id)
    html = bs(response.read(), features="lxml")
    return html.get_text()


def get_fasta(gene_id):

    URL = 'https://rest.kegg.jp'
    FUN = '/get/'
    FUN2 = '/aaseq'
    try:
        response = urllib.request.urlopen(URL + FUN + gene_id + FUN2)
        html = bs(response.read(), features="lxml")
    except urllib.error.HTTPError:
        print("No such sequence")
        return "\n"

    return html.get_text()


if __name__ == '__main__':
    main()
