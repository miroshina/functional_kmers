#! /usr/bin/python

from datetime import datetime
from sklearn import preprocessing
import Bio
from Bio import SeqIO, SeqUtils,bgzf
from Bio.SeqUtils import lcc
import pandas as pd
import os
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.metrics import mean_squared_error, r2_score


class Kmers(object):

    def __init__(self, pathToSamples,pathToTest,kmerLength,disease,path_to_kmerset):

        self.pathabundanceDir = ""
        self.testabundaceDir = ""
        self.path_to_kmerset = path_to_kmerset
        self.kmerLength = kmerLength
        self.disease = disease
        self.pathwayMatrix = pd.DataFrame()
        self.pathToSamples = pathToSamples
        self.kmerMatrix = pd.DataFrame()
        self.pathToTest = pathToTest
        self.pathwayTest = pd.DataFrame()
        self.kmerTest = pd.DataFrame()
        self.nsamples=0
        self.ntest=0
        self.model_map={}
        self.kmerset = set()

    def count_kmers(self, listOfFiles):
        sample = 0
        kmerDB = {}
        new = pd.DataFrame()
        countFiles=[]
        self.nsamples= len(listOfFiles)
        print(len(self.kmerset))
        for filepath in sorted(listOfFiles):
            kmerCountFile = 'kmers_in_sample_' + str(sample)
            os.system('kmc -cs1000000000 -k' + str(self.kmerLength) + ' ' + filepath + ' ' + kmerCountFile + ' . >>kmc.log')
            kmerTextFile = kmerCountFile + '.txt'
            os.system('kmc_tools transform ' + kmerCountFile + ' dump ' + kmerTextFile)
            os.system('rm ./*_suf')
            os.system('rm ./*_pre')
            kmPath = './' + kmerTextFile
            countFiles.append(kmPath)
            with open(kmPath, "r") as tmpFile:
                print(kmPath)
                for line in tmpFile:
                    kmer = line.split("\t")[0].rstrip()
                    count = line.split("\t")[1].rstrip()
#complexity=Bio.SeqUtils.lcc.lcc_simp(kmer)
                    if kmer in self.kmerset: #complexity >= 1.6:

                        if kmer not in kmerDB.keys():
                            counts = [0] * self.nsamples
                            counts[sample] = count
                            kmerDB[kmer] = counts
                        else:
                            kmerDB[kmer][sample] = count
                        #if kmer not in kmerDB.keys():
                        #    kmerDB[kmer] = 1
                        #else:
                        #    count=kmerDB[kmer]
                        #    kmerDB[kmer] = count + 1
            sample+=1
            tmpFile.close()
            print(str(sample)+" processed")
            os.system("rm "+kmPath)
        #resultDict = {}
#        for key in kmerDB.keys():
#            value = kmerDB[key]
#            if value >= self.nsamples*80/100 :
#                resultDict[key]=1
        #kmerDB.clear()
        #kmerDB={}
        #sample=0
        #for filepath in countFiles:
        #    print(filepath)
        #    with open(filepath, "r") as tmpFile:
        #        for line in tmpFile:
        #            kmer = line.split("\t")[0].rstrip()
        #            count = line.split("\t")[1].rstrip()
        #            if kmer in kmer.kmerset:
        #                if kmer not in kmerDB.keys():
        #                    counts = [0] * self.nsamples
        #                    counts[sample] = count
        #                    kmerDB[kmer] = counts
        #                else:
        #                    kmerDB[kmer][sample] = count
        #    os.system('rm '+ filepath)
        #    sample+=1

        for i in range(0,len(listOfFiles)):
            tmp=listOfFiles[i].split("/")
            listOfFiles[i]=tmp[-1].split(".")[0]

        pf1 = pd.DataFrame.from_dict(kmerDB, columns=sorted(listOfFiles), orient='index')
        ind=pf1.index.tolist()
        pf1.insert(0, "", [i for i in range(0,len(kmerDB))], True)
        pf1.insert(0, "kmer", ind, True)
        pf1=pf1.set_index(pf1[""])
        del pf1[""]
        pf1=pf1.set_index(pf1["kmer"])
        del pf1["kmer"]
        print(pf1)
        self.kmerMatrix=pf1
        pf1.to_csv('countmatrix_train_1206.tsv', sep="\t")
        print("----- done counting -----\n")




    def read_pathcoverage(self,path,fileName):
        print("----- start reading path abundance files-----")
#        os.system('rm -r '+path+ '/*_temp')
#        os.system('rm '+path+ '/*_genefamilies.tsv')
#        os.system('rm '+path+ '/*_pathcoverage.tsv')
        os.system('rm '+path+'/*_genefamilies_cpm.tsv')
        os.system('rm '+path+'/*_level4ec.tsv')
        os.system('rm '+path+'/*.log')
        dirList = os.listdir(path)
        sample = 0
        count = 1
        new = pd.DataFrame()
        for inputfile in sorted(dirList):
            sample += 1
            ids = []
            func = {}
            name = path + "/" + os.path.basename(inputfile)
            with open(name, "r") as file:
                for line in file:
                    if not "PWY" in line:
                        continue
                    else:
                        splitted = line.split("\t")
                        firstIndex = line.index(":")
                        function = str(splitted[0][firstIndex + 2:]).rstrip()
                        abundance = splitted[1].rstrip()
                        id = splitted[0][:firstIndex]
                        if id not in ids:
                            func[function] = round(float(abundance), 2)
                            ids.append(id)
                tmp = str(os.path.basename(name)).split("_")[0]
                pf1 = pd.DataFrame(func.items(), columns=["function", tmp])
                if sample == 1:
                    new = pd.concat([new, pf1], axis=1)
                else:
                    new = new.merge(pf1, on='function', how="outer").fillna(0)
            print(str(count) + " file done")
            count += 1
        new = new.set_index(new["function"])
        del new["function"]
        new = new.T
        new = new.loc[:, ((new == 0).sum(axis=0) <= (sample*20) // 100)]
        new = new.T
        file.close()
        new.to_csv(fileName, sep="\t")
        return new


    def concatFiles(self):
#        dirList = os.listdir(self.pathToSamples)
#        names = []
#        os.system('mkdir concatFiles')
#        for file in dirList:
#            if file.split('_')[0] not in names:
#                names.append(file.split('_')[0])
#            else:
#                os.system('cat '+self.pathToSamples+'/'+ file.split('_')[0] + '*.fastq.gz > concatFiles/' + file.split('_')[0] + '.fastq.gz')
#                os.system('rm '+self.pathToSamples+'/'+ file.split('_')[0] + '*.fastq.gz')
        self.pathToSamples="/nfs/data/functional_k_mers/concatFiles"
#
#        dirList = os.listdir(self.pathToTest)
#        names = []
#        os.system('mkdir concatFiles_test')
#        for file in dirList:
#            if file.split('_')[0] not in names:
#                names.append(file.split('_')[0])
#            else:
#                os.system('cat ' + self.pathToTest + '/' + file.split('_')[0] + '*.fastq.gz > concatFiles_test/' +
#                          file.split('_')[0] + '.fastq.gz')
        self.pathToTest = "/nfs/data/functional_k_mers/concatFiles_test"
        print("----- paired-end fastq files successfully concatenated -----")


    def get_files_fromdir(self,inputDir):
        inputDisease = self.disease
        kmerLength = self.kmerLength
        listWithGzFiles = []
        dirList = os.listdir(inputDir)
        for file in dirList:
            name = os.path.basename(file)
            listWithGzFiles.append(inputDir + "/" + name)
        return listWithGzFiles


    def run_humann(self,pathToFiles, pathToAbundance):
        print("----- starting humann-----")
        start_time = datetime.now()
        dirList = os.listdir(pathToFiles)
        for file in dirList:
            pathToFile = pathToFiles+"/"+os.path.basename(file)
            os.system('humann --input '+pathToFile+' --output '+pathToAbundance+' --resume --threads 8 --bypass-translated-search')
        end_time = datetime.now()
        print('Duration humann: {}'.format(end_time - start_time))
        return pathToAbundance


    def prepare_test_set(self, listOfFiles,listOfKmers):
        sample = 0
        kmerDB = {}
        new = pd.DataFrame()
        countFiles = []
        self.ntest = len(listOfFiles)
        print(listOfFiles)
        for filepath in sorted(listOfFiles):
            kmerCountFile = 'kmers_in_sample_' + str(sample)
            os.system('kmc -cs1000000000 -k' + str(self.kmerLength) + '  ' + filepath + ' ' + kmerCountFile + ' . >>kmc.log')
            kmerTextFile = kmerCountFile + '.txt'
            os.system('kmc_tools transform ' + kmerCountFile + ' dump ' + kmerTextFile)
            os.system('rm ./*_suf')
            os.system('rm ./*_pre')
            kmPath = './' + kmerTextFile
            countFiles.append(kmPath)
            with open(kmPath, "r") as tmpFile:
                for line in tmpFile:
                    kmer = line.split("\t")[0].rstrip()
                    count = line.split("\t")[1].rstrip()
                    if kmer in listOfKmers:
                      #  print("found")
                        if kmer not in kmerDB.keys():
                            counts = [0] * self.ntest
                            counts[sample] = count
                            kmerDB[kmer] = counts
                        else:
                            kmerDB[kmer][sample] = count
            sample += 1
            tmpFile.close()
            os.system('rm ' + kmPath)
            print(str(sample) + " processed")
        for entry in listOfKmers:
            if entry not in kmerDB.keys():
                kmerDB[entry]=[0] * self.ntest

        for i in range(0, len(listOfFiles)):
            tmp = listOfFiles[i].split("/")
            listOfFiles[i] = tmp[2].split(".")[0]

        pf1 = pd.DataFrame.from_dict(kmerDB, columns=sorted(listOfFiles), orient='index')
        ind = pf1.index.tolist()
        pf1.insert(0, "", [i for i in range(0, len(kmerDB))], True)
        pf1.insert(0, "kmer", ind, True)
        pf1 = pf1.set_index(pf1[""])
        del pf1[""]
        pf1 = pf1.set_index(pf1["kmer"])
        del pf1["kmer"]
        pf1 = pf1.loc[listOfKmers]
        self.kmerTest = pf1
        pf1.to_csv('countmatrix_test_1206.tsv', sep="\t")


    def train_model(self):

        model_map={}

        kmers=self.kmerMatrix.index.tolist()
        pathways = self.pathwayMatrix.index.tolist()

        tmp = self.kmerMatrix #.apply(np.random.permutation, axis=1)
        #tmp = pd.DataFrame(tmp.tolist(), index=self.kmerMatrix.index, columns=self.kmerMatrix.columns)

        pathMatrixTrans = pd.DataFrame(preprocessing.normalize(self.pathwayMatrix),columns = self.pathwayMatrix.columns,index=self.pathwayMatrix.index)
        kmerMatrixTrans = pd.DataFrame(preprocessing.normalize(tmp),columns = tmp.columns,index=tmp.index)

        kmers = kmerMatrixTrans.index.tolist()

        rows = pathMatrixTrans.shape[0]
        for i in range(0, rows):

            pathMatrixTrans_one = pathMatrixTrans.iloc[[i]].T
            kmerMatrixTrans_one = kmerMatrixTrans.T

            kmerMatrixTrans_one.insert(0, "Sample", kmerMatrixTrans_one.index, True)
            pathMatrixTrans_one.insert(0, "Sample", pathMatrixTrans_one.index, True)

            merged = pathMatrixTrans_one.merge(kmerMatrixTrans_one, on='Sample', how="outer")
            merged = merged.set_index(merged["Sample"])

            del merged["Sample"]
            X, y = merged.iloc[:, 1:], merged.iloc[:, 0]

            alphas = [0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1]
            elastic_cv = ElasticNetCV(alphas=alphas, cv=10)
            model = elastic_cv.fit(X, y)
            self.model_map[pathways[i]]=model

        return kmers


    def prepare_train(self):
        self.concatFiles()
#        self.pathabundanceDir = self.run_humann(self.pathToSamples,"humann_out_train")
#        self.pathwayTrain = self.read_pathcoverage(self.pathabundanceDir, "pathway_train.tsv")
#        print(self.pathwayTrain)
        file_list = self.get_files_fromdir(self.pathToSamples)
        self.count_kmers(file_list)
#        self.kmerMatrix.columns=self.pathwayMatrix.columns



    def predict_test(self):
       # kmers = self.train_model()
#        print("done training")
 #       self.testbundanceDir = self.run_humann(self.pathToTest, "humann_out_test")
 #       self.testabundaceDir = "humann_out_test"
  #      self.pathwayTest = self.read_pathcoverage("/nfs/home/students/a.miroshina/bachelor/tmpfiles/humann_out_test", "pathway_ibd.tsv")
        file_list = self.get_files_fromdir(self.pathToTest)
    #    print(sorted(file_list))
    #    for i in sorted(file_list):
    #        print(i.split("/")[-1])
        self.prepare_test_set(file_list,self.kmerset)
        print('done test')




if __name__ == "__main__":
    start_time = datetime.now()
    kmer=Kmers("/nfs/home/students/a.miroshina/arbeit/train_data/","/nfs/home/students/a.miroshina/arbeit/test_data/",23,"demo","/nfs/home/students/a.miroshina/arbeit/kmer_matrix.tsv")
    #kmer.pathToSamples = "/nfs/home/students/a.miroshina/bachelor/concatFiles/"
    #kmer.pathabundanceDir = "/nfs/home/students/a.miroshina/arbeit/humann3_HMP2/"
#    kmer.kmerMatrix=pd.read_csv("/nfs/home/students/a.miroshina/bachelor/countmatrix_train_22.tsv", sep='\t')
#    kmer.pathwayMatrix=pd.read_csv("/nfs/home/students/a.miroshina/bachelor/pathway_train.tsv", sep='\t')
#    kmer.kmerTest=pd.read_csv("/nfs/home/students/a.miroshina/bachelor/countmatrix_test.tsv", sep='\t')
#    kmer.pathwayTest=pd.read_csv("/nfs/home/students/a.miroshina/bachelor/pathway_test.tsv", sep='\t')

#    kmer.pathToTest = "/nfs/home/students/a.miroshina/bachelor/concatFiles_healthy"
#    kmer.pathToSamples = "/nfs/home/students/a.miroshina/bachelor/concatFiles"

#    kmer.pathToTest = "/nfs/home/students/a.miroshina/bachelor/concatFiles_test"
#    kmer.kmerMatrix = kmer.kmerMatrix.set_index(kmer.kmerMatrix["kmer"])
#    kmer.kmerTest = kmer.kmerTest.set_index(kmer.kmerTest["kmer"])
#    kmer.pathwayMatrix = kmer.pathwayMatrix.set_index(kmer.pathwayMatrix["function"])
#    kmer.pathwayTest = kmer.pathwayTest.set_index(kmer.pathwayTest["function"])

#    del kmer.pathwayMatrix["function"]
#    del kmer.pathwayTest["function"]
#    del kmer.kmerMatrix["kmer"]
#    del kmer.kmerTest["kmer"]

 #   kmer.kmerMatrix.columns = kmer.pathwayMatrix.columns
 #   print(kmer.pathwayMatrix)
 #   print(kmer.kmerMatrix)
    with open(kmer. path_to_kmerset,"r") as file:
        for line in file:
            kmer.kmerset.add(line.rstrip())
    #print(len(kmer.kmerset))
    kmer.prepare_train()
    kmer.predict_test()

    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

