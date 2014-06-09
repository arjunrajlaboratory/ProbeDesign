import bowtie_search
import fasta
import seq
import re
import math
import probe_design

lgr = fasta.Fasta('Lgr5.txt')
shark_tank = fasta.Fasta('shark_tank.fa')

for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'humanReference')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'humanReference')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)
 

for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'humanPseudo')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'humanPseudo')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)
 

for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'humanMito')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'humanMito')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)


for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'humanMito_rRNA')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'humanMito_rRNA')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)


for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'hg19')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'hg19')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)

for i in range(5,20):
   hts = bowtie_search.align_for_hits(shark_tank,i,'human_rDNA')
   sm1 = sum(hts)*1.0/len(hts)
   hts = bowtie_search.align_for_hits(lgr,i,'human_rDNA')
   sm2 = sum(hts)*1.0/len(hts)
   print(i,sm1,sm2,sm1/sm2)


