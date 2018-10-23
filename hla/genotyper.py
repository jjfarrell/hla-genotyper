#!/usr/bin/env python
# Author: John J. Farrell,Ph.D <farrell@bu.edu>.
# program: hla-genotyper.py
# purpose: Predicts 4 digit HLA genotype from next gen sequencing data
# Copyright John Farrell 2014 All Rights Reserved
import os
import sys
import math
import pysam
import sys
import version
import urllib2
# import hla
import numpy
import csv
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from optparse import OptionParser


# utility functions
def is_url(url):
  try:
    ret = urllib2.urlopen(url)
    if ret.code == 200:
       return True
    else:
       return False
  except:
       return False
def good_bq(read,cutoff):
    for b in read:
        if ord(b)-33 < cutoff:
           return  False
    return True

def doses(hla_genes_to_call,final,hla_want):
    hla_name=[]
    hla_dose=[]
    for g in hla_genes_to_call:
        hla_call=final[g]
        for h in sorted(hla_want):
           [hla_gene, hla_allele_type]=h.split("*",1)     
           if hla_gene==g:
               hla_name.append(h)
    for g in hla_genes_to_call:
        hla_call=final[g]
        [h1,h2]=hla_call
        for h in sorted(hla_want):
            [hla_gene, hla_allele_type]=h.split("*",1)
            if hla_gene==g:
                hla_dose.append(dose(h,h1,h2))
    return hla_name,hla_dose
    
def dose(h,h1,h2):
    d=0
    if h==h1:
        d+=1
    if h==h2:
        d+=1
    return d
def read_length(filename):
#
# determine read length of bam files reads (Take max of first 200 records)
#
    mylength=0
    bamfile=pysam.Samfile(filename,'rb')
    n = 0
    for alignedread in bamfile.fetch(until_eof = True):
              mylength = max(alignedread.rlen,mylength)
              n = n+1
              if n > 200:
                  break
    bamfile.close()
    return mylength

def is_hg_ref(filename):
#
# determine read length of bam files reads (Take max of first 200 records)
#
    mylength=0
    bamfile=pysam.Samfile(filename,'rb')
    n = 0
    for alignedread in bamfile.fetch(until_eof = True):
              n+=1
              tid=alignedread.tid
              chr=bamfile.getrname(tid)
              if n > 20:
                  break
    hg=False
    if chr[0:3]=="chr":
        hg=True
    bamfile.close()
    return hg

def SM(filename):
# get sample id from bam file
    bamfile=pysam.Samfile(filename,'rb')
    header=bamfile.header.copy()
    sm_id="0"
    if "RG" in header.keys():
       sm_id= header["RG"][0]["SM"]
    # if no sample id use bamfile name
    if sm_id=="0":
        sm_id=os.path.basename(filename)
        sm_id=sm_id.replace(".bam","")
    bamfile.close()

    return sm_id 

def CN(filename):
# get Center id from bam file
# get Center id from bam file (1000 Genomes Center use different Exome Kits)
    bamfile=pysam.Samfile(filename,'rb')
    header=bamfile.header.copy()
#    print header["RG"]
    cn_id="NA"
    if "CN" in header.keys():
        cn_id= header["RG"][0].get("CN","NA")
    bamfile.close()
    return cn_id

                                                                                                                            
def containsAny(str, set):
    """Check whether 'str' contains ANY of the chars in 'set'"""
    return 1 in [c in str for c in set]

class AutoVivification(dict):
      """Implementation of perl's autovivification feature."""
      def __getitem__(self, item):
        try:
              return dict.__getitem__(self, item)
        except KeyError:
              value = self[item] = type(self)()
              return value


# import os
# import csv
# import operator


def hla_freq(ethnicity,gene):
 this_dir, this_filename = os.path.split(__file__)
 ETHNIC_PRIORS= os.path.join(this_dir,  "data/ethnic_priors.txt")
 prior={}
# with open('data/ethnic_priors.txt', 'rb') as f:
 with open(ETHNIC_PRIORS, 'rb') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE,fieldnames=["ethnicity","gene","allele","freq"])
        for row in reader:
	 if(ethnicity == row["ethnicity"]):
#                         print row["ethnicity"], row["gene"], row["allele"],row["freq"]
#                  if float(row["freq"]) > 0.001:
                         if gene == row["gene"] :
                          prior[row["allele"]]=float(row["freq"])
                         elif  gene=="all" :
                          prior[row["allele"]]=float(row["freq"])
#                          prior[get_hla_gene(row["allele"])][row["allele"]]=float(row["freq"])
#                          print  row["allele"],prior[row["allele"]]
 f.close()
 sorted_prior = sorted(prior.iteritems(), key=operator.itemgetter(0))
 return prior
def hla_loci(ethnicity,gene):
 this_dir, this_filename = os.path.split(__file__)
 ETHNIC_PRIORS= os.path.join(this_dir,  "data/ethnic_priors.txt")
 hla_loci=set([])
 with open(ETHNIC_PRIORS, 'rb') as f:
        reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE,fieldnames=["ethnicity","gene","allele","freq"])
        for row in reader:
	 if(ethnicity == row["ethnicity"]):
#                  if float(row["freq"]) > 0.001:
                           hla_loci=hla_loci|set([get_hla_gene(row["allele"])])
 f.close()
 return sorted(hla_loci)

def hla_typing_exons(gene):
    e={"HLA-DRB1": ["2"],
       "HLA-DQA1": ["1","2","3","4"],
       "HLA-DQB1": ["1","2","3","4"],
       "HLA-A": ["1","2","3","4"],
       "HLA-B": ["1","2","3","4"],
       "HLA-C": ["1","2","3","4"]
       }
    return e[gene]
def hla_alleles_in_gene(alleles,want_gene):
    hla_alleles=[]
    for a in alleles:
        allele_gene=get_hla_gene(a)
        if  allele_gene==want_gene:
            hla_alleles.append(a)
    return sorted(hla_alleles)
def hla_exon_region(gene,exon,EXON_INFO):
    start=0;
    stop=0;
    with open(EXON_INFO, 'rb') as f:
         reader = csv.DictReader(f, delimiter='\t', quoting=csv.QUOTE_NONE,fieldnames=["gene","exon","chr","start","stop"])
         for row in reader:
             if gene == row["gene"]:
                 if exon == row["exon"]:
#                     region= [row["start"],row["stop"]]
                     start= int(row["start"])
                     stop=int(row["stop"])
    f.close()
    return start,stop
def get_hla_gene(hla_allele):
    [gene,allele_number]=hla_allele.split("*")
    return gene
def n_hla_genes(alleles):
    gene_set=set([])
    for a in alleles:
        gene_set=gene_set|set([get_hla_gene(a)])
    return len(gene_set)
def main(argv):
      this_dir, this_filename = os.path.split(__file__)
      HLA_DATA = os.path.join(this_dir,  "data/hla.dat")
      if len(argv) == 1:
          argv.append("--help")
#      print "EP:",ETHNIC_PRIORS
      usage = "usage: %prog [options] BAMFILE"
      desc="hla-genotyper "+version.version+" predicts HLA genotypes from RNA-Seq and DNA-Seq bam files."
      parser = OptionParser(usage=usage,description=desc)
#      parser.add_option("-b", "--bam", action="store", type="string", dest="bamfile",help="BAM file (sorted and indexed by position)")
      parser.add_option("-u", "--unmapped_bam", action="store", type="string", default="None",dest="unmapped_bamfile",help="BAM file of unmapped READs (optional)")
      parser.add_option("-e", "--ethnicity", action="store", dest="ethnicity",type='choice',choices=['EUR','AFA','HIS','API','UNK'],help="Ethnicity of sample: EUR,AFA,HIS,API,UNK (required)")
      parser.add_option("-o", "--outdir", action="store", dest="outdir",type="string",help="Output directory")
      parser.add_option("-r", "--reference", action="store", dest="ref",type='choice',choices=['37','38'],default='37',help="Reference used for BAM alignment: 37 or 38 [default: %default])")
      parser.add_option("-m", "--mapq", action="store", type="int", dest="mapq",default=20,help="Minimum mapping quality of aligned reads  [default: %default]")
      parser.add_option("-q", "--qual", action="store", type="int", dest="bq",default=5,help="Min base quality for read  [default: %default]")
      # parser.add_option("-v", "--verbose", action="store_true", dest="verbose")
      parser.add_option("--debug", action="store_true", dest="debug")
      parser.add_option("--rnaseq", action="store_true", dest="rnaseq",help="Indicates RNA-Seq BAM file")
      parser.add_option("--exome", action="store_true", dest="exome",help="Indicates Whole Exome BAM file")
      parser.add_option("--genome", action="store_true", dest="genome",help="Indicates Whole Genome BAM file")
      parser.add_option("-s", "--sample", action="store", type="string", dest="sample",help="To specifiy Sample ID for reports")
      parser.add_option("--cn", action="store", type="string", dest="cn",default="None",help="Sequencing Center")
      parser.add_option("-l", "--len", action="store", type="int", dest="readlen", default=0,help="READ Length (optional)  [default: %default]")
      parser.add_option("-g", "--gene", action="store", dest="gene",default="all")

      (options, args) = parser.parse_args()
# check options 
      if len(args)==1:
         options.bamfile=args[0]
      else:
         parser.print_help()
         print "ERROR: BAMFILE argument is required"
         print "usage: hla-genotyper -e [EUR,AFA,API,HIS,UNK] [--genome,--rnaseq,--exome]  [options] myfile.bam"
         print "For help, hla-genotyper --help"     
         exit(-1)
      print "HLA-Genotyper Version:",version.version
      if options.bamfile == None:
#         parser.print_help()
          print "ERROR: BAMFILE is a required argument"
          print "usage: hla-genotyper -e [EUR,AFA,API,HIS,UNK] [--genome,--rnaseq,--exome]  [options] myfile.bam"
          print "For help, type: hla-genotyper --help"     
          exit(-1)
      else:
        # check for file
        if os.path.isfile(options.bamfile) is False and is_url(options.bamfile) is False:
          print
          print "ERROR: The bam file specified ", options.bamfile, "was not found."
          print "usage: hla-genotyper -e [EUR,AFA,API,HIS,UNK] [--genome,--rnaseq,--exome]  [options] myfile.bam"
          print "For help, type: hla-genotyper --help"     
          print
          exit(-1)
      if options.ethnicity==None:
          print
          print "ERROR: ethnicity is a required parameter"
          print "usage: hla-genotyper -e [EUR,AFA,API,HIS,UNK] [--genome,--rnaseq,--exome]  [options] myfile.bam"
          exit(-1)
      if options.unmapped_bamfile!="None":
          if os.path.isfile(options.unmapped_bamfile) is False and is_url(options.unmapped_bamfile) is False:
             print
             print "ERROR: The unmapped bam file specified ", options.unmapped_bamfile, "was not found."
             print
             parser.print_help()
             exit(-1)
      if options.sample == None:
          study_id=SM(options.bamfile)
      else:
          study_id=options.sample
      print "Sample=",study_id            
      if options.cn == "None":
          cn_id=CN(options.bamfile)
      else:
          cn_id=options.cn
      print "Sequencing Center=",cn_id
      # open log and output files
      prefix=""
      if options.rnaseq:
          prefix="hla."+options.ethnicity+"."+study_id+".rnaseq"
      elif options.exome:
          prefix="hla."+options.ethnicity+"."+study_id+".exome"
      elif options.genome:
          prefix="hla."+options.ethnicity+"."+study_id+".genome"
      else:
          prefix="hla."+options.ethnicity+"."+study_id
 
      if prefix=="":
          print 
          print "ERROR: Required option missing. Please specify sequencing type with --rnaseq, --exome, or --genome"
          print
          exit(-1)
      else:
          if options.outdir != None: 
              if os.path.isdir(options.outdir):
                  prefix=options.outdir+"/"+prefix
              else:
                  print "ERROR: Output directory,",options.outdir,", does not exist."
      fout=open(prefix+".txt", "w")
      out = csv.writer(fout, delimiter="\t")
      out.writerow(["#CN","BAM file","sample","ethnicity","gene","a1","a2","pp","qual","a1 mr","a1 ur","a2 mr","a2 ur","a1+a2"])
      flog=open(prefix+".log", "w")
      fdose=open(prefix+".dose", "w")
      print >>flog,"HLA-Genotyper Version:",version.version
      print >>flog,"Sample=",study_id
      print >>flog,"BAM input file to Scan (mapped reads)  :",options.bamfile
      print >>flog,"BAM input file to Scan (unmapped reads):",options.unmapped_bamfile
      print >>flog,"Ethnicity:",options.ethnicity
      print >>flog,"Base Quality Cutoff:",options.bq
      
      chr6="6"
      if is_hg_ref(options.bamfile):
        chr6="chr6"  
      EXON_INFO= os.path.join(this_dir,  "".join(["data/exon_info_",options.ref,".txt"]))
      
      hla_prior=hla_freq(ethnicity=options.ethnicity,gene=options.gene)
      if options.gene=="all":
        hla_genes_to_call=hla_loci(options.ethnicity,gene=options.gene)
      else:
        hla_genes_to_call=["HLA-"+options.gene]
      print >>flog,"HLA Loci to call: ", hla_genes_to_call
      want_exons={}
      for g in hla_genes_to_call:
          want_exons[g]=hla_typing_exons(g)
          print >>flog,"Exons read for genotyping ",g,": ",want_exons[g]
      # want_exons=hla.typing_exons("HLA-"+options.gene)
      
      
      
      
      
      # print "SM:",SM(options.bamfile)
      # 
      basename=os.path.basename(options.bamfile)
      print "Opening",options.bamfile
      
      if options.readlen==0:
          options.readlen=read_length(options.bamfile)
      print >>flog,"Read lengths:",options.readlen
      read_len=options.readlen
      # print >>flog,"Verbose setting:",options.verbose
      this_dir, this_filename = os.path.split(__file__)
      HLA_DATA = os.path.join(this_dir,  "data/hla.dat")
      print >>flog,"IMGT Database=",HLA_DATA
      
      # init variables
      hla_want=hla_prior.keys()
      hla_total={}
      n_hla_alleles=0;
      n_hla_reads=0;
      hla_read = AutoVivification()
      print "Reading HLA.dat database"
      # read in hla database and create hash of all possible reads whose length match bam file
      for seq_record in SeqIO.parse(HLA_DATA, "imgt"):
      # parse hla info from record description
      #  print  seq_record.description
        if "HLA" in str(seq_record.description) :
           [hla_allele, mhc_region] = seq_record.description.split(", ",1)
           [hla_gene, hla_allele_type]=hla_allele.split("*",1)
      #    not all alleles are 4 digit NN:NN some are NN:NNN
           hv=hla_allele_type.split(":")
           hla_4digit=hla_gene+"*"+hv[0]+":"+hv[1]
      #     print hla_4digit
           if (hla_4digit in hla_want):
      #     if hla_gene in ["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DQA1","HLA-DQB1"]:
              n_hla_alleles +=1
              hla_total[hla_4digit]=hla_total.get(hla_4digit,0)+1
      #
      # For rna-seq, the exons are spliced together.
      #
              if options.rnaseq:
                exons = [f for f in seq_record.features if f.type == "exon"]
                exon_seq=""
      # splice exons together and then split into potential reads
                last_position=0
                for e in exons :
                  if e.qualifiers["number"][0] in want_exons[hla_gene]:
      # splice exons together depending on whether on positive or negative strand                
                     if  e.location.start.position >last_position:
                          exon_seq=exon_seq+str(seq_record.seq[e.location.start.position : e.location.end.position])
                     else:
                          exon_seq=str(seq_record.seq[e.location.start.position : e.location.end.position])+exon_seq
                     last_position=e.location.start.position
                for i in range(0,len(exon_seq)-read_len+1):
                       myread=exon_seq[i:i+read_len]
      #                 if hla_4digit == "HLA-B*50:01" :
      #                     print "debug e=",e.qualifiers["number"][0],i,i+read_len,str(seq_record.description),hla_allele,hla_4digit,myread
      #                 hla_read[myread][hla_4digit]=hla_read[myread].get(hla_4digit,0)+1
                       hla_read[myread][hla_4digit]=1
                       n_hla_reads += 1
              else:
      # loop through exon  sequences used to predict HLA reads
                exons = [f for f in seq_record.features if f.type == "exon"]
                introns = [f for f in seq_record.features if f.type == "intron"]
                extra=0
                if  hla_gene in ["HLA-C","HLA-DQA1","HLA-DQB1","HLA-DRB1"] and len(introns) >2 :
      #              extra=50  # caused problems when one hla had introns and another did not
                    extra=0  #test hla-drb1*08:01 calls
      #              print "intron len",hla_gene,len(introns),len(exons)
                for e in exons :
                  if e.qualifiers["number"][0] in want_exons[hla_gene]:
                     exon_seq=seq_record.seq[e.location.start.position-extra: e.location.end.position+extra]
                     for i in range(0,len(exon_seq)-read_len+1):
                       myread=str(exon_seq[i:i+read_len])
                       hla_read[myread][hla_4digit]=1
                       n_hla_reads += 1
      
      #                 for k,v in hla_read.iteritems():
      #                   print k,v
      
                                                     
      print >>flog,"Number of HLA Alleles=",n_hla_alleles
      print >>flog,"Number of 4 Digit HLA Alleles=",len(hla_total)
      print >>flog,"Number of HLA Reads=",n_hla_reads
      print >>flog,"N unique Reads=",len(hla_read)
      print >>flog,"-"*80
      print >>flog,"Scanning Mapped Reads"
      # Read sequences from bam file
      import pysam
      bamfile=pysam.Samfile(options.bamfile,'rb')
      # init variables
      bam_hla_reads=[]
      mapped_read_total= AutoVivification()
      unmapped_read_total= AutoVivification()
      bam_hla_allele_set=set([])
      
      # start scanning Mapped Reads looking for exact reads found in  exons of specified HLA Genes
      print "Scanning mapped reads"
      gene_read_total={}
      for g in hla_genes_to_call:
       gene_read_total[g]=0
       for exon_number in want_exons[g]:
        [start,stop]=hla_exon_region(g,exon_number,EXON_INFO)
        print >>flog,"Scanning mapped reads for exon ",exon_number," gene:",g," start:",start," stop:",stop,
        exon_reads=0
        for alignedread in bamfile.fetch(chr6,start,stop) :
         my_seq = Seq(alignedread.seq, IUPAC.unambiguous_dna)
      # for archaic genomes-trim to some standard length
         my_seq =my_seq[0:read_len]
      #  filtering out reads with very poor quality 
      #  many filtered since they are not an exact match
         if good_bq(str(alignedread.qual),options.bq) and  alignedread.mapq >options.mapq and not alignedread.is_duplicate and not alignedread.is_unmapped :
           if str(my_seq) in hla_read:
              if n_hla_genes(hla_read[str(my_seq)])==1:
                  if options.debug: # and len(hla_read(str(my_seq))) ==1:
                        print "my seq ", str(my_seq), "gene info",n_hla_genes(hla_read[str(my_seq)]),"qual=",alignedread.qual,"Alleles=",hla_read[str(my_seq)]
                  exon_reads+=1 
                  bam_hla_reads.append(str(my_seq))
                  read_hla_alleles=hla_read.get(str(my_seq))            
                  bam_hla_allele_set= bam_hla_allele_set|set(read_hla_alleles)
                  
                  for a in read_hla_alleles:
                    mapped_read_total[a]=mapped_read_total.get(a,0)+1
           reverse_comp=str(my_seq.reverse_complement())
      #      print "my seqr", reverse_comp, "gene info",hla.n_hla_genes(hla_read[str(my_seq)]),hla_read[str(reverse_comp)]      
           if reverse_comp in hla_read:
              if n_hla_genes(hla_read[reverse_comp])==1:
                  if options.debug: # and  len(hla_read(reverse_comp)) ==1:
                        print "my seqr", reverse_comp, "gene info",n_hla_genes(hla_read[reverse_comp]),"qual=",alignedread.qual,"Alleles=",hla_read[reverse_comp]
                  exon_reads+=1 
                  bam_hla_reads.append(reverse_comp)
                  read_hla_alleles=hla_read.get(reverse_comp)
                  bam_hla_allele_set= bam_hla_allele_set|set(read_hla_alleles)
                  for a in read_hla_alleles:
                    mapped_read_total[a]=mapped_read_total.get(a,0)+1
        print >>flog,"N mapped reads=",exon_reads
        gene_read_total[g]= gene_read_total[g]+exon_reads
      bamfile.close()
      bam_hla_reads_set=set(bam_hla_reads)
      mapped_reads=len(bam_hla_reads_set)
      
      if options.unmapped_bamfile=="None":
        print "Skipping unmapped read scan"
      else:
      
      # Now scan through all unmapped reads
      # Finds exact match to known HLA gene exon region
       unmapped_bamfile=pysam.Samfile(options.unmapped_bamfile,'rb')
       print "Scanning unmapped reads"
       for alignedread in unmapped_bamfile.fetch(until_eof = True) :
      #  if not "#" in str(alignedread.qual):
        my_seq = Seq(alignedread.seq, IUPAC.unambiguous_dna)
        my_seq =my_seq[0:read_len]
        if good_bq(str(alignedread.qual),options.bq) and not alignedread.is_duplicate and alignedread.is_unmapped:
            
      #      print "my seq ", str(my_seq)
            reverse_comp=str(my_seq.reverse_complement())
            if str(my_seq) in hla_read:
               if n_hla_genes(hla_read[str(my_seq)])==1:
                  if options.debug:
                        print "my seq ", str(my_seq), "gene info",n_hla_genes(hla_read[str(my_seq)]),"qual=",alignedread.qual,hla_read[str(my_seq)]
                  bam_hla_reads.append(str(my_seq))
                  read_hla_alleles=hla_read.get(str(my_seq))
                  bam_hla_allele_set= bam_hla_allele_set|set(read_hla_alleles)
                  for a in read_hla_alleles:
                    unmapped_read_total[a]=unmapped_read_total.get(a,0)+1
                    g=get_hla_gene(a)
                  gene_read_total[g]= gene_read_total[g]+1
            if reverse_comp in hla_read:
              if n_hla_genes(hla_read[reverse_comp])==1:
                  if options.debug:
                        print "my seqr", reverse_comp, "gene info",n_hla_genes(hla_read[reverse_comp]),"qual=",alignedread.qual,hla_read[reverse_comp]
                  bam_hla_reads.append(reverse_comp)
                  read_hla_alleles=hla_read.get(reverse_comp)
                  bam_hla_allele_set= bam_hla_allele_set|set(read_hla_alleles)
                  for a in read_hla_alleles:
                    unmapped_read_total[a]=unmapped_read_total.get(a,0)+1
                    g=get_hla_gene(a)
                  gene_read_total[g]= gene_read_total[g]+1
       unmapped_bamfile.close()
      
      bam_hla_reads_set=set(bam_hla_reads)
      unmapped_reads=len(bam_hla_reads_set)-mapped_reads
      print >>flog,"-"*80
      print >>flog,"Total Read Summary"
      print >>flog,"HLA Genes to Call:",hla_genes_to_call
      print >>flog,"Total Exact match reads from HLA Genes:",len(bam_hla_reads)
      print >>flog,"Unique Exact match reads from HLA Genes",len(bam_hla_reads_set)
      print >>flog,"Reads mapped to "+ str(len(bam_hla_allele_set))+" HLA Alleles"
      # print >>flog,"bam_hla_allele_set:",bam_hla_allele_set
      
      hla_genotypes=[]
      hla_prob = AutoVivification()
      #   print hla_read[r]
      final={}
      print >>flog,"-"*80
      print >>flog,"HLA Genotyping Results"
      print >>flog,"#CN","bam file","sample","ethnicity","gene","a1","a2","pp","qual","a1 mr","a1 ur","a2 mr","a2 ur","a1+a2"
      min_prob=0.01
      total_g_reads={}
      total_h_reads={}
      allele_list={}
      for g in hla_genes_to_call:
           allele_list[g]=[]
           total_g_reads[g]=0
           bam_gene_alleles=hla_alleles_in_gene(bam_hla_allele_set,g)
           for h in sorted(bam_gene_alleles):
               total_g_reads[g]=total_g_reads[g]+mapped_read_total.get(h,0)+unmapped_read_total.get(h,0)
               total_h_reads[h]=mapped_read_total.get(h,0)+unmapped_read_total.get(h,0)
           max_freq=0.0
           for h in sorted(bam_gene_alleles):
                  freq=float(total_h_reads[h])/float(total_g_reads[g])
                  if  freq >max_freq:
                      max_freq=freq
           for h in sorted(bam_gene_alleles):
                  freq=float(total_h_reads[h])/float(total_g_reads[g])
                  # Filter out low freq not close to 0.5 (not likely to be one of 2 heterozygous alleles)
                  # Reduces search space
                  if  freq > (max_freq/5.0):      
                      allele_list[g]=allele_list[g]+[h]
                  # else:
                  # For debugging
                  #   print "Dropping allele",g,h,total_h_reads[h],total_g_reads[g],"percent",float(total_h_reads[h])/float(total_g_reads[g])
           
      #     print allele_list[g]     
      for g in hla_genes_to_call:
       n=0
      # bam_gene_alleles=hla.alleles_in_gene(bam_hla_allele_set,g)
       bam_gene_alleles=allele_list[g]
       print "Calculating Probabilities for",g,"Genotype Calls"
       for h1 in sorted(bam_gene_alleles):
         for h2 in sorted(bam_gene_alleles):
            if (h1 <= h2):
              n=n+1
      #        prior_val=math.log(1-((1-hla_prior[h1])*(1-hla_prior[h2])))
              prior_val=math.log(hla_prior[h1]*hla_prior[h2])
      #        hla_prob[h1][h2]=prior_val
              for r  in  bam_hla_reads_set:
                   num_hits=float(len(hla_alleles_in_gene(hla_read[r],g)))
                   if num_hits > 0:
                      p=(hla_read[r].get(h1,0))/(num_hits)
                      q=(hla_read[r].get(h2,0))/(num_hits)
      # other genotype prob that were looked at
      #                genotype_prob=p*q
      #                genotype_prob=max(p,q)
                      genotype_prob=1-((1-p)*(1-q))
                      genotype_prob=max([genotype_prob,min_prob])
                      hla_prob[h1][h2]=hla_prob[h1].get(h2,prior_val)+numpy.log(genotype_prob)   
                   else:
                       genotype_prob=min_prob             
      #                hla_prob[h1][h2]=hla_prob[h1].get(h2,prior_val)+numpy.log(genotype_prob) 
      # print >>flog,g,"Genotypes probablities examined:",n
      
      # find maximum genotype likelihood
       max_prob=-99999
       hla_call=["NA","NA"]
      
       total_prob=-99999;
       for h1 in sorted(bam_gene_alleles):
         for h2 in sorted(bam_gene_alleles):
                  if ((h1 <= h2) and (hla_prob[h1][h2] >max_prob)):
                        max_prob=hla_prob[h1][h2]
                        hla_call=[h1,h2]
                  if (h1 <= h2):
                        total_prob=numpy.logaddexp(total_prob,hla_prob[h1][h2])
      #                  print h1,h2,hla_prob[h1][h2],total_prob
      # Printing Results to log file and hla results file
       h1= hla_call[0]
       h2= hla_call[1]
       total_reads=mapped_read_total.get(h1,0)+unmapped_read_total.get(h1,0)+mapped_read_total.get(h2,0)+unmapped_read_total.get(h2,0)
       delta=max_prob-total_prob
       pval=numpy.exp(delta)
      
       qc="Pass"
       if  pval < 0.66:  #Ambiguous set for below 0.66
          qc="Ambiguous"
       if gene_read_total[g] <10:  
          qc="Low Coverage"
          pval="NA"
       
       print  >>flog,basename,study_id,options.ethnicity,g,h1,h2,str(pval)[0:4], qc,mapped_read_total.get(h1,0),unmapped_read_total.get(h1,0),mapped_read_total.get(h2,0),unmapped_read_total.get(h2,0),total_reads
      # write out hla results
       line=[cn_id,basename,study_id,options.ethnicity,g,h1,h2, str(pval)[0:4], qc , mapped_read_total.get(h1,0),unmapped_read_total.get(h1,0),mapped_read_total.get(h2,0),unmapped_read_total.get(h2,0),total_reads]
      # line=[basename,study_id,options.ethnicity,g,h1,h2, str(pval)[0:4], qc , mapped_read_total.get(h1,0),unmapped_read_total.get(h1,0),mapped_read_total.get(h2,0),unmapped_read_total.get(h2,0),total_reads,gene_read_total[g]]

       out.writerow(line)
      # store final call  
       final[g]=hla_call
      
      
      # print detailed results to log
      
       if pval < 0.80:
         for h1 in sorted(bam_gene_alleles):
            for h2 in sorted(bam_gene_alleles):
                 if h1 <= h2:
                        pval_other=numpy.exp(hla_prob[h1][h2]-total_prob)
                        if(pval_other > 0.05 and pval_other <pval):
                            print  >>flog, basename,study_id,options.ethnicity,g, h1,h2  , str(pval_other)[0:4], qc,mapped_read_total[h1],unmapped_read_total[h1],mapped_read_total[h2],unmapped_read_total[h2],total_reads
      # Print Ped file
       
      
      flog.close()
      fout.close()
      # print dose file
      header=["bam","id","ethnicity"]
      line=[basename,study_id,options.ethnicity]
      [names,values]=doses(hla_genes_to_call,final,hla_want)
      header=header+names
      line=line+values
      out = csv.writer(fdose, delimiter="\t")
      out.writerow(header)
      out.writerow(line)
      fdose.close()
      
      print "HLA Genotypes written to file:",prefix+".txt"
      print "HLA Run Log written to file:  ",prefix+".log"
      print "HLA Doses written to file:    ",prefix+".dose"
      exit(0)      
if __name__ == "__main__":
     main(sys.argv)      
