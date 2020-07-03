#!/usr/bin/env python3

# provide a COI list based on publication authors using Pubmed query
# Author: Jason Stajich (jason.stajich[AT]ucr.edu) / @hyphaltip

import csv, re, sys, os, argparse
import xml.etree.ElementTree as ET
from Bio import Entrez

Entrez.email = 'PROVIDE_EMAIL_AT_gmail.com'
outsamples="coi.tsv"

parser = argparse.ArgumentParser("pubmed_authors_to_COA.py",add_help=True)
parser.add_argument("-v","--verbose", help="Print out Verbose processing",action='store_true')
parser.add_argument("--debug", help="Debug only printing one record",action='store_true')
parser.add_argument("-e","--email", help="Email address for Entrez query",required=True)
parser.add_argument("-q","--query", help="Pubmed Query",required=True)
parser.add_argument("-m","--max", type=int, help="Maximum number of records to retrieve",default=100)
parser.add_argument("-o","--out", type=argparse.FileType('w'), help="Output file name",default="coi.tsv")
args = parser.parse_args()

Entrez.email = args.email

handle = Entrez.esearch(db="pubmed",term=args.query,retmax=args.max)
record = Entrez.read(handle)

splitline=re.compile(r'^\S+\s+\-\s+')
authors = {}
for publication in record["IdList"]:
    if args.verbose:
        print("processing PMID %s"%(publication))
    recordhandle = Entrez.efetch(db="pubmed", id=publication, retmode="text", rettype="medline")
    author = ""
    inst   = ""
    inAD = False
    for line in recordhandle:
        line = line.strip()
        if line.startswith("FAU"):
            author = splitline.sub('',line)
        elif line.startswith("AD"):
            inst = splitline.sub('',line)
            inAD = True
            if author not in authors:
                authors[author] = inst
        elif inAD and not re.match(r'^\S+\s+\-\s+',line):
            if args.debug:
                print("[DEBUG] address is '%s'"%(line))
            authors[author] += " " + line
        else:
            if args.debug:
                print("[DEBUG] inAD=%s skipping line %s"%(inAD,line))
            inAD = False
    if args.debug:
        break

outcsv    = csv.writer(args.out,delimiter="\t")
outcsv.writerow(["COLLAB","AUTHOR","INSTITUTION"])
for author in sorted(authors.keys()):
    outcsv.writerow(["A:",author,authors[author]])
