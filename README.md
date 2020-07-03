# pubmed_to_COA
Write a Collaborating Authors table suitable for USDA/NSF COA/COI files. Note this will not cleanup issues where authors do or don't have middle initials included so the list may have some redundancies. This should be easy to cleanup manually once you have the list.

The prefix "A:" is following the [NSF COA policy](https://nsf.gov/bfa/dias/policy/coa.jsp) - if you need to simplify this in USDA format should be easy enough to do manually.  Arguably it should do better to pull out institutional affiliations to separate department name out, but this seems a bit hard to do well right now so it is just the whole institutional address from the publication.

EXAMPLES
=====
__*It is best to replace yourname@email.com with your actual email address!*__

Run a  search and print out the record IDs for the Biologist [Tyrone B Hayes](http://ib.berkeley.edu/people/faculty/hayest), write to the file `coi.tsv`:
`./pubmed_authors_to_COA.py -v -q '(Hayes TB OR Hayes Tyrone)' -e yourname@email.com`

Quietly run all and generate report in file `coi.tsv` for the Neuroscientist [Erich Jarvis](https://jarvislab.net/):

```./pubmed_authors_to_COA.py -q '(Jarvis Erich or Jarvis ED) AND (("2017/01/01"[Date - Publication] : "3000"[Date - Publication]))' -e yourname@email.com```

Quietly run all and generate report in file `Blackwell_M_co_authors.tsv` for the Mycologist Meredith Blackwell:
```./pubmed_authors_to_COA.py --query 'Blackwell Meredith' -o Blackwell_M_co_authors.tsv -e yourname@email.com```

USAGE
=====
```
usage: pubmed_authors_to_COA.py [-h] [-v] [--debug] -e EMAIL -q QUERY [-m MAX]
                                [-o OUT]

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print out Verbose processing
  --debug               Debug only printing one record
  -e EMAIL, --email EMAIL
                        Email address for Entrez query
  -q QUERY, --query QUERY
                        Pubmed Query
  -m MAX, --max MAX     Maximum number of records to retrieve
  -o OUT, --out OUT     Output file name
  ```
  
  AUTHOR
  =====
  [Jason Stajich](http://lab.stajich.org) [@hyphaltip](https://github.com/hyphaltip)
  
  REQUIREMENTS
  ============
  * [BioPython](https://biopython.org/)
