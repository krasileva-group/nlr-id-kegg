## Mapping IDs to KEGG pathways

### This pipeline uses following programs
- [KEGG database](https://www.genome.jp/kegg/kegg1.html) 
- [BLAST software suite](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [HMMER software suite](http://hmmer.org)
- [Bedtools toolset](https://bedtools.readthedocs.io/en/latest/)
- [NLRparser tools (including meme)](https://github.com/steuernb/NLR-Parser)

### The basic steps in the pipeline are following

0. Prepare all data files, including NLR database.
1. Subtract genes with canonical NLR domains from a plant proteome(s) of interest.
2. Run blast of a plant proteome(s) against NLRs.
3. Intersect canonical NLR domains with blast results to filter them out.
4. Run the reciprocal blast.
5. Map reciprocal blast hits to KEGG pathways.
