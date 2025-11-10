#This is a script to add a note from a blast table for Reciprocal best hit (RBH) or, if no RBH is availible, Blast best hit (BBH) to the corresponding gene in a gff3 file.

#The following packages need to be installed beforehand

#conda install bioconda::gffutils
#conda install anaconda::pandas

# Version of the used packages
# Gffutils 0.13
# Pandas 2.3.0
# Python 3.13.3

#inputs: .gff3 output file from EVidenceModeler, RBH and BBH tables


import gffutils as gf
import pandas as pd
from copy import deepcopy

gf.create_db("h2.EVM.gff3", dbfn="h2.EVM.db", keep_order=True, force=True, merge_strategy="create_unique")
gff = gf.FeatureDB('h2.EVM.db', keep_order=True)


columns = [
    "query", "target", "pident", "qstart", "qend", "qlen",
    "tstart", "tend", "tlen", "qcov", "tcov", "evalue", "theader"
]

rbh = pd.read_csv("h2Uniref100RBH", sep="\t", header=None, names=columns)

rbh_dict = rbh.drop_duplicates('query').set_index('query').to_dict(orient='index')

bbh = pd.read_csv("h2Uniref100BBH", sep="\t", header=None, names=columns)

bbh_dict = bbh.drop_duplicates('query').set_index('query').to_dict(orient='index')


temp = []
n_annotated = 0

all_mrnas = list(gff.features_of_type("mRNA"))
print(f"\nTotal mRNAs in GFF: {len(all_mrnas)}")

all_genes = list(gff.features_of_type("gene"))
print(f"\nTotal genes in GFF: {len(all_genes)}")

print("\nAnnotating genes...")


for mrna in gff.features_of_type("mRNA"):
    for queryid in mrna.attributes.get('ID', []):
        parent_genes = list(gff.parents(mrna, featuretype='gene'))
        for gene in parent_genes:
            attr = dict(gene.attributes) 
            if queryid in rbh_dict:
                note = rbh_dict[queryid]["theader"].replace(";", ",").replace("=", "-")              
                attr["RBH"] = [note]

            elif queryid not in rbh_dict and queryid in bbh_dict:
                note = bbh_dict[queryid]["theader"].replace(";", ",").replace("=", "-")              
                attr["BBH"] = [note]
            
            elif queryid not in rbh_dict and queryid not in bbh_dict:        
                attr["RBH"] = ["No"]
                attr["BBH"] = ["No"]

            gene_copy = deepcopy(gene)                   
            gene_copy.attributes = attr
            temp.append(gene_copy)
            n_annotated += 1  
            break  
              


print(f"\nAnnotated {n_annotated} genes with Note tags.")    

if temp:
    gff.update(temp, merge_strategy="replace")
    print("Database updated successfully.")
else:
    print("No gene annotations found — database not updated.")


with open("annotated_h2_gene.gff3", "w") as out:
    for feature in gff.all_features(order_by=(None)):
        out.write(str(feature) + "\n")
