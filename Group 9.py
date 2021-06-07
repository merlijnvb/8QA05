"""
Course code: 8QA05
group number: 9
members: Iris Almekinders, Merlijn van Benthem, Dimo Devetzis, Elizabeth Ninh, Niels van Noort, Sam Reijs en Femke Schapendonk
"""


# libraries
import Fase1 as F1
import Fase2 as F2
import Fase3 as F3

# FASE 1
outf_name1 = F1.main()
print("The data of Phase 1 has been written to",outf_name1)

# FASE 2
F2.main() # can be assigned to variables still
print("The data of Phase 2 has been written to",outf_name2) # give some user feedback that phse 2 has been completed successfully one way or another

#FASE 3
#variables
verwijderen = ['protein','similar', 'acidic', '4-like', '8-like', '-like', 'ESTs',
                      'like', 'acid-rich', 'adhesion', 'affinity', 'activity', 'alpha-like', 'anion', 'assembly',
                      'association', 'basic', 'candidate', 'carbon','Coenzyme', 'cofactor', 'coiled-coil',
                      'coiled-coil-helix-coiled-coil-helix', 'cold', 'double', 'domain-containing', 'enabled',
                      'fast', 'four', 'glucose', 'half', 'hand', 'inner', 'insert', 'inorganic', 'isoenzyme',
                      'molecule', 'mouse', 'never', 'neighbor', 'nitrogen', 'omega', 'only', 'organic', 'outer', 'paired',
                      'partner', 'region', 'ring', 'slow', 'similarity', 'system', 'very', '3-like', 'beta-like',
                      'coenzyme', 'complex', 'constant', 'component', 'dependent', 'early', 'light', 'long',
                      'protein-like', 'short', '1-like', 'activated', 'group', 'high', 'nine', 'small', 'cell',
                      'chain', 'heavy', 'with', 'acid', 'alpha', 'beta', 'associated', 'containing', 'gamma',
                      'gene', 'inter', 'rich', 'type', 'repeat']
min_frequency = 2
in_nr_clusters = 1
min_length_substring = 4
f_clus_res = 'kmca_results.txt'
f_desc = 'GenDescription2.txt'
f_exprs = 'Voorbeeld_clusterdata.txt'
f_fam = 'CloneIdFamily.txt'

F3.main(f_clus_res, f_desc, f_exprs, f_fam, min_frequency, in_nr_clusters, min_length_substring, verwijderen,outf_name1[1])