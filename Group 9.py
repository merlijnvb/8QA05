"""
Course code: 8QA05
group number: 9
members: Iris Almekinders, Merlijn van Benthem, Dimo Devetzis, Elizabeth Ninh, Niels van Noort, Sam Reijs en Femke Schapendonk
"""


# libraries
import Fase1 as F1
# import Fase2 as F2
import Fase3 as F3

# FASE 1
outf_name_fase1 = F1.main()
print("The data of Phase 1 has been written to",outf_name_fase1)

# FASE 2
# F2.main() # can be assigned to variables still
# print("The data of Phase 2 has been written to",outf_name2) # give some user feedback that phse 2 has been completed successfully one way or another

#FASE 3
#variables
min_frequency = 2
min_length_substring = 4
f_clus_res = 'kmca_results.txt'
f_desc = 'GenDescription2.txt'
f_exprs = outf_name_fase1[0]
f_fam = 'CloneIdFamily.txt'

F3.main(f_clus_res, f_desc, f_exprs, f_fam, min_frequency, min_length_substring, outf_name_fase1[1])
