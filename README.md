# portal_files
INPUT FILE
To run the perl script you need an input files with the names of the projects

===================================
TO RUN ON COMMAND LINE
===================================
perl make_file_portal.pl Projects 


User inputs the given parameters in the JSON file: 
signatures => Signature file for the whole project
wgd_delivery => The FACETS file. The default is "/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_comb_mafs/WES_facets_estimates_WGD.v2.txt"
maf_file => MAF file for the project
gene_level => Gene level calls for the whole project
purity_seg => Segment file for the whole project
msi_score=> Standard file is "/ifs/res/taylorlab/chavans/roslin_2.4_deliveries/all_comb_mafs/exome_msi.txt" unless you want to change it.

Tumor_id_mapping=> Tab delimited mapping file:
P-0000059-T01-WES       s_C_000154_T001_d
P-0000067-T01-WES       s_C_000138_T001_d
P-0000086-T01-WES       s_C_58HT7F_M001_d
P-0000092-T01-WES       s_C_000002_T001_d


