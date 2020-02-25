# portal_files
INPUT FILE

To run the perl script you need an input files with the names of the projects


TO RUN ON COMMAND LINE

perl make_file_portal.pl Projects 


User inputs the given parameters in the JSON file: 

signatures => Signature file for the whole project

wgd_delivery => The FACETS file. The default is 

maf_file => MAF file for the project

gene_level => Gene level calls for the whole project

purity_seg => Segment file for the whole project

msi_file=>MSI file


Tumor_id_mapping=> Tab delimited mapping file:

XXXX-X01-WES       XXXX_d

XXXX-X01-WES       XXXX_d



Header/no header for the Tumor_id_mapping file  doesn’t matter.  
