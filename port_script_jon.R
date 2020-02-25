########################################################################################################################
library(googlesheets)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
"%ni%" = Negate("%in%")

#######################################################################################################################

### Load cohort data

args = commandArgs(trailingOnly=TRUE)

#Arguments 8 AND 9 ARE data_clinical and data_sample
output_folder = args[10]
#output_folder = paste(args[6],"/output/", sep="")
#msi<-fread('XXXX_msi.txt')
msi<-fread(args[7])
names(msi)[1] = "Tumor_Sample_Barcode"
names(msi)[4] = "MSIscore"

#msi = msi %>% mutate(CMO_Sample_ID_fixed = str_replace_all(Tumor_Sample_Barcode,"s_","")) %>% 
 # mutate(CMO_Sample_ID_fixed = str_replace_all(CMO_Sample_ID_fixed,"_","-")) %>%
  #mutate(DMP_noIM = str_replace(CMO_Sample_ID_fixed,"-WES",""))

  length(unique(msi$Tumor_Sample_Barcode))

  #uniq_cmo_ids = unique(msi$DMP_noIM)

  data_clinical_sample = fread(args[8],skip = 4, fill=TRUE) %>% 
    filter(.,substr(SAMPLE_ID,1,13) %in% substr(msi$Tumor_Sample_Barcode,1,13)) %>% 
    select(PATIENT_ID,
	            SAMPLE_ID,
		             CANCER_TYPE,
			              CANCER_TYPE_DETAILED,
				               SAMPLE_TYPE,
					                METASTATIC_SITE,
							         PRIMARY_SITE,
								          ONCOTREE_CODE) %>%
      mutate(Pool=paste("Proj_",args[1],sep=""))
      #mutate(Pool="Proj_06049_U")

      data_clinical_patient = fread(args[9],skip = 4, fill=TRUE) %>% 
        filter(.,substr(PATIENT_ID,1,9) %in% substr(msi$Tumor_Sample_Barcode,1,9)) %>% 
        select(PATIENT_ID,
	                SEX,
			         PARTC_CONSENTED_12_245) 
	  #%>% mutate(Sex=str_replace_all(SEX,c("F"="Female","M"="Male"))) 

	head(data_clinical_patient)
	samplesheet=inner_join(data_clinical_patient,data_clinical_sample, by=c(PATIENT_ID = 'PATIENT_ID')) %>% 
	mutate(SAMPLE_ID_WES = str_c(substr(SAMPLE_ID,1,13),'-WES'))
	head(samplesheet); dim(samplesheet) 


	print("MASTERSHEET")
	master=samplesheet 
	dim(master);head(master)

	########################################################################################################################
	#FACETS CNA data GENELEVEL CALLS
	  cn_calls = fread(args[2]) %>% 
	  mutate(FACETS_CNA = ifelse(is.na(FACETS_CNA), NA, FACETS_CNA)) %>%
	  #print (cn_calls) 
	  dcast(Hugo_Symbol~Tumor_Sample_Barcode, value.var='FACETS_CNA')
	  #head(cn_calls)  #length(names(cn_calls)) #names(cn_calls) = str_replace(names(cn_calls),'IM6','WES')
	  write.table(cn_calls, paste(output_folder,'/Proj_',args[1],'_data_CNA.txt', sep=""), sep="\t",quote=F,row.names=F)
         ########################################################################################################################
	      ########################################################################################################################
	      ### Mutations
	      #cut -f1-109 
	      maf_ = fread(args[3]) %>%
	      mutate(Tumor_Sample_Barcode = str_replace_all(Tumor_Sample_Barcode,"s_","")) %>% 
	      mutate(Tumor_Sample_Barcode = str_replace_all(Tumor_Sample_Barcode,"_","-"))
		  
	      exonic_mutations = c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Silent',
				                          'Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Targeted_Region','Translation_Start_Site')

		maf = filter(maf_, Variant_Classification %in% exonic_mutations)
		#dim(maf)
		#head(maf)
		unique(maf$Tumor_Sample_Barcode) 

		write.table(maf, paste(output_folder,'/Proj_', args[1],'_data_mutations_extended.txt', sep=''), sep="\t",quote=F,row.names=F)

###################################################################################
	    #FACETS WGD
	    wgd_facets<-fread(args[4]) %>% 
	    mutate(Sample_ID = str_replace_all(Sample_ID,"s_","")) %>% 
	    mutate(Sample_ID = str_replace_all(Sample_ID,"_","-"))
	    #head(wgd_facets)
	    unique(wgd_facets$Sample_ID) 
	    
	    
	    
########################################################################################################################
	     #FACETS Segmentation

	     seg = fread(args[5]) %>% 
	     mutate(V1 = str_replace_all(V1,"_","-"))
	     setnames(seg,c('ID','chrom','loc.start','loc.end','num.mark','seg.mean'))
	     #dim(seg) head(seg)
             unique(seg$ID)
             #seg = seg %>% mutate(ID = plyr::mapvalues(ID, master$CMO_Sample_ID_fixed, master$DMP_Sample_ID)
	     write.table(seg, paste(output_folder,'/Proj_', args[1],'_data_cna_hg19.seg', sep=''), sep="\t",quote=F,row.names=F)
	     
	     
	########################################################################################################################


		#Mut_Sig
	     	  sig  = fread(args[6]) %>%
		  mutate(signatures=str_replace_all(signatures,",","|"),mean=str_replace_all(mean,",","|"),conf=str_replace_all(conf,",","|"))

		  head(sig)
		  unique(sig$sample_name) #rest did not have enough mutations

		  #Combine all

		  print("MASTER_MSI")
		  head(msi)
		  master_msi<-left_join(master,msi,by=c(SAMPLE_ID_WES='Tumor_Sample_Barcode'))
		  head(master_msi); dim(master_msi) #
		  print("MASTER_MSI_WGD")
		  master_msi_wgd<-left_join(master_msi,wgd_facets,by=c(SAMPLE_ID_WES='Sample_ID'))
		  head(master_msi_wgd); dim(master_msi_wgd) #

		  print(str(sig$Tumor_Sample_Barcode))	
		  #print(str(master_msi_wgd$SAMPLE_ID_WES))

		  if(length(sig$Tumor_Sample_Barcode) == 0){
	
				master_msi_wgd_sig = master_msi_wgd %>% mutate(mean=NA, conf =NA, signatures=NA)

		  }else{
		  	master_msi_wgd_sig = left_join(master_msi_wgd,sig,by=c(SAMPLE_ID_WES='Tumor_Sample_Barcode'))
		  	head(master_msi_wgd_sig); dim(master_msi_wgd) 

		  }

		  clin_sample = mutate(master_msi_wgd_sig) %>%
		      select(PATIENT_ID, SAMPLE_ID = SAMPLE_ID_WES, CANCER_TYPE, CANCER_TYPE_DETAILED, SAMPLE_TYPE, METASTATIC_SITE, PRIMARY_SITE, ONCOTREE_CODE,
			                FACETS_PURITY = Purity,
					           FACETS_PLOIDY = Ploidy,
						              GENOME_DOUBLED = WGD,
							                 FRACTION_CNA = FGA1,
									            MSI_Score = MSIscore,
										               Mutational_Signatures = signatures,
											                  Mutationa_Signatures_Proportions = mean,
													             Mutational_Signatures_Pval = conf
					           
					    ) %>%
		      mutate(PROJECT_CODE = paste('Proj_',args[1],sep=""),
			                INSTITUTE = 'MSKCC',
					           CUD_CATEGORY = 'N/A',
						              SAMPLE_CLASS = 'Tumor',
							                 LST = 'N/A',
									            NTAI = 'N/A',
										               HRD_LOH = 'N/A',
											                  BRCA_TYPE = 'N/A',
													             BRCA_VARIANT = 'N/A',
														                SOMATIC_BRCA_WT_STATUS = 'N/A'
					    )
		      head(clin_sample); dim(clin_sample) #


		      clin_patient=mutate(master_msi_wgd_sig) %>%
		                 select(PATIENT_ID = PATIENT_ID,
					                  SEX = SEX,
							                    `12_245_PARTC_CONSENTED` = PARTC_CONSENTED_12_245)
		      head(clin_patient); dim(clin_patient)


		      #write.table(clin_sample, '/ifs/res/taylorlab/pmishra/PORTAL/Proj_06049_U_data_clinical_sample_n3.txt',sep='\t',row.names=F,quote=F)
		      data_clinical_sample_file   = paste(output_folder,"/Proj_",args[1],'_data_clinical_sample.txt', sep="")
		      write.table(unique(clin_sample, decreasing=FALSE), data_clinical_sample_file ,sep='\t',row.names =F,quote=F)
		      data_clinical_patient_file  = paste(output_folder,"/Proj_",args[1],'_data_clinical_patient.txt', sep="")
		      write.table(unique(clin_patient, decreasing=FALSE), data_clinical_patient_file ,sep='\t',row.names=F,quote=F)

		      ### Sample map

		      sample_map = master_msi_wgd_sig %>% 
		        mutate(PATIENT_ID = PATIENT_ID, DMP_ID=SAMPLE_ID, SAMPLE_ID = SAMPLE_ID_WES) %>%
		        select(PATIENT_ID, SAMPLE_ID, DMP_ID)
			head(sample_map); dim(sample_map)
			sample_map_file = paste(output_folder,'/Proj_',args[1],'_sample_map.txt',sep="")
			write.table(unique(sample_map), sample_map_file, sep="\t",quote=F,row.names=F)

			##################################### END
