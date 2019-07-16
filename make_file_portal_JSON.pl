use strict;
use File::Basename;
use JSON::Parse 'json_file_to_perl';
use List::MoreUtils qw(uniq);
use Data::Dumper;

my @projects   = read_file(\$ARGV[0]);


#read_config_file(\"Exome_info.json");

########## ALL HARD CODED ####################
#my $cs_delivery      = $json_exome->{"input_dir"}[0]->{"cs_delivery"};
#my $cs_delivery_g    = $json_exome->{"input_dir"}[0]->{"cs_delivery_g"};
#my $mutation_file    = $json_exome->{"input_dir"}[0]->{"maf_file"};
#my $gene_level_calls = $json_exome->{"input_dir"}[0]->{"gene_level"};
#my $purity_seg       = $json_exome->{"input_dir"}[0]->{"purity_seg"};
#my $wgd_delivery     = $json_exome->{"input_dir"}[0]->{"wgd_delivery"};
#my $signature_calls  = $json_exome->{"input_dir"}[0]->{"signatures"};
#my $redacted         = $json_exome->{"input_dir"}[0]->{"redacted"};
#my $msi_score        = $json_exome->{"input_dir"}[0]->{"msi_score"};
#my $data_sample      = $json_exome->{"input_dir"}[0]->{"data_sample"};
#my $data_patient     = $json_exome->{"input_dir"}[0]->{"data_patient"};



#NEW WAYS
#my @new_checks = read_file(\"latest_roslin_runs_v241_WES_20190306.txt");

#This is where I am 
my $pwd =`pwd`;chomp $pwd;


open(LOG2, ">ALL_PROJECT_LOG_FINAL");
#Loop through each project and collect the input files and then run the R script for the output file.
#Also write into the log
foreach my $input_line (@projects){
	my ($proj, $JSON_file) = split("\t", $input_line);
	print "$JSON_file\n";
	my $json_exome = read_config_file(\$JSON_file);

	########## from JSON file  ####################
	my $cs_delivery      = $json_exome->{"input_dir"}[0]->{"cs_delivery"};
	my $cs_delivery_g    = $json_exome->{"input_dir"}[0]->{"cs_delivery_g"};
	my $mutation_file    = $json_exome->{"input_dir"}[0]->{"maf_file"};
	my $gene_level_calls = $json_exome->{"input_dir"}[0]->{"gene_level"};
	my $purity_seg       = $json_exome->{"input_dir"}[0]->{"purity_seg"};
	my $wgd_delivery     = $json_exome->{"input_dir"}[0]->{"wgd_delivery"};
	my $signature_calls  = $json_exome->{"input_dir"}[0]->{"signatures"};
	my $redacted         = $json_exome->{"input_dir"}[0]->{"redacted"};
	my $msi_score        = $json_exome->{"input_dir"}[0]->{"msi_score"};
	my $data_sample      = $json_exome->{"input_dir"}[0]->{"data_sample"};
	my $data_patient     = $json_exome->{"input_dir"}[0]->{"data_patient"};
	my $tumor_mapping    = $json_exome->{"input_dir"}[0]->{"Tumor_id_mapping"};
	my $roslin_mapping    = $json_exome->{"input_dir"}[0]->{"Roslin_mapping"};


	my @new_checks = read_file(\$roslin_mapping);
	chomp $proj;
	my $new_p = "Proj_".$proj;
	my %hash_log = ();
	my @type;
	#my $tumor_mapping = "Proj_07250_X_mapping";
	#my $tumor_mapping = "Tumor_id_mapping";
	#$hash_log{'Project'} = $proj;
	
	my @redacted = redaction_list(\$new_p, \$redacted );
	#Make directories:
	my $log     = "log";
	my $input   = "input";
	my $output  = "output";
	
	my $log_dir    = $pwd."/Proj_".$proj."/".$log;
	my $input_dir  = $pwd."/Proj_".$proj."/".$input;
	my $output_dir = $pwd."/Proj_".$proj."/".$output;


	make_directories(\$log_dir, \$input_dir, \$output_dir);
	my $log_file = mkdir_log(\$log_dir, \$proj);

	
	#Define files that need to be copied 
	print "================================= Proj_$proj ===============================\n";
	#my $mutation_file    = define_input_files(\$cs_delivery, \$proj, \"final_comb_ccs_022119.maf");
	#my $gene_level_calls = define_input_files(\$cs_delivery_g,\$proj,\"ccs_filtered_genelevel_cna_043019.txt");
	#my $purity_seg       = define_input_files(\$cs_delivery,\$proj,\"purity_seg_061119.txt");
	#my $signature_calls  = define_input_files(\$signatures,\$proj,\"fil_top3list_022119.txt");
	#my $msi_calls        = define_input_files(\$wgd_delivery,\$proj,\"exome_msi.txt");

	#Only for WGD give the  direct file
	#my $wgd_file         = $wgd_delivery."/WES_facets_estimates_WGD.v2.txt";
	
	my $wgd_file  = $wgd_delivery;
	#just say it exists:
	if (-e $wgd_file){
		print "$wgd_file exists\n";

	}else{
		print "$wgd_file does not exists\n";
	}


	if (( -e $mutation_file) && (-e $gene_level_calls) && (-e $purity_seg)  && (-e $signature_calls) && (-e $wgd_file)){
		my $new_proj = "Proj_".$proj."/input";
		#`mkdir $new_proj`;
		print "All the input files exist in project Proj_$proj for Rscript port_script_jon.R\n";
	
		#FIND ROSLIN CMO IDS 
		my @roslin_orig_cmo_ids   =  find_ids_from_original(\$proj, \@new_checks);
		#print "ROSLIN ORIGINAL IDS :@roslin_orig_cmo_ids\n";
		#MUTATION CALL FILES 
		`cp $mutation_file $new_proj`;
		 print "Copied mutation file:$mutation_file to $input_dir\n";
		 write_log (\$log_file, \"=================== Mutation file processing =====================\n");
		 write_log (\$log_file, \"Copied  MAF/mutation file from $mutation_file to  $pwd/$new_proj");


		 #remove_redactioons(\@redacted);

		 #Change the IDS
		 my $new_mut_file                    = find_input_files(\$new_proj,\"final_comb");
		
		 #my $new_mut_file                    = $new_proj."/Proj_".$proj."_final_comb_ccs_022119.maf";
		 my ($new_mut_mapped, $changed_muts, $maf_count) = change_mutation(\$tumor_mapping, \$new_mut_file, \$output_dir, \$log_file);
		 write_log (\$log_file, \"Replaced in $$new_mut_mapped:CMO IDs with DMP IDs");
		 #print "COMPARING CMO IDS\n";
		 my %hash1 = compare_CMO_IDs($changed_muts,\@roslin_orig_cmo_ids,\"MAF"); 
		 push (@type, "MAF");
		 $hash_log{"MAF"} = {%hash1};

		 #SUBSAMPLE data_clinical and data_patient
		 my $dmp_ids_in_proj  = extract_dmpid($changed_muts);
		 #print "DMP  ids in the project:  @$dmp_ids_in_proj\n";
		 my $new_patient_file = subset_data_clinical(\$data_patient,$dmp_ids_in_proj, \"patient",\$output_dir );
		 my $new_sample_file  = subset_data_clinical(\$data_sample,$dmp_ids_in_proj, \"sample",\$output_dir );
		
		 write_log(\$log_file, \"Made new data_clinical_sample.txt file: $new_sample_file");
		 write_log(\$log_file, \"Made new data_clinical_patient.txt file: $new_patient_file");


		 #GENE LEVEL CALLS
		 `cp $gene_level_calls $new_proj`;
		 print "Copied genelevel facets calls: $gene_level_calls to $pwd/$new_proj\n";
		 write_log (\$log_file, \"=================== FACETS Gene Level call  processing =====================\n");
		 write_log (\$log_file, \"Copied  Gene level calls from $gene_level_calls to  $pwd/$new_proj");


		 #find the ids that correspond to this project
		 my @cmd_proj_ids   =  find_ids_for_project(\$gene_level_calls);
		 write_log(\$log_file, \"This project has $#cmd_proj_ids CMO IDS:");
		 my $all_ids = join ("\t",@cmd_proj_ids);
		 write_log(\$log_file, \$all_ids);
		 #$hash_log{'CMO_IDS'} = $#cmd_proj_ids;

		 #Change the IDS
		 #my $new_gene_file   		       = $new_proj."/Proj_".$proj."_ccs_filtered_genelevel_cna.txt";
		 my $new_gene_file                     = find_input_files(\$new_proj,\"genelevel_cna");
		 #print "NEW GENEVEL $new_gene_file\n";
		 my ($new_gene_mapped, $changed_genes, $genelevel_count) = change_mutation(\$tumor_mapping, \$new_gene_file, \$output_dir, \$log_file);
		 print "MAPPED : $new_gene_mapped\n";
		 write_log (\$log_file, \"Replaced in $gene_level_calls:CMO IDs with DMP IDs");
		 print "REPLACE ALL IDS $$new_gene_mapped  $$new_gene_mapped $$changed_genes \n";
		 replace_all_ids($changed_genes, $new_gene_mapped); 

		 my %hash2 = compare_CMO_IDs($changed_genes,\@roslin_orig_cmo_ids,\"Gene_level");
		 push (@type, "Gene_level");
		 $hash_log{"Gene_level"} = {%hash2};
		#EXOME MSI FILE
		#Subset the original MSI score files from CMO IDs and then replace with DMP ids
		print "SUBSETTING MSI SCORES: $msi_score\n\n";
                my $new_msi_score            = subset_file(\$msi_score, \@cmd_proj_ids );
		print "NEW MSI SCORE MADE : $new_msi_score\n";
		my ($msi_file, $msi_changed_samples,$msi_count) = change_mutation(\$tumor_mapping, \$new_msi_score, \$output_dir, \$log_file);
		my  %hash3 = compare_CMO_IDs($msi_changed_samples,\@roslin_orig_cmo_ids,\"MSI");
		push (@type, "MSI");
		$hash_log{"MSI"} = {%hash3};

		 #FACETS PURITY SEGMENTATION FILE
		 `cp $purity_seg  $new_proj`;
		  my $time                                   = getLogTime();
		  print "$time:Copied purity facets calls : $purity_seg to $input_dir\n";
		  write_log (\$log_file, \"=================== FACETS Purity segment call processing =====================\n");
		  write_log (\$log_file, \"$time:Copied  segment file from $purity_seg to  $pwd/$new_proj");

		  #Change the IDS
		  my $time                                   = getLogTime();
		  my $new_purityseg_file                     = find_input_files(\$new_proj,\"purity_seg");
		  #my $new_purityseg_file                     = $new_proj."/Proj_".$proj."_purity_seg.txt";
		  my ($new_purityseg_mapped, $changed_genes, $purity_count) = change_mutation(\$tumor_mapping, \$new_purityseg_file, \$output_dir, \$log_file);
		  write_log (\$log_file, \"$time:Replaced in $purity_seg:CMO IDs with DMP IDs");
		  my %hash4 = compare_CMO_IDs($changed_genes,\@roslin_orig_cmo_ids,\"Purity");
		  push (@type, "Purity");
		  $hash_log{"Purity"} = {%hash4};

		  #FACETS WGD CALLS
		  `cp $wgd_file $new_proj`;
		  my $time                                   = getLogTime();
		  print "$time:Copied WGD calls : $wgd_file to $input_dir\n";

		  #find the ids that correspond to this project
		  write_log (\$log_file, \"=================== FACETS WGD call  processing ==========================\n");
		  write_log (\$log_file, \"$time:Copied ALL WGD calls from $wgd_file to  $pwd/$new_proj");

		  #Change the IDS
		  my $new_wgd_file		      = $new_proj."/WES_facets_estimates_WGD.v2.txt";
		  my $new_wgd_changed_file            = subset_file(\$new_wgd_file, \@cmd_proj_ids );
		  my ($new_wgd_mapped, $changed_genes)= change_mutation(\$tumor_mapping, \$new_wgd_changed_file, \$output_dir, \$log_file);
		  my %hash5 = compare_CMO_IDs($changed_genes,\@roslin_orig_cmo_ids,\"WGD");
		  push (@type, "WGD");
		  $hash_log{"WGD"} = {%hash5};
		 #SIGNATURE CALLS
		 `cp $signature_calls $new_proj`;
		 my $time                                   = getLogTime();
		 print "$time:Copied signature for each project : $signature_calls to $input_dir\n";
		 write_log (\$log_file, \"=================== Mutational Signature call  processing ==========================\n");
		 write_log (\$log_file, \"$time:Copied mutational signature calls from $signature_calls to  $pwd/$new_proj");

		#Change the IDs
		#my $new_signature_file                     = $new_proj."/Proj_".$proj."_fil_top3list_022119.txt";
		my $new_signature_file                      = find_input_files(\$new_proj,\"fil_top3list");
		my ($new_signature_mapped, $changed_genes) = change_mutation(\$tumor_mapping, \$new_signature_file, \$output_dir, \$log_file);
		my %hash6 = compare_CMO_IDs($changed_genes,\@roslin_orig_cmo_ids,\"Signature");	
		#put_dummy_variables($new_signature_mapped, \$output_dir);
		push (@type, "Signature");
		$hash_log{"Signature"} = {%hash6};


		my $cmd_r = run_portal_script(\$proj, $new_mut_mapped,$new_gene_mapped, $new_purityseg_mapped, $new_wgd_mapped, $new_signature_mapped, $msi_file, \$new_sample_file, \$new_patient_file, \$output_dir);

		print Dumper \%hash_log;


		        print LOG2 "PROJECT_ID\t";
			print LOG2 "CMO_ID\t";
			print LOG2 join ("\t", @type);
			print LOG2 "\n";
			foreach my $keys (sort @roslin_orig_cmo_ids){
				print LOG2 "LOG::$proj\t$keys\t";
				for(my $i =0; $i<=$#type; $i++){
					print LOG2 "$hash_log{$type[$i]}{$keys}\t";
		#ddMSI scores
			}
			print LOG2 "\n";
		}

		write_log(\$log_file, \"==================== Rscript command =========================\n");
		write_log(\$log_file, \"Rscript command : $cmd_r");


	}else{

		print "\nDoes not exists $mutation_file or $gene_level_calls or $purity_seg\n";
	}

}
close(LOG2);


#Read the sample_pairing file from the  roslin inputs folder
sub find_ids_from_original{
	my $proj   = shift;
	my $new_ar = shift;
	my @cmo_orig_paths;
	my @cmo_tumor_ids; 

		

	chomp $$proj;
	my $new_proj = "Proj_".$$proj;
	foreach my $line(@$new_ar){
		chomp $line;
		my @a = split("\t",$line);
		chomp $a[0];
		if ($new_proj eq $a[0]){
			my $pairing = $a[4]."/inputs/*_sample_pairing.txt";
			my @path = qx{ ls $pairing}; chomp $path[0]; 
			if (-e $path[0]){
				#print "CHECK :$path[0]\n";
				push(@cmo_orig_paths,$path[0]);
				}
			
		}
	}

	#UNIQUE IT 
	my @paths = uniq(@cmo_orig_paths);
	my @sample_pairing = read_file(\$paths[0]);

	#Change this code
	foreach my $i (@sample_pairing){
		chomp $i;
		my @a = split("\t", $i);
		push (@cmo_tumor_ids, $a[1]);
	}

	return (@cmo_tumor_ids);
}

#Subset data_sample and data_clinical files
#Rememeber the sample ids differ
sub subset_data_clinical{
	my $dfile    = shift;
	my $dmp_ids  = shift;
	my $type     = shift;
	my $out      = shift;

	
	my $out_file = $$out."/data_".$$type;	
	print "Original files : $$dfile\n";
	print "Making new data_clinical_sample and patient files: $out_file\n";
	open(WRITE, ">$out_file");
	foreach my $dmp (@$dmp_ids){
		chomp $dmp;
		my @a      = split("-", $dmp);
		my $substr = $a[0]."-".$a[1]."-".$a[2];
	

		my @dcfile = read_file($dfile);
		#Look through the data clinical file
		foreach my $data( @dcfile){
			chomp $data;
			if ($$type eq "patient"){
				if (($data =~ m/^PATIENT_ID/) || ($data =~ m/^#/)){
				       print WRITE "$data\n";
				}

			}elsif($$type eq "sample") {
				if (($data =~ m/^SAMPLE_ID/) || ($data =~ m/^#/)){ 
					print WRITE "$data\n"; 			
				}
			}

			#redo this
			if (($data =~ m/PATIENT_ID/) || ($data =~ m/SAMPLE_ID/) || ($data =~ m/^#/)){
				#print WRITE "$data\n";	
			}else{
				my $new_substr =  _substr_(\$substr, $type);
				my @dc         = split("\t", $data);
				my $id_match   = _substr_(\$dc[0],$type);
				
				if ($id_match eq $new_substr){
					print WRITE  "$data\n";
				}
			}

		}

	}
	close(WRITE);
	return $out_file;
}

sub _substr_{
	my $id   = shift;
	my $type = shift;

	if( $$type eq "sample"){
		my ($id1, $id2, $id3, $id4) = split("-",$$id);
		return ($id1."-".$id2."-".$id3);
	}
	
	if ($$type eq "patient"){
		my ($id1, $id2, $id3, $id4) = split("-",$$id);
		return ($id1."-".$id2);
	}
}


#run the portal script
sub run_portal_script{
	my $proj        = shift;
	my $mut_file    = shift;
	my $gene_file   = shift;
	my $purity_file = shift;
	my $wgd_file    = shift;
	my $sig_file    = shift;
	my $msi_file    = shift;
	my $sample_f	= shift;
	my $patient_f	= shift;
	my $out         = shift;

	my $cmd = "Rscript  port_script_jon.R $$proj $$gene_file $$mut_file  $$wgd_file  $$purity_file $$sig_file $$msi_file $$sample_f $$patient_f  $$out" ;
	print "$cmd\n";
        `$cmd`;	

	return $cmd;
}

#Pick the right seg file depending on the 

#Write into a log file
sub write_log{
	my $log     = shift;
	my $message = shift;

 	open(WRITE, ">>$$log");
	print WRITE "$$message\n";

}

#Make a log file
sub mkdir_log{
	my $log   = shift;
	my $proj  = shift;

	my $log_file = $$log."/Proj_".$$proj."_LOG_FILE";
	open (WRITE, ">>$log_file");
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	print WRITE "Running Proj_$$proj at $mday/$mon/$year:$hour:$min:$sec\n";

	return $log_file;
}

#Call to make directories
sub make_directories{
	my $dir_l = shift;
	my $dir_i = shift;
	my $dir_o = shift;
	
	if (! mkdir_proj($dir_l)){ print "Something went wrong in creating log directory $$dir_l\n"};
	if (! mkdir_proj($dir_i)){ print "Something went wrong in creating input directory $$dir_i\n"};
	if (! mkdir_proj($dir_o)){ print "Something went wrong in creating output directory $$dir_o\n"};
}


sub define_input_files{
	my $delivery = shift;
	my $proj     = shift;
	my $file     = shift;

	my $new_file;

	if ( $$file =~ m/fil_top3list_/){
		#$new_file = $$delivery.$$proj."_".$$file;
		$new_file = $$delivery;
            if (-e "$new_file"){
		 	print "Exists input file : $new_file\n";
		}else{

			print "$new_file file  does not exists\n";
			exit;
		}
	}elsif($$file =~ m/genelevel_cna/){
	 	$new_file = $$delivery."/Proj_".$$proj."_".$$file;

		
	}else{
		if (-e "$$delivery$$proj"){
	
			$new_file = $$delivery.$$proj."/Proj_".$$proj."_".$$file;
			print "Exists input file : $new_file\n";
		}else{
			print "$$delivery$$proj  directory does not exists\n";
			exit;
		}
	}
      return $new_file;
}

sub mkdir_proj{
 	my $dir = shift;

	`mkdir -p $$dir`;
	if ( -e $$dir){
		return 1;
	}else{
		return 0;
	}	
}



#Extract DMP Ids from change_mutation function to subset data_clinical and patient.
sub extract_dmpid{
	my $string = shift;
	my @all    = split("\n",$$string);
	my @all_dmp_ids;

	foreach my $a(@all){
		chomp $a;
		my ($cmo_id, $dmp_id) = split("\t", $a);
		push (@all_dmp_ids,$dmp_id);
	}
	return (\@all_dmp_ids);
}

sub change_mutation{
	my $tumor_file = shift;
	my $mut_file   = shift;
	my $output_dir = shift;	
	my $log_file   = shift;

	my @replaced_values;
	my %hash   =();


	#Replace the ids in the mutation file
	my $filename  = basename ($$mut_file);
	my @tumor_ids = read_file ($tumor_file);	
	my $new_mapped_file = $$output_dir."/".$filename."_MAPPED";
	open(WRITE, ">$new_mapped_file");
	open(FILE, $$mut_file);
	while(<FILE>){
		chomp;
		if (($_ =~ m/Tumor/) || ($_ =~ m/Sample_ID/)  || ($_ =~ m/CMO/)){
			#$header =$_;
			print WRITE "$_\n";
		}else{

		my @a = split("\t", $_);
		foreach my $h(@tumor_ids){
		      chomp $h;
		      my ($wes, $tumor) = split("\t", $h);
		      #IF ITS A GENELEVEL CALL THEN CHANGE JUST 1ST WORD
		      if ($$mut_file =~ m/final_comb_ccs/){
			 $a[15] =~ s/^\s*(.*?)\s*$/$1/;
		      	if ($a[15] eq $tumor){
				#print $a[15]."\n";
			 s/$a[15]/$wes/g;
			 #s/$a[16]/$wes/g;
			 print WRITE $_."\n";
			 push (@replaced_values, $a[15]);
			 $a[15] =~ s/^\s*(.*?)\s*$/$1/;
			 if ($a[15] =~ m/s_/){
			 	$hash{$a[15]} = $wes;
			 }
		         }
		       }elsif($$mut_file =~ m/genelevel_cna/){
			  $a[0] =~ s/^\s*(.*?)\s*$/$1/;
			  #print $a[1]."--";
			  if ($a[0] eq $tumor){
				s/$a[0]/$wes/g;
				#my $new_str = make_string(\$_,\$wes); 
			   	print WRITE $_."\n";
				$a[0] =~ s/^\s*(.*?)\s*$/$1/;
				push (@replaced_values, $a[0]);
				if ($a[0] =~ m/s_/){
					$hash{$a[0]} = $wes;
				}
		          }		  
		       }elsif($$mut_file =~ m/purity_seg/){
			        $a[0] =~ s/^\s*(.*?)\s*$/$1/;
		       		if ($a[0] eq $tumor){
			           s/$a[0]/$wes/g;
				   print WRITE $_."\n";
				   push (@replaced_values, $a[0]);
				   $a[0] =~ s/^\s*(.*?)\s*$/$1/;
				   if ($a[0] =~ m/s_/){
					   $hash{$a[0]} = $wes;
			           }
				}
		       }elsif ($$mut_file =~ m/estimates_WGD.v2.txt/){
			      
			      #Add header in only this file
			      $a[0] =~ s/^\s*(.*?)\s*$/$1/;
			      if ($a[0] eq $tumor){
				 s/$a[0]/$wes/g;
			         print WRITE $_."\n";
				 $hash{$a[0]} = $wes;
				 push (@replaced_values, $a[0]);
                              }
		       }elsif($$mut_file =~ m/fil_top3list/){
			       $a[0] =~ s/^\s*(.*?)\s*$/$1/;
			       if  ($a[0] eq $tumor){
				s/$a[0]/$wes/g;
				print WRITE $_."\n";
				$hash{$a[0]} = $wes;
				push (@replaced_values, $a[0]);
			       }
		       
		      
	       		}elsif($$mut_file =~ m/exome_msi/){
				$a[0] =~ s/^\s*(.*?)\s*$/$1/;
				#Print the header

				if  ($a[0] eq $tumor){
					s/$a[0]/$wes/g;
					print WRITE $_."\n";
					$hash{$a[0]} = $wes;
					push (@replaced_values, $a[0]);
				}	
			}else{ }
		      }
		}
		}

close(WRITE);
close(FILE);

#Make uniq and make a string of all the replaced values
my @uniq_array = uniq(@replaced_values);
my $string_ids = join ("\n", @uniq_array);
my $count      = scalar(@uniq_array);
my @arr_write;
foreach my $keys (sort keys %hash){
	my $l = "$keys\t$hash{$keys}";
	push(@arr_write,$l);

}
my $string = join("\n",@arr_write);

#write into log
write_log($log_file, \"Replaced $count CMO ids in the new file\n");
write_log($log_file,\$string);
return (\$new_mapped_file, \$string, \$count);

}


#Get log time for
sub getLogTime {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $nice_timestamp = sprintf ( "%04d%02d%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec);
	return $nice_timestamp;
}

#Subset files  when ever needed by supplying an array
sub subset_file{
	my $file =  shift;
	my $ar   =  shift;
	
	my $new_file  = $$file."_subset";
	open(WRITE, ">$new_file");
	my @old_file  = read_file($file);
	foreach my $line(@old_file){
		chomp $line;
		if (($line =~ m/Sample_ID/) || ($line =~ m/CMO/)){
			print  WRITE "$line\n";
		}else{
		my @temp = split("\t", $line);
		foreach my $l (@$ar){
			chomp $l;
			if ($l eq $temp[0]){
				print  WRITE "$line\n";
			}
		}
		}
	}
	close(WRITE);

return $new_file;

}

#Find IDS for  the project
sub find_ids_for_project{
	my $file = shift;
	my @cmo_ids = "";

	my @arr  = read_file($file);
	foreach my $line (@arr){
		chomp $line;
		if ($line =~ m/Tumor_Sample_Barcode/){
		}else{
			my @a = split("\t", $line);
			push (@cmo_ids, $a[0]);
		}
	}

	return (uniq(@cmo_ids));
}

sub replace_all_ids{
	my $ids = shift;
	my $file  = shift;

	my $rand_file = $$file."_HELL_EXCHANGE";
	#OF THE FORMAT 
	#s_C_006550_M001_d       P-0010301-T02-WES
	#s_C_006611_M001_d       P-0011466-T02-WES
	
	my @all_ids = break_each_line(\$$ids);
	foreach my $cmo_ids (@all_ids){
		chomp $cmo_ids;
		my ($cmo, $wes) = split("\t", $cmo_ids);
		my $cmd = "sed -i  \'s/$cmo/$wes/g\' $$file ";
		print "$cmd\n";
		`$cmd`;

	}
}


#Compare the Ids from tumor_id_mapping and latest roslin run
sub compare_CMO_IDs{
	my $string_proj = shift;
	my $arr_ros     = shift;
	my $type        = shift;

	my %hash;
	my @arr_proj = break_each_line($string_proj);

	foreach my $ros_id(@$arr_ros){
		  $hash{$ros_id} = 0;
	}

	foreach my $ros_id (@$arr_ros){
		foreach my $proj_id (@arr_proj){
			chomp $proj_id; chomp $ros_id;
			my ($cmo, $dmp) = split("\t",$proj_id);
			if ($ros_id eq $cmo){
				if (exists $hash{$ros_id}){
			    		$hash{$ros_id} = 1;
				}

			}
			     
		}
	}

  return %hash;
}




sub break_each_line{
	my $line = shift;
	my @arr  = split("\n", $$line);

	return @arr;
}

#Read file and return an array 
sub read_file{
	my $file = shift;
	open(FILE, $$file);
	my @arr = <FILE>;
	close(FILE);

return @arr;
}


#Make a new string after replacing ids. TRIAL!!
sub make_string{
	my $line = shift;
	my $wes  = shift;

	my $str = '';
	my @arr = split("\t", $$line); shift @arr;
	#foreach my $element(@arr){ 
	#	chomp $element;
		$str = join ("\t", @arr);
		my $new_str = $$wes."\t".$str;
		if ($new_str =~ m/BC/){
			print $new_str."\n";
		}
		#}
	return $new_str;
}
#Read the JSON file as a config file
sub read_config_file{
	my $file = shift;  chomp $$file;
	my $json_file = json_file_to_perl($$file);
	return ($json_file); 
}

#Redacted list;
sub redaction_list{
	my $proj   = shift;
	my $file   = shift;

	my @redact = read_file($file);
	my @redaction = ();
	foreach my $line(@redact){
		chomp $line;
		my ($proj_r, $cmo_r, $wes_r,$wes_r2) = split("\t", $line);

		if ($$proj eq $proj_r){
			push (@redaction, $cmo_r);
		}
	
	}
	print "Project REDACTED:  $$proj  @redaction\n";

	return  @redaction;

}

sub remove_redactions{
	my $redact_arr = shift;
	my $file       = shift;

	my @contents =  read_file($file);

	for my $cmo_id (@$redact_arr){ 
		chomp $cmo_id;
		#@$redacted_array = grep ! /$cmo_id/, @$redacted_arr;
	}

	#print @$redacted_array;
}

sub find_input_files{
	my $input = shift;
	my $type  = shift;

	my @a = qx{ ls $$input/*$$type*};
	chomp $a[0];
	return $a[0];
}
