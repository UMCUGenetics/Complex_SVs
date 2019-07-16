#! /usr/bin/perl

use strict;

my $delly_folder = shift; # this should be the link to the folder of Delly output
my $run_id = shift;
my $output_folder = shift;

my $work_folder = $output_folder."Raw/".$run_id."/";

print($work_folder);

mkdir($work_folder);

my @SV_types = ("DUP", "INS", "INV","TRA");

my $outputs = "";

# first remove the headers (except for deletions, this file is used for the headers)
foreach my $n (@SV_types) {

	my $output = join "",$work_folder, $run_id, "_", $n, "_noheader.vcf";
	my $remove_header = join "", "tail -n +37 ",$delly_folder, $run_id, "_", $n, ".vcf > ", $output;
	
	print $remove_header, "\n";
	
	$outputs = join " ", $output, $outputs;
	system($remove_header);
	
}

#print $outputs, "\n";

# merge the VCFs without headers

my $delly_del = join "", $delly_folder, $run_id, "_DEL.vcf";
my $merge_vcfs = join "", "cat ", $delly_del, " ", $outputs, "> ", $work_folder, $run_id, "_Delly_merged.vcf";
print $merge_vcfs, "\n";
system($merge_vcfs);

my $delly = join "", $work_folder, $run_id,"_Delly_merged.vcf";

# split the VCFs by sample

my @column_names;
open DELLY, $delly;
while(<DELLY>) {
  
  if ($_ =~ /^#CHROM/) {
  	@column_names = split/\t/, $_;

  	#@samples = @column_names[9...25];
  }
}

my $number_columns = (scalar @column_names) - 1;

#print join ",", @column_names, "\n";
#print join "", "Number of columns: ", scalar @column_names, "\n";

close DELLY;

foreach my $column (9...$number_columns){
	
	my $sample = (split /_/, @column_names[$column])[0];
	print $sample, "\n";

	my $output_file = $output_folder.$sample."_delly.vcf";
	
	print $output_file, "\n";
	
	open DELLY, $delly;
	open OUT, ">", $output_file;
		
  	while(<DELLY>){
  		if ($_ =~ /^##/) {
  			print OUT $_;	
  		} elsif ($_ =~ /^#CHROM/){
  			my @header = split/\t/, $_;
  			foreach my $column_name (@header[0...8]){
  				print OUT $column_name, "\t"; 
  			}
  			foreach my $column_name (@header[$column]){
  				print OUT $column_name; 
  			}
  			
  			if ($column != $number_columns){
  				print OUT "\n";
  			}
  		} else {
  				
  		my @row = split/\t/, $_;
  		
  		my $genotype = @row[$column];
  	 	next if $genotype =~ /^0\/0/ ; # skip the row if the genotype is 0/0
    
    	next if $genotype =~ /^\.\/\./ ; # skip the row if the genotype is ./.
  		
  		foreach my $cell (@row[0...8]){
			print OUT $cell, "\t";
  		}
  		foreach my $cell (@row[$column]){
  			print OUT $cell;
  		}
  		if ($column != $number_columns){
  			print OUT "\n";
  			
  		}
  		
  		} 					
	}
	close DELLY;
	close OUT;
		
}





