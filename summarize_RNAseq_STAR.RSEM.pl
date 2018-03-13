#!/usr/bin/perl

use strict;
use warnings;

#script to report RNA-seq summary from outputs after running Ting's pipeline "run_RNAseq_STAR.RSEM.sh"
#Ting Wang
#Jan 2017

if(@ARGV != 2){
	die "wrong inputs\nUSAGE :: perl summarize_RNAseq_RSEM.pl <directory> <outfile>\n";
}else{
	open my $outfile, ">$ARGV[1]" or die "wrong outfile\n";

	print $outfile "Sample\tTotal_Reads\tTrimmed_Reads\tRemaining_Reads\tUniquely-mapped_Rate\tMulti-mapped_Rate\tUnmapped_Rate\tDuplicate_Rate\tTotal_Assigned_Tags\tCDS_Exon_Tags\t5'_UTR_Exon_Tags\t3'_UTR_Exon_Tags\tIntron_Tags\tTSS_Up10kb\tTES_Down10kb\tSample_InTable\t#Gene(count>0)\tTotal_Count(>0)\t#Gene_Locus(FPKM>0)\tTotal_FPKM(>0)\n";
	
	my $dir = $ARGV[0];

	#17 and #18 -- get count info
	open my $counts, "$dir/counts.txt " or die "counts file error\n";
	
	my @tot_count;
	my @gene_count;
	while(<$counts>){
		if($. > 1){
			chomp;
			my @line = split /\s+/;

			for(my $i = 1; $i < @line ; $i++){
				$tot_count[$i-1] += $line[$i];
			
				if($line[$i] > 0){
					$gene_count[$i-1]++;
				}
			}
		}
	}

	#get fpkm info
	open my $fpkm, "$dir/rsem_fpkm_genes.txt " or die "counts file error\n";
			
	my @num_gene_locus;
	my @tot_fpkm;
	my @sample_name;
	while(<$fpkm>){
		if($. == 1){
			chomp;
			@sample_name = split /\s+/;
			shift @sample_name;
			shift @sample_name;
		}
		if($. > 1){
			chomp;
			my @line = split /\s+/;

			for(my $i = 2; $i < @line ; $i++){
				$tot_fpkm[$i-2] += $line[$i];
			
				if($line[$i] > 0){
					$num_gene_locus[$i-2]++;
				}
			}
		}
	}

	#do rest of stuff
	opendir(DIR, $dir) or die $!;
	my $sample_num = 0;
	foreach my $file (sort readdir(DIR)){
		next unless (-d "$dir/$file");

		#tried "ne" and it didn't work (?)
		if($file eq "." || $file eq ".."){}else{
			#1 name
			print $outfile "$file\t";

			#2 and #3 and #4 total reads and trimmed reads and remaining
			opendir(SUBDIR, "$dir/$file");
			my @trim_report_file = grep(/trimming_report.txt$/, sort readdir(SUBDIR));
			if(@trim_report_file == 1){
				open my $trim_report, "$dir/$file/$trim_report_file[0]" || die "trimming_report.txt open error\n";
				my $numlines = () = <$trim_report>;
				seek($trim_report, 0, 0);

				my $tot_reads = 0;
				my $trimmed_reads = 0;

				my $lc = 1;
				while(<$trim_report>){
					if($lc == $numlines - 2){
						my @line = split /\s+/;
						$tot_reads = $line[0];
					}

					if($lc == $numlines -1 ){
						my @line = split /\s+/;
						$trimmed_reads = $line[-2];
					}	
					$lc++;
				}						
				my $remain_reads = $tot_reads - $trimmed_reads;
				print $outfile "$tot_reads\t$trimmed_reads\t$remain_reads\t";	
				close($trim_report);
			}elsif(@trim_report_file == 2){
				open my $trim_report, "$dir/$file/$trim_report_file[1]" || die "trimming_report.txt open error\n";
				my $numlines = () = <$trim_report>;
				seek($trim_report, 0, 0);
			
				my $tot_reads = 0;
				my $trimmed_reads = 0;

				my $lc = 1;
				while(<$trim_report>){
					if($lc == $numlines - 4){
						my @line = split /\s+/;
						$tot_reads = $line[0];
					}
					if($lc == $numlines ){
						my @line = split /\s+/;
						$trimmed_reads = $line[-2];
					}	
					$lc++;
				}

				my $remain_reads = $tot_reads - $trimmed_reads;						
				print $outfile "$tot_reads\t$trimmed_reads\t$remain_reads\t";	
				close($trim_report);
			}else{
				die "number of trimming_report.txt file is not 1 (if single-end data) or 2 (if paired-end data)\n";
			}

			#5, #6, #7 -- uniquely-mapped, multi-mapped, unmapped
			open my $log_final, "$dir/$file/Log.final.out" or die "log.final.out error\n";
			my $line25 = 0;
			my $line27 = 0;

			my $line29 = 0;
			my $line30 = 0;
			my $line31 = 0;
			while (<$log_final>){
				#5
				if($. == 10){
					my @line = split /\s+/;
					print $outfile "$line[6]\t";
				}

				#6
				if($. == 25){
					my @line = split /\s+/;
					$line25 = $line[9];
				}

				if($. == 27){
					my @line = split /\s+/;
					$line27 = $line[10];
                                        
					$line25 =~ s/%//g;
					$line27 =~ s/%//g;

					my $mult_map = $line25 + $line27;
					print $outfile "$mult_map"."%\t";
				}

				#7
				if($. == 29){
					my @line = split /\s+/;
					$line29 = $line[9];
				}

				if($. == 30){
					my @line = split /\s+/;
					$line30 = $line[8];
				}

				if($. == 31){
					my @line = split /\s+/;
					$line31 = $line[7];

					$line29 =~ s/%//g;
					$line30 =~ s/%//g;
					$line31 =~ s/%//g;

					my $unmapped  = $line29 + $line30 + $line31;
					print $outfile "$unmapped"."%\t"
				}
			}

			#8 dup rate
			open my $rmdup_file, "$dir/$file/Aligned.sortedByCoord.out.dedupped.flagstat" or die "rmdup file error\n";
			my $read_numb = 0;
			my $dup_numb = 0;
			while(<$rmdup_file>){
				chomp;
				my @line = split /\s+/;
				if($. == 1){
					$read_numb = $line[0];
				}
				if($. == 4){
					$dup_numb = $line[0];
				}
			}
			my $dup_rate = ($dup_numb / $read_numb) * 100;
			my $rounded = sprintf("%.2f",$dup_rate);
			print $outfile "$rounded"."%\t";

			#9-14, total assigned tags, other stuff
			open my $dist, "$dir/$file/Aligned.sortedByCoord.out.bam.distribution" or die "dist file error\n";
			my $tags;
			while(<$dist>){
				chomp;
				my @line = split /\s+/;

				if($. == 3){
					$tags = $line[3];
					print $outfile "$tags\t";
				}
				
				if($. == 6 || $. == 7 || $. == 8 || $. == 9 || $. == 12 || $. == 15){
					my $tmp = $line[2] / $tags * 100;
					my $rounded = sprintf("%.2f",$tmp);
					print $outfile "$rounded"."%\t";
				}
			}

			#16 and #17, # gene and total count
			print $outfile "$sample_name[$sample_num]\t";
			print $outfile "$gene_count[$sample_num]\t$tot_count[$sample_num]\t";	

			#18 and #19 gene locus and fpkm
			print $outfile "$num_gene_locus[$sample_num]\t";
			my $tmp_tot_fpkm = sprintf("%.2f",$tot_fpkm[$sample_num]);
			print $outfile "$tmp_tot_fpkm";

			#all done
			print $outfile "\n";
			$sample_num++;
		}
	}

	closedir(DIR);
}

#transform matrix
open(FD, "$ARGV[1]") || die;
my $line = <FD>;
my @out = split(/\s+/, $line);
while($line = <FD>){
	my @arr = split(/\s+/, $line);
	for (my $i = 0; $i < @arr; $i++){
		$out[$i].="\t".$arr[$i];
	}
}
close(FD);
open(FD,">$ARGV[1].tmp") || die;
for (my $i = 0; $i < @out; $i++){
	print FD $out[$i]."\n";
}
close(FD);
system("mv $ARGV[1].tmp $ARGV[1]");
exit 0;

