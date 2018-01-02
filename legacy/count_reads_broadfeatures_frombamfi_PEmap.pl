use warnings;
use strict;

## Note 2015/11/12 - this script takes in ucsc table format data from gencode (now manually generated by me based off of gencode GTF files), which are 0-based open ended
# &read_gencode converts those features to 1-based, closed ended features for overlaps with reads from sam format (also 1-based, closed ended)

my $verbose_flag = 0;
my $hashing_value = 10000;

my %all_features;
my %enst2type;
my %enst2ensg;
my %ensg2name;
my %ensg2type;

my $clip_fi1 = $ARGV[0];
my $output_fi = $ARGV[1];

#defaults to hg19
my $species = "hg19";
if (exists $ARGV[2]) {
    $species = $ARGV[2];
}

if (-e $output_fi) {
    print STDERR "$output_fi already exists, skipping...\n";
    exit;
}

open(OUT,">$output_fi");
my $gencode_gtf_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
#my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencodev19_comprehensive";

if ($species eq "mm9") {
    $gencode_gtf_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf.parsed_ucsc_tableformat";
}
&read_gencode_gtf($gencode_gtf_file);
&read_gencode($gencode_tablebrowser_file);

my %readcounts_by_region;


#my $clip_fi1 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/IGF2BP1-205-INPUT_S4_L001_R2.polyATrim.adapterTrim.rmRep.rmDup.sorted.bam.5.pos.bg";
#my $clip_fi2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/IGF2BP1-205-INPUT_S4_L001_R2.polyATrim.adapterTrim.rmRep.rmDup.sorted.bam.5.neg.bg";
#my $clip_fi1 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/IGF2BP1_HepG2_205_02.bam.5.pos.bg";
#my $clip_fi2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/IGF2BP1_HepG2_205_02.bam.5.neg.bg";

#my $clip_fi1 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/RBFOX2_HepG2_204_02.bam.5.pos.bg";
#my $clip_fi2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/RBFOX2_HepG2_204_02.bam.5.neg.bg";
#my $clip_fi1 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/RBFOX2-204-INPUT_S2_L001_R2.polyATrim.adapterTrim.rmRep.rmDup.sorted.bam.5.pos.bg";
#my $clip_fi2 = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/RBFOX2-204-INPUT_S2_L001_R2.polyATrim.adapterTrim.rmRep.rmDup.sorted.bam.5.neg.bg";


&read_bam_fi($clip_fi1);

#&read_bgfi($clip_fi1,"+");
#&read_bgfi($clip_fi2,"-");



my @final_feature_types = ("CDS","5utr","3utr","5utr|3utr","intron","intergenic","noncoding_exon","noncoding_intron");
print OUT "ENSG";
for my $feature_type (@final_feature_types) {
    print OUT "\t$feature_type";
}
print OUT "\n";

for my $key (keys %readcounts_by_region) {
    print OUT "$key";
    for my $feature_type (@final_feature_types) {
	my $toprint = 0;
	if (exists $readcounts_by_region{$key}{$feature_type}) {
	    $toprint = $readcounts_by_region{$key}{$feature_type};
	}
	print OUT "\t$toprint";
    }
    print OUT "\n";
    
}


sub read_bam_fi {
    my $bamfile = shift;

    print STDERR "now doing $bamfile\n";
    if ($bamfile =~ /\.bam/) {
        open(B,"samtools view -h $bamfile |") || die "no $bamfile\n";
    } elsif ($bamfile =~ /\.sam/) {
        open(B,$bamfile) || die "no sam $bamfile\n";
    } else {
        print STDERR "file format error not .sam or .bam \n";
        exit;
    }
    while (<B>) {
        my $r1 = $_;

	$verbose_flag = 0;

        if ($r1 =~ /^\@/) {
            next;
        }
        my @tmp_r1 = split(/\t/,$r1);
        my @read_name = split(/\:/,$tmp_r1[0]);
        my $randommer = $read_name[0];

        my $r1_cigar = $tmp_r1[5];
        my $r1sam_flag = $tmp_r1[1];
        my $mismatch_flags = $tmp_r1[14];

        my $r1_chr = $tmp_r1[2];
        my $r1_start = $tmp_r1[3];

        my $frag_strand;
	if ($r1sam_flag == 147) {
#        if ($r1sam_flag == 16) {
            $frag_strand = "-";
        } elsif ($r1sam_flag == 163) {
            $frag_strand = "+";
        } else {
	    next;
            print STDERR "R1 strand error $r1sam_flag\n";
        }

	my @read_regions = &parse_cigar_string($r1_start,$r1_cigar,$r1_chr,$frag_strand);

	my $feature_flag = 0;
        my %tmp_hash;
        for my $region (@read_regions) {
	    
	    my ($rchr,$rstr,$rpos) = split(/\:/,$region);
            my ($rstart,$rstop) = split(/\-/,$rpos);
	    
            my $rx = int($rstart / $hashing_value);
            my $ry = int($rstop  / $hashing_value);
            for my $ri ($rx..$ry) {
		
		
		for my $feature (@{$all_features{$rchr}{$rstr}{$ri}}) {
		    my ($feature_enst,$feature_type,$feature_region) = split(/\|/,$feature);
		    my ($feature_start,$feature_stop) = split(/\-/,$feature_region);

#		    print STDERR "feature $feature region $region $_\n" if ($rstart == $feature_stop);
#		    print STDERR "feature $feature region $region $_\n" if ($rstop == $feature_start);

		    next if ($rstop < $feature_start);
		    next if ($rstart > $feature_stop);
		    my $overlap_start = &max($feature_start,$rstart);
		    my $overlap_stop  = &min($feature_stop,$rstop);
		    
		    


#		    if ($feature_start <= $rstart && $rstop <= $feature_stop) {
			# read start is within region                                                                                     
			
#                   if ($feature_start == $rstart_i) {                                                                                
#			if ($feature_stop == $rstop) {
#			if ($feature_start == $rstart) {      
#			    print "fi $bamfile r1 $r1 region $region strand $rstr\nrstarti $rstart $_ feature $feature\n";
#			}
                             

# If any part of read is intronic -> call intronic for that transcript
# else - use whatever majority of read is in
          
		    my $feature_ensg = $enst2ensg{$feature_enst};
#		    print "rstart $rstart rstop $rstop featurestart $feature_start featurestop $feature_stop overlapstart $overlap_start overlapstop $overlap_stop\n";
		    $tmp_hash{$feature_ensg}{$feature_enst}{$feature_type}= ($overlap_stop - $overlap_start+1)/($feature_stop-$feature_start+1);

		    $feature_flag = 1;
#		    $verbose_flag = 1 if ($feature_ensg eq "ENSG00000140988.11");       
#                   $verbose_flag = 1 if ($feature_ensg eq "ENSG00000169683.3");                                                      
		    print "found feature $_ $feature $feature_ensg $feature_type\n" if ($verbose_flag == 1);
		}
	    }
	}
	
	if ($feature_flag == 1) {
	    for my $feature_ensg (keys %tmp_hash) {
		
		my %enst_feature_type_hash;
		for my $feature_enst (keys %{$tmp_hash{$feature_ensg}}) {
		    my $enst_feature_type = "intergenic";
		    if (exists $tmp_hash{$feature_ensg}{$feature_enst}{"intron"}) {
			$enst_feature_type = "intron";
		    } elsif (exists $tmp_hash{$feature_ensg}{$feature_enst}{"5utr"} || exists $tmp_hash{$feature_ensg}{$feature_enst}{"CDS"} || exists $tmp_hash{$feature_ensg}{$feature_enst}{"3utr"}) {
			$enst_feature_type = &max_from_hash(\%{$tmp_hash{$feature_ensg}{$feature_enst}},"3utr|CDS|5utr");
			print "chose feature by max percent $feature_ensg $feature_enst $enst_feature_type\n" if ($verbose_flag == 1);
		    } elsif (exists $tmp_hash{$feature_ensg}{$feature_enst}{"noncoding_intron"}) {
			$enst_feature_type = "noncoding_intron";
		    } elsif (exists $tmp_hash{$feature_ensg}{$feature_enst}{"noncoding_exon"}) {
			$enst_feature_type = "noncoding_exon";
		    }

		    $enst_feature_type_hash{$enst_feature_type} = 1;
		}

		my $final_feature_type = "intergenic";
		if (exists $enst_feature_type_hash{"CDS"}) {
		    $final_feature_type = "CDS";
		} elsif (exists $enst_feature_type_hash{"3utr"} || exists $enst_feature_type_hash{"5utr"}) {
		    if (exists $enst_feature_type_hash{"3utr"} && exists $enst_feature_type_hash{"5utr"}) {
			$final_feature_type = "5utr|3utr";
		    } elsif (exists $enst_feature_type_hash{"3utr"}) {
			$final_feature_type = "3utr";
		    } elsif (exists $enst_feature_type_hash{"5utr"}) {
			$final_feature_type = "5utr";
		    } else {
			print STDERR "weird shouldn't hit this\n";
		    }
		    
		} elsif (exists $enst_feature_type_hash{"intron"}) {
		    $final_feature_type = "intron";
		} elsif (exists $enst_feature_type_hash{"noncoding_exon"}) {
		    $final_feature_type = "noncoding_exon";
		} elsif (exists $enst_feature_type_hash{"noncoding_intron"}) {
		    $final_feature_type = "noncoding_intron";
		}
		
		print "adding to $feature_ensg $final_feature_type \n" if ($verbose_flag == 1);
		$readcounts_by_region{$feature_ensg}{$final_feature_type} ++;
		$readcounts_by_region{$feature_ensg}{all} ++;
		$readcounts_by_region{all}{$final_feature_type} ++;
		$readcounts_by_region{all}{all} ++;
	    }
	} else {
	    my $final_feature_type = "intergenic";
	    print "adding to all $final_feature_type \n" if ($verbose_flag == 1);
	    $readcounts_by_region{all}{$final_feature_type} ++;
	    $readcounts_by_region{all}{all} ++;
	}
	
    }
    close(B);
}

sub max_from_hash {
    my $ref = shift;
    my %hashref = %$ref;

    my $desired = shift;
    my @desired_array = split(/\|/,$desired);
    
    my $current_max = "NA";
    my $max_element;
    for my $des (@desired_array) {
	next unless (exists $hashref{$des});
	if ($current_max eq "NA" || $hashref{$des} > $current_max) {
	    $max_element = $des;
	    $current_max = $hashref{$des};
	}
    }
    
    return($max_element);
}

sub max {
    my $x = shift;
    my $y = shift;
    if ($x > $y) {
	return($x);
    } else {
	return($y);
    }
}

sub min {
    my $x = shift;
    my $y = shift;
    if ($x < $y) {
        return($x);
    } else {
        return($y);
    }
}


sub parse_cigar_string {
    my $region_start_pos = shift;
    my $flags = shift;
    my $chr = shift;
    my $strand = shift;

    my $current_pos = $region_start_pos;
    my @regions;

    while ($flags =~ /(\d+)([A-Z])/g) {
       
        if ($2 eq "N") {
            #read has intron of N bases at location
            
            push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);

            $current_pos += $1;
            $region_start_pos = $current_pos;
        } elsif ($2 eq "M") {
            #read and genome match
            $current_pos += $1;
        } elsif ($2 eq "S") {
            #beginning of read is soft-clipped; mapped pos is actually start pos of mapping not start of read
        } elsif ($2 eq "I") {
            #read has insertion relative to genome; doesn't change genome position
        } elsif ($2 eq "D") {
#           push @read_regions,$chr.":".$current_pos."-".($current_pos+=$1);
            $current_pos += $1;
            #read has deletion relative to genome; genome position has to increase
        } else {
            print STDERR "flag $1 $2 $flags\n";

        }
    }
    push @regions,$chr.":".$strand.":".$region_start_pos."-".($current_pos-1);

    return(@regions);
}



sub read_gencode {
    ## eric note: this has been tested yet for off-by-1 issues with ucsc brower table output 5/29/15
    my $fi = shift;
#    my $fi = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencodev19_comprehensive";
    print STDERR "reading in $fi\n";
    open(F,$fi);
    while (<F>) {
	chomp($_);
	my @tmp = split(/\t/,$_);
	my $enst = $tmp[1];
	next if ($enst eq "name");
	my $chr = $tmp[2];
	my $str = $tmp[3];
	my $txstart = $tmp[4];
	my $txstop = $tmp[5];
	my $cdsstart = $tmp[6];
	my $cdsstop = $tmp[7];
	
	my @starts = split(/\,/,$tmp[9]);
	my @stops = split(/\,/,$tmp[10]);
	
	my @tmp_features;

	my $transcript_type = $enst2type{$enst};
	unless ($transcript_type) {
	    print STDERR "error transcript_type $transcript_type $enst\n";
	}
	if ($transcript_type eq "protein_coding") {
	
	    for (my $i=0;$i<@starts;$i++) {
		if ($str eq "+") {
		    if ($stops[$i] < $cdsstart) {
			# exon is all 5' utr
			push @tmp_features,$enst."|5utr|".($starts[$i]+1)."-".$stops[$i];
		    } elsif ($starts[$i] > $cdsstop) {
			#exon is all 3' utr
			push @tmp_features,$enst."|3utr|".($starts[$i]+1)."-".$stops[$i];
		    } elsif ($starts[$i] > $cdsstart && $stops[$i] < $cdsstop) {
			#exon is all coding
			push @tmp_features,$enst."|CDS|".($starts[$i]+1)."-".$stops[$i];
		    } else {
			my $cdsregion_start = $starts[$i];
			my $cdsregion_stop = $stops[$i];
			
			if ($starts[$i] <= $cdsstart && $cdsstart <= $stops[$i]) {
			    #cdsstart is in exon
			    my $five_region = ($starts[$i]+1)."-".$cdsstart;
			    push @tmp_features,$enst."|5utr|".$five_region;
			    $cdsregion_start = $cdsstart;
			}
			
			if ($starts[$i] <= $cdsstop && $cdsstop <= $stops[$i]) {
			    #cdsstop is in exon
			    my $three_region = ($cdsstop+1)."-".$stops[$i];
			    push @tmp_features,$enst."|3utr|".$three_region;
			    $cdsregion_stop = $cdsstop;
			}
			
			my $cds_region = ($cdsregion_start+1)."-".$cdsregion_stop;
			push @tmp_features,$enst."|CDS|".$cds_region;
		    }
		} elsif ($str eq "-") {
		    if ($stops[$i] < $cdsstart) {
			# exon is all 5' utr                                                                                                  
			push @tmp_features,$enst."|3utr|".($starts[$i]+1)."-".$stops[$i];
		    } elsif ($starts[$i] > $cdsstop) {
			#exon is all 3' utr                                                                                                   
			push @tmp_features,$enst."|5utr|".($starts[$i]+1)."-".$stops[$i];
		    } elsif ($starts[$i] > $cdsstart &&$stops[$i] < $cdsstop) {
			#exon is all coding                                                                                                   
			push @tmp_features,$enst."|CDS|".($starts[$i]+1)."-".$stops[$i];
		    } else {
			my $cdsregion_start = $starts[$i];
			my $cdsregion_stop = $stops[$i];
			
			if ($starts[$i] <= $cdsstart && $cdsstart <= $stops[$i]) {
			    #cdsstart is in exon                                                                                              
			    my $three_region = ($starts[$i]+1)."-".$cdsstart;
			    push @tmp_features,$enst."|3utr|".$three_region;
			    $cdsregion_start = $cdsstart;
			}
			
			if ($starts[$i] <= $cdsstop && $cdsstop <= $stops[$i]) {
			    #cdsstop is in exon                                                                                               
			    my $five_region = ($cdsstop+1)."-".$stops[$i];
			    push @tmp_features,$enst."|5utr|".$five_region;
			    $cdsregion_stop = $cdsstop;
			}
			
			my $cds_region = ($cdsregion_start+1)."-".$cdsregion_stop;
			push @tmp_features,$enst."|CDS|".$cds_region;
		    }
		}
	    }
	    for (my $i=0;$i<scalar(@starts)-1;$i++) {
		push @tmp_features,$enst."|intron|".($stops[$i]+1)."-".$starts[$i+1];
	    }
	} else {
	    
	    for (my $i=0;$i<@starts;$i++) {
		push @tmp_features,$enst."|noncoding_exon|".($starts[$i]+1)."-".$stops[$i];
	    }
	    for (my $i=0;$i<scalar(@starts)-1;$i++) {
		push @tmp_features,$enst."|noncoding_intron|".($stops[$i]+1)."-".$starts[$i+1];
	    }
	}

	
	for my $feature (@tmp_features) {
	    my ($enst,$type,$region) = split(/\|/,$feature);
	    my ($reg_start,$reg_stop) = split(/\-/,$region);
	    my $x = int($reg_start/$hashing_value);
	    my $y = int($reg_stop /$hashing_value);

	    if ($reg_start > $reg_stop) {
		next;
	    }
	    
	    for my $j ($x..$y) {
		push @{$all_features{$chr}{$str}{$j}},$feature;
	    }
	}
    }
    close(F);

}




sub read_gencode_gtf {

#my %gene_types;
    my $file = shift;
#    my $file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
    print STDERR "Reading in $file\n";
    open(F,$file);
    for my $line (<F>) {
	chomp($line);
	next if ($line =~ /^\#/);
	my @tmp = split(/\t/,$line);
	
	my $stuff = $tmp[8];
	my @stufff = split(/\;/,$stuff);
	my ($ensg_id,$gene_type,$gene_name,$enst_id,$transcript_type);
	
	for my $s (@stufff) {
	    $s =~ s/^\s//g;
	    $s =~ s/\s$//g;
	    
	    if ($s =~ /gene_id \"(.+?)\"/) {
		if ($ensg_id) {
		    print STDERR "two ensg ids? $line\n";
		}
		$ensg_id = $1;
	    }
	    if ($s =~ /transcript_id \"(.+?)\"/) {
		if ($enst_id) {
		    prinst DERR "two enst ids? $line\n";
		}
		$enst_id = $1;
	    }
	    if ($s =~ /gene_type \"(.+?)\"/) {
		if ($gene_type) {
		    print STDERR "two gene types $line\n";
		}
		$gene_type = $1;
		
	    }
	    if ($s =~ /transcript_type \"(.+?)\"/) {
		$transcript_type = $1;
	    }
	    if ($s =~ /gene_name \"(.+?)\"/) {
		$gene_name = $1;
	    }
	}

	if (exists $enst2ensg{$enst_id} && $ensg_id ne $enst2ensg{$enst_id}) {
	    print STDERR "error two ensgs for enst $enst_id $ensg_id $enst2ensg{$enst_id}\n";
	}
	$enst2ensg{$enst_id} = $ensg_id;
	$ensg2name{$ensg_id}{$gene_name}=1;
	$ensg2type{$ensg_id}{$gene_type}=1;
	$enst2type{$enst_id} = $transcript_type;
    }
    close(F);
    
}
