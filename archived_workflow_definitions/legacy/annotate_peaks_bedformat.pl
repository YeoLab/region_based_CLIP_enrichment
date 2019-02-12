use warnings;
use strict;

# 2015/11/12 - fixed to deal with ucsc table format (0-based, open ended) vs peaks (0-based, open ended)
# & add genome (mm9) option

my $verbose_flag = 0;
my $hashing_value = 100000;


my %all_features;
my %enst2type;
my %enst2ensg;
my %ensg2name;
my %ensg2type;


#defaults to hg19
my $species = "hg19";
if (exists $ARGV[1]) {
    $species = $ARGV[1];
}


my $gencode_gtf_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf";
my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
#my $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/hg19/gencode_v19/gencodev19_comprehensive"; 

if ($species eq "mm9") {
    $gencode_gtf_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/mm9/gencode.vM1.annotation.gtf.parsed_ucsc_tableformat";
} elsif ($species eq "hg38") {
    $gencode_gtf_file = "/projects/ps-yeolab/genomes/GRCh38/gencode/v24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf";
    $gencode_tablebrowser_file = "/projects/ps-yeolab/genomes/GRCh38/gencode/v24/gencode.v24.chr_patch_hapl_scaff.annotation.gtf.parsed_ucsc_tableformat";
}

print STDERR "using $species $gencode_tablebrowser_file\n";
&read_gencode_gtf($gencode_gtf_file);
&read_gencode($gencode_tablebrowser_file);

my $peak_fi = $ARGV[0];
my $output = $peak_fi.".annotated";
open(OUT,">$output");
&read_peak_fi($peak_fi);
close(OUT);



sub read_peak_fi {
    my $peakfi = shift;
    open(PEAK,$peakfi) || die "no peakfi $peakfi\n";
    for my $line (<PEAK>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $chr = $tmp[0];
	my $str = $tmp[5];
	my $start = $tmp[1];
	my $stop = $tmp[2];
	
	my %tmp_hash;
	my $feature_flag = 0;

	my $rx = int($start / $hashing_value);
	my $ry = int($stop  / $hashing_value);
	for my $ri ($rx..$ry) {
	    for my $feature (@{$all_features{$chr}{$str}{$ri}}) {
		my ($feature_enst,$feature_type,$feature_region) = split(/\|/,$feature);
		my ($feature_start,$feature_stop) = split(/\-/,$feature_region);
		if ($feature_start <= $start && $stop <= $feature_stop) {
                        # peak is entirely within region                                                                                                                                                                                                         
		    my $feature_ensg = $enst2ensg{$feature_enst};
		    $tmp_hash{$feature_type}{$feature_ensg}="contained";
		    $feature_flag = 1;
		    print "found feature $_ $feature $feature_ensg $feature_type\n" if ($verbose_flag == 1);
		} elsif (($feature_start >= $start && $feature_start < $stop) || ($feature_stop > $start && $feature_stop <= $stop)) {
		    # peak is partly overlapping feature
                    my $feature_ensg = $enst2ensg{$feature_enst};
		    $tmp_hash{$feature_type}{$feature_ensg}="partial";
		    $feature_flag = 1;
		} elsif ($start < $feature_start && $stop > $feature_stop) {
		    # feature is contained within peak
		    my $feature_ensg = $enst2ensg{$feature_enst};
		    $tmp_hash{$feature_type}{$feature_ensg}="featurewithin";
		    $feature_flag = 1;
		}

	    }
	}

	my $final_feature_type = "intergenic";
	my $final_feature_ensg = "NA";
	if ($feature_flag == 1) {
	    if (exists $tmp_hash{"CDS"}) {
		my $feature_type_flag = &get_type_flag(\%tmp_hash,"CDS");
		$final_feature_type = "CDS|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"CDS"}});
		$final_feature_ensg = join("|",keys %{$tmp_hash{"CDS"}});
	    } elsif (exists $tmp_hash{"3utr"} || exists $tmp_hash{"5utr"}) {
		if (exists $tmp_hash{"3utr"} && exists $tmp_hash{"5utr"}) {
		    my $feature_type_flag = &get_type_flag(\%tmp_hash,"3utr");
		    my $feature_type_flag2 = &get_type_flag(\%tmp_hash,"5utr");

		    unless ($feature_type_flag eq "contained" && $feature_type_flag2 eq "contained") {
			$feature_type_flag = "partial";
		    }

		    $final_feature_type = "5utr_and_3utr|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"5utr"}}).":".join(":",keys %{$tmp_hash{"3utr"}});
		    $final_feature_ensg = join("|",keys %{$tmp_hash{"3utr"}},keys %{$tmp_hash{"5utr"}});
		} elsif (exists $tmp_hash{"3utr"}) {
		    my $feature_type_flag = &get_type_flag(\%tmp_hash,"3utr");
		    $final_feature_type = "3utr|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"3utr"}});
		    $final_feature_ensg = join("|",keys %{$tmp_hash{"3utr"}});
		} elsif (exists $tmp_hash{"5utr"}) {
		    my $feature_type_flag = &get_type_flag(\%tmp_hash,"5utr");
		    $final_feature_type = "5utr|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"5utr"}});
		    $final_feature_ensg = join("|",keys %{$tmp_hash{"5utr"}});
		} else {
		    print STDERR "weird shouldn't hit this\n";
		}  
	    } elsif (exists $tmp_hash{"intron"}) {
		my $feature_type_flag = &get_type_flag(\%tmp_hash,"intron");
		$final_feature_type = "intron|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"intron"}});
		$final_feature_ensg = join("|",keys %{$tmp_hash{"intron"}});
	    } elsif (exists $tmp_hash{"noncoding_exon"}) {
		my $feature_type_flag = &get_type_flag(\%tmp_hash,"noncoding_exon");
		$final_feature_type = "noncoding_exon|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"noncoding_exon"}});
		$final_feature_ensg = join("|",keys %{$tmp_hash{"noncoding_exon"}});
	    } elsif (exists $tmp_hash{"noncoding_intron"}) {
		my $feature_type_flag = &get_type_flag(\%tmp_hash,"noncoding_intron");
		$final_feature_type = "noncoding_intron|".$feature_type_flag."|".join(":",keys %{$tmp_hash{"noncoding_intron"}});
		$final_feature_ensg = join("|",keys %{$tmp_hash{"noncoding_intron"}});
	    }
	}
	
	print OUT "$line\t$final_feature_type\t$final_feature_ensg\n";

    }
    close(PEAK);
}





sub get_type_flag {
    my $ref = shift;
    my %feature_hash = %$ref;
    my $feature_type = shift;

    my $feature_type_final = "contained";
    for my $ensg (keys %{$feature_hash{$feature_type}}) {
	$feature_type_final = "partial" unless ($feature_hash{$feature_type}{$ensg} eq "contained");
    }
    return($feature_type_final);
}


sub read_gencode {
    ## eric note: this has been tested for off-by-1 issues with ucsc brower table output!                                                                                                                                                              
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

            for my $j ($x..$y) {
                push @{$all_features{$chr}{$str}{$j}},$feature;
            }
        }
    }
    close(F);
    
}


sub read_gencode_gtf {

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
                    print STDERR "two enst ids? $line\n";
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



	$ensg2name{$ensg_id}{$gene_name}=1;
	$ensg2type{$ensg_id}{$gene_type}=1;

	if ($enst_id) {
	    if (exists $enst2ensg{$enst_id} && $ensg_id ne $enst2ensg{$enst_id}) {
		print STDERR "error two ensgs for enst $enst_id $ensg_id $enst2ensg{$enst_id}\n";
	    }
	    
	    $enst2ensg{$enst_id} = $ensg_id;
	    $enst2type{$enst_id} = $transcript_type;
	}
    }
    close(F);

}
