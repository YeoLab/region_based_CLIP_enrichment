use warnings;
use strict;

my $fi1 = $ARGV[0];
my $species = $ARGV[1];

my @file_split = split(/\//,$fi1);
my @file_split2 = split(/\_/,$file_split[$#file_split]);
my $UID = $file_split2[0]."_".$file_split2[1];


my %repinfo_hash;
my $l2cutoff = 2;
my @rep_listing;
my $num_replicates;
my %labels;
my %data;
&parse_l2fcwithpval_fi($fi1);


if (scalar(keys %repinfo_hash) == 2) {
    @rep_listing = ($UID."_01",$UID."_02","01_02_overlap");
} elsif (scalar(keys %repinfo_hash) == 1) {
    @rep_listing = ($UID."_01");
} else {
    print STDERR "Rep number doesn't make sense ".join("|",keys %repinfo_hash)."\n";
    exit;
}

#$data{$ensg}{$repinfo}{$type} =$tmp[$i];

my $n_file = $fi1.".cutoff_counts.".$l2cutoff;
open(N,">$n_file");

for my $rep (@rep_listing) {
#    my $n_file = $fi1."_and_".$fi2short.".cutoff_counts.".sprintf("%.4f",$l2cutoff);
 
    for my $desired_label ("CDS","3utr","intron","noncoding_exon","5utr") {
#    for my $desired_i (0,2,4,6) {
	my $file_label = $fi1.".".$rep."_l2cutoff_".$l2cutoff.".".$desired_label;
	my $bgdfile_label = $fi1.".".$rep."_l2cutoff_".$l2cutoff.".".$desired_label.".bgd";
	my $all_count=0;
	my %cutoff_count;

    
	open(F,">$file_label");
	open(BGD,">$bgdfile_label");
	print "rep $rep\n";

	if ($rep eq "01_02_overlap") {
	    
	    for my $ensg (keys %data) {
		
		next unless (exists $data{$ensg}{$rep_listing[0]}{$desired_label} && exists $data{$ensg}{$rep_listing[1]}{$desired_label});
		my ($l2fc1,$l10pv1) = split(/\|/,$data{$ensg}{$rep_listing[0]}{$desired_label});
		my ($l2fc2,$l10pv2) = split(/\|/,$data{$ensg}{$rep_listing[1]}{$desired_label});

		if ($l2fc1 eq "NaN" || $l2fc2 eq "NaN") {
		    # not quantified in one of the datasets - skip
		} else {
		    print BGD "$ensg\n";
		    $all_count++;
		    
		    my ($fc1_flag,$fc2_flag) = (0,0);
		    my ($pv1_flag,$pv2_flag) = (0,0);
		    if ($l10pv1 >= 5) {
			$pv1_flag = 1;
		    }
		    if ($l10pv2 >= 5) {
			$pv2_flag = 1;
		    }
		    if ($l2fc1 >= $l2cutoff) {
			$fc1_flag=1;
		    } 
		    if ($l2fc2 >= $l2cutoff) {
			$fc2_flag = 1;
		    }
			
		    if ($fc1_flag == 1 && $fc2_flag == 1 && $pv1_flag == 1 && $pv2_flag == 1) {
			print F "$ensg\n";
		    }
		    $cutoff_count{$pv1_flag."|".$fc1_flag."_".$pv2_flag."|".$fc2_flag}++;
		}
	    }
	    print N "$fi1\t$rep\t$desired_label\t$all_count";
	    for my $pv1_flag (0,1) {
		for my $pv2_flag (0,1) {
		    for my $fc1_flag (0,1) {
			for my $fc2_flag (0,1) {
			    my $counted_type = $pv1_flag."|".$fc1_flag."_".$pv2_flag."|".$fc2_flag;
			    $cutoff_count{$counted_type} = 0 unless (exists $cutoff_count{$counted_type});
			    print N "\t$counted_type\t$cutoff_count{$counted_type}";
			}
		    }
		}
	    }
	    print N "\n";
	} elsif ($rep =~ /\_01$/ || $rep =~ /\_02$/) {
	    for my $ensg (keys %data) {
		next unless (exists $data{$ensg}{$rep}{$desired_label});
		my ($l2fc1,$l10pv1) = split(/\|/,$data{$ensg}{$rep}{$desired_label});
		if ($l2fc1 eq "NaN") {
		    next;
		}
		$all_count++;

		my $fc1_flag = 0;
		my $pv1_flag = 0;

		print BGD "$ensg\n";
		if ($l10pv1 >= 5) {
		    $pv1_flag = 1;
		    if ($l2fc1 >= $l2cutoff) {
			$fc1_flag=1;
			print F "$ensg\n";
		    }
		}
		$cutoff_count{$pv1_flag."|".$fc1_flag}++;
		
	    }
	    print N "$fi1\t$rep\t$desired_label\t$all_count";
	    for my $types ("1|1","1|0","0|1","0|0") {
		$cutoff_count{$types} = 0 unless (exists $cutoff_count{$types});
		print N "\t$types\t$cutoff_count{$types}";
	    }
	    print N "\n";
	}
	
	close(F);
	close(BGD);
    
	my $go_outputfi = $file_label.".go_output";
#	system("perl ensembl_GO_enrichmentpvalue.pl $bgdfile_label $file_label $go_outputfi $species");
	system("perl GOparsed_GO_enrichmentpvalue.pl $bgdfile_label $file_label $go_outputfi $species");
    }
}
close(N);






sub parse_l2fcwithpval_fi {
    my $fi = shift;
    open(F,$fi);
    for my $line (<F>) {
	chomp($line);
	my @tmp = split(/\t/,$line);
	my $ensg = shift(@tmp);
	if ($ensg eq "ENSG") {
	    for (my $i=0;$i<@tmp;$i++) {
		my @label_i = split(/\|/,$tmp[$i]);
		my $repinfo = shift(@label_i);
		my $type = join("|",@label_i);
#		print "i $i rep $repinfo type $type\n";
		$labels{$i}{repinfo} = $repinfo;
		$labels{$i}{type} = $type;
	    }
	} else {
	    for (my $i=0;$i<@tmp;$i++) {
		$data{$ensg}{$labels{$i}{repinfo}}{$labels{$i}{type}} = $tmp[$i];
#		print "$i label $labels{$i}{repinfo} type $labels{$i}{type} $tmp[$i]\n";
		$repinfo_hash{$labels{$i}{repinfo}} = 1;
	    }
	}
    }
    close(F);
}
