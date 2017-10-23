use warnings;
use strict;
use Statistics::Basic qw(:all);

# 0 = do all pairwise comparisons
# 1 only do peak comparisons vs its own peak calls, not vs peak calls of other datasets
# 2 only do peak comparisons within a uID across replicates                                                                          

my $within_dataset_only_flag = 1;

my %rbp_list;
my @peak_filelist;
my %clipshort_to_peakshort;
my @clip_bam_fi_shorts;
my @peak_filelist_shorts;

my $filist_fi = $ARGV[0];
my $output_folder = $ARGV[1];
my $species = $ARGV[2];

unless ($filist_fi && $output_folder) {
    print STDERR "Usage Error\n";
    print STDERR "Usage: perl compare_vs_input_peakbased_wrapper2_TwoRepVersion_PEbamfileversion_submit.pl filist_fi output_folder\n";
    exit;
}

my %mapped_num;
my $filist_fi_mappedreadnums = $filist_fi.".mapped_read_num";
if (-e $filist_fi_mappedreadnums) {
    print STDERR "reading $filist_fi_mappedreadnums\n";
    open(READNUM,"$filist_fi_mappedreadnums");
    for my $line (<READNUM>) {
	chomp($line);
	my ($fi,$num) = split(/\t/,$line);
	next unless ($num);
	$mapped_num{$fi} = $num;
    }
    close(READNUM);
}

unless ($output_folder =~ /\/$/) {
    $output_folder = $output_folder."/";
}

my @rep_listing;
my $type_flag;

my %peakfi2uid;

#my $filist_fi = "Stefan_filelist_033015";
open(F,$filist_fi);
for my $line (<F>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\s+/,$line);
    next unless ($tmp[0]);
    next if ($tmp[0] eq "uID");
    if (scalar(@tmp) == 6) {
	$type_flag = "two_replicate_ENCODEstyle";
    } elsif (scalar(@tmp) == 5) {
	$type_flag = "one_replicate";
    }

    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);
    my %CLIP;
    $CLIP{"_01"} = shift(@tmp);
    $CLIP{"_02"} = shift(@tmp) if ($type_flag eq "two_replicate_ENCODEstyle");
    my $input = shift(@tmp);

    $CLIP{"_01"} =~ s/\.bam$//;
    $CLIP{"_02"} =~ s/\.bam$// if ($type_flag eq "two_replicate_ENCODEstyle");
    $input =~ s/\.bam$//;

#    print "uid $uid, rbp $rbp, cell $cellline, clip $CLIP{_01}, input $input\n";
#    exit;


    if ($type_flag eq "two_replicate_ENCODEstyle") {
	@rep_listing = ("_01","_02");
    } elsif ($type_flag eq "one_replicate") {
	@rep_listing = ("_01");
    } else {
	print STDERR "TYPE flag is not set properly!!!!\n";
    }

    next if ($CLIP{"_01"} eq "NA" || $input eq "NA");
    if ($type_flag eq "two_replicate_ENCODEstyle") {
	next if ($CLIP{"_02"} eq "NA");
    }
#    print "line $line X\n";
    for my $rep (@rep_listing) {

	my @clip_fi_split = split(/\//,$CLIP{$rep});
	my $clip_fi_short = $clip_fi_split[$#clip_fi_split];

	my $clip_bam_fi_short = $clip_fi_short.".bam";
	my $clip_bam_fi_softlink = $output_folder.$clip_bam_fi_short;
	my $clip_bam_fi = $CLIP{$rep}.".bam";
	push @clip_bam_fi_shorts,$clip_bam_fi_short;
	
	my @input_fi_split = split(/\//,$input);
	my $input_fi_short = $input_fi_split[$#input_fi_split];

	my $input_bam_fi_short = $input_fi_short.".bam";
	my $input_bam_fi_softlink = $output_folder.$input_bam_fi_short;
	my $input_bam_fi = $input.".bam";
	
	my $peak_fi_short = $clip_fi_short.".peaks.bed";
	my $peak_fi_softlink = $output_folder.$peak_fi_short;
	my $peak_fi = $CLIP{$rep}.".peaks.bed";
#	my $peak_fi = "/projects/ps-yeolab3/encode/analysis/encode_v13_recalled_peaks/".$peak_fi_short;

	$peakfi2uid{$peak_fi} = $uid.$rep;
	push @peak_filelist,$peak_fi;
	push @peak_filelist_shorts,$peak_fi_short;
	$clipshort_to_peakshort{$clip_bam_fi_short} = $peak_fi_short;
	
#	my $output = $output_folder.$clip_bam_fi_short."_vs_".$input_bam_fi_short."_basedon_".$peak_fi_short;
#	my $output = $output_folder.$uid.$rep.".peaks.l2inputnormnew.bed";
	my $output = $output_folder.$uid.$rep.".basedon_".$uid.$rep.".peaks.l2inputnormnew.bed";
	
#	print STDERR "clip $clip_bam_fi\n";
#	print STDERR "input $input_bam_fi\n";
#	print STDERR "peak $peak_fi\n";
	unless (-e $clip_bam_fi && -e $input_bam_fi) {
	    print STDERR "ERROR ERROR one of these files doesn't exists $clip_bam_fi $input_bam_fi\n";
	    next;
	}
	system("ln -s $clip_bam_fi $clip_bam_fi_softlink") unless (-e $clip_bam_fi_softlink);
	system("ln -s $input_bam_fi $input_bam_fi_softlink") unless (-e $input_bam_fi_softlink);
	system("ln -s $peak_fi $peak_fi_softlink") unless (-e $peak_fi_softlink);

	unless (exists $mapped_num{$clip_bam_fi_short}) {
	    my $mapped_read_num = `samtools view -c -F 4 $clip_bam_fi_softlink`;
	    $mapped_num{$clip_bam_fi_short} = $mapped_read_num;
	}
	unless (exists $mapped_num{$input_bam_fi_short}) {
	    my $input_read_num = `samtools view -c -F 4 $input_bam_fi_softlink`;
	    $mapped_num{$input_bam_fi_short} = $input_read_num;
	}

	if (-e $output) {
	    print STDERR "skipping $output alreadydone\n";
	    next;
	}
    }
}
close(F);


open(READNUMOUT,">$filist_fi_mappedreadnums");
for my $fi (keys %mapped_num) {
    chomp($mapped_num{$fi});
    print READNUMOUT "$fi\t$mapped_num{$fi}\n";
}
close(READNUMOUT);


my $l2foldenriched_filelist = $filist_fi.".l2foldenr_files.txt";
open(LOUT,">$l2foldenriched_filelist");
my %foldenriched;
open(F,$filist_fi);
for my $line (<F>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\s+/,$line);
    next unless ($tmp[0]);
    next if ($tmp[0] eq "uID");



    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);
    my %CLIP;
    $CLIP{"_01"} = shift(@tmp);
    $CLIP{"_02"} = shift(@tmp) if ($type_flag eq "two_replicate_ENCODEstyle");
    my $input = shift(@tmp);

    $CLIP{"_01"} =~ s/\.bam$//;
    $CLIP{"_02"} =~ s/\.bam$// if ($type_flag eq "two_replicate_ENCODEstyle");
    $input =~ s/\.bam$//;

    
    my $next_flag = 0;
    for my $rep (@rep_listing) {
	$next_flag = 1 if ($CLIP{$rep} eq "NA");
    } 
    $next_flag = 1 if ($input eq "NA");
    next if ($next_flag == 1);
    
    print LOUT "$line";

#    print "line $line X\n";                                                                                                             
    for my $rep (@rep_listing) {

        my @clip_fi_split = split(/\//,$CLIP{$rep});
        my $clip_fi_short = $clip_fi_split[$#clip_fi_split];

        my $clip_bam_fi_short = $clip_fi_short.".bam";
        my $clip_bam_fi_softlink = $output_folder.$clip_bam_fi_short;
        my $clip_bam_fi = $CLIP{$rep}.".bam";

	my @input_fi_split = split(/\//,$input);
        my $input_fi_short = $input_fi_split[$#input_fi_split];

        my $input_bam_fi_short = $input_fi_short.".bam";
        my $input_bam_fi_softlink = $output_folder.$input_bam_fi_short;
	my $input_bam_fi = $input.".bam";

        my $peak_fi_short = $clip_fi_short.".peaks.bed";
        my $peak_fi_softlink = $output_folder.$peak_fi_short;
	my $peak_fi = $CLIP{$rep}.".peaks.bed";
#        my $peak_fi = "/projects/ps-yeolab3/encode/analysis/encode_v13_recalled_peaks/".$peak_fi_short;

	unless (-e $peak_fi && -e $input_bam_fi && -e $clip_bam_fi) {
	    print STDERR "ERROR ERROR ERROR one of these doesn't exist $peak_fi $clip_bam_fi $input_bam_fi\n";
	    next;
	}

	my @desired_peak_filelist = @peak_filelist;
	if ($within_dataset_only_flag == 1) {
	    @desired_peak_filelist = ();
	    $desired_peak_filelist[0] = $peak_fi;
	} elsif ($within_dataset_only_flag == 2) {
	    @desired_peak_filelist = ();
            $desired_peak_filelist[0] = $peak_fi;
	    $desired_peak_filelist[1] = $CLIP{&opposite_rep($rep)}.".peaks.bed" if ($type_flag eq "two_replicate_ENCODEstyle");
	}

	my $sh_out = "sort_".$clip_bam_fi_short.".sh";
	open(SH,">$sh_out");
	print SH "\#\!\/bin\/sh\n";
	print SH "#PBS -N ".$sh_out."\n";
	print SH "#PBS -o ".$sh_out.".out\n";
	print SH "#PBS -e ".$sh_out.".err\n";
	print SH "#PBS -V\n";
	print SH "#PBS -l walltime=4:00:00\n";
	print SH "#PBS -l nodes=1:ppn=1\n";
	print SH "#PBS -A yeo-group\n";
	print SH "#PBS -q condo\n";
#	print SH "#PBS -q home-scrm\n";
#	print SH "#PBS -q home-yeo\n";
	my $userid = `whoami`;
	chomp($userid);
	print SH "#PBS -M ".$userid."\@ucsd.edu\n";
	print SH "#PBS -m a\n";
	print SH "cd $output_folder\n";
	my $needs_to_be_submitted = 0;

	for my $compare_peakfi (@desired_peak_filelist) {
	    
	    my @compare_peakfii = split(/\//,$compare_peakfi);
	    my $compare_peakfi_short = $compare_peakfii[$#compare_peakfii];
	    my $compare_peakfi_softlink = $output_folder.$compare_peakfi_short;
#	    my $output = $output_folder.$clip_bam_fi_short."_vs_".$input_bam_fi_short."_basedon_".$compare_peakfi_short."_l2inputnormnew";
	    my $output = $output_folder.$uid.$rep.".basedon_".$peakfi2uid{$compare_peakfi}.".peaks.l2inputnormnew.bed";
	    print LOUT "\t$output";


	    my $output_compressed = $output.".compressed.bed";
	    
	    if (-e $output) {
		print STDERR "skipping $output alreadydone\n";
		if (-e $output_compressed) {
		} else {
		    system("perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl $output");
#		    $needs_to_be_submitted = 1;

		}
	    } else {
#		system("perl compare_vs_input_peakbased_l2inputandrpm.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output $output2 $outputchi");
#		print STDERR "command perl compare_vs_input_peakbased_newfoldenr.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output $output2 $outputchi\n";
#		print STDERR "command perl overlap_peakfi_with_bam_PE.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output\n";
		$needs_to_be_submitted = 1;
		
		print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/overlap_peakfi_with_bam_PE.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output\n";
#		system("perl overlap_peakfi_with_bam_PE.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output");
		print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/compress_l2foldenrpeakfi_for_replicate_overlapping_bedformat.pl $output\n";
	    }

	    my $annotated_output = $output_compressed.".annotated";
	    if (-e $annotated_output) {
		print STDERR "skipping $annotated_output\n";
	    } else {
		print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/annotate_peaks_bedformat.pl $output_compressed $species\n";
		$needs_to_be_submitted = 1;

	    }
	}
#	exit;
	close(SH);
	if ($needs_to_be_submitted == 1) {
	    system("qsub $sh_out");
	}
    }
    print LOUT "\n";
#    exit;
}
close(F);

exit if ($within_dataset_only_flag == 1 || $within_dataset_only_flag == 2);


for my $clip_bam1 (@clip_bam_fi_shorts) {
    print "$clip_bam1";
    for my $clip_bam2 (@clip_bam_fi_shorts) {
	
	my $peak_fi_1 = $clipshort_to_peakshort{$clip_bam1};
	my $output1 = $clip_bam1."_basedon_".$peak_fi_1;
	my $output2 = $clip_bam2."_basedon_".$peak_fi_1;

	my $correlation = correlation(\@{$foldenriched{$output1}},\@{$foldenriched{$output2}});
#	print "$clip_bam1\t$clip_bam2\t$correlation\n";
	print "\t$correlation";
    }
    print "\n";
}



sub opposite_rep {
    my $x = shift;
    if ($x eq "_01") {
	return("_02");
    } elsif ($x eq "_02") {
	return("_01");
    } else {
	print STDERR "opposite rep subroutine error $x\n";
    }
}
