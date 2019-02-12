use warnings;
use strict;

my $filist_fi = $ARGV[0];
my $working_directory = $ARGV[1];
my $species = $ARGV[2];

unless ($working_directory =~ /\/$/) {
    $working_directory = $working_directory."/";
}


my @rep_listing;
my $type_flag;
my %alreadydone;
#my $working_directory = "/home/elvannostrand/data/clip/CLIPseq_analysis/IGF2BP_analysis/";
#my $filist_fi = $working_directory."filelist_042415";
#my $filist_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/Stefan_filelist_033015";
#my $filist_fi = "/home/elvannostrand/data/clip/CLIPseq_analysis/input_comparison/ENCODE_current_file_list.txt";
open(F,$filist_fi) || die "no $filist_fi\n";
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
    my %output_fis;
    $CLIP{"_01"} = shift(@tmp);
    $CLIP{"_02"} = shift(@tmp) if ($type_flag eq "two_replicate_ENCODEstyle");
#    my $input = shift(@tmp);
    $CLIP{"input"} = shift(@tmp);
    next if ($CLIP{"_01"} eq "NA");
    
    $CLIP{"_01"} =~ s/\.bam$//;
    $CLIP{"_02"} =~ s/\.bam$// if ($type_flag eq "two_replicate_ENCODEstyle");
    $CLIP{"input"} =~ s/\.bam$//;

    if ($type_flag eq "two_replicate_ENCODEstyle") {
        @rep_listing = ("_01","_02","input");
    } elsif ($type_flag eq "one_replicate") {
        @rep_listing = ("_01","input");
    } else {
        print STDERR "TYPE flag is not set properly!!!!\n";
    }

    my $clip_id = $uid."_".$rbp;
    my $sh_out = $clip_id.".reads_by_loc.sh";
    open(SH,">$sh_out");
    print SH "\#\!\/bin\/sh\n";
    
    print SH "#PBS -N reads_loc_".$clip_id."\n";
    print SH "#PBS -o reads_loc_".$clip_id.".out\n";
    print SH "#PBS -e reads_loc_".$clip_id.".err\n";
    print SH "#PBS -V\n";
    print SH "#PBS -l walltime=6:00:00\n";
    print SH "#PBS -l nodes=1:ppn=1\n";
    print SH "#PBS -A yeo-group\n";
    print SH "#PBS -q home-yeo\n";
#    print SH "#PBS -M elvannostrand\@ucsd.edu\n";
    my $userid = `whoami`;
    chomp($userid);
    print SH "#PBS -M ".$userid."\@ucsd.edu\n";

    print SH "#PBS -m a\n";
    print SH "cd /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/\n";
    
    my $missing_dataset_flag = 1;
    my $needs_to_submit_flag = 0;
#    for my $rep ("_01","_02","input") {
    for my $rep (@rep_listing) {
	next if ($CLIP{$rep} eq "NA");
	my @clip_fi_split = split(/\//,$CLIP{$rep});
	my $clip_fi_short = $clip_fi_split[$#clip_fi_split];

	my $clip_bam_fi_short = $clip_fi_short.".bam";
	my $clip_bam_fi = $CLIP{$rep}.".bam";
        my $clip_bam_fi_softlink = $working_directory.$clip_bam_fi_short;

	unless (-e $clip_bam_fi) {
	    $missing_dataset_flag .= $clip_bam_fi."_";
	    next;
	}
	
	system("ln -s $clip_bam_fi $clip_bam_fi_softlink") unless (-e $clip_bam_fi_softlink);
	
	my $output_fi = $working_directory.$clip_bam_fi_short.".reads_by_loc.csv";
	$output_fis{$rep} = $output_fi;

	if (exists $alreadydone{$output_fi}) {
	    print STDERR "This is a critical error - this dataset shares an output with a previous dataset $output_fi\n";
	}
	next if (exists $alreadydone{$output_fi});
	next if (-e $output_fi);
	print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/count_reads_broadfeatures_frombamfi_PEmap.pl $clip_bam_fi_softlink $output_fi $species\n";
	$needs_to_submit_flag = 1;
	$alreadydone{$output_fi} = 1;
    }
    
    if ($missing_dataset_flag ne "1") {
	print STDERR "ERROR ERROR dataset is missing $missing_dataset_flag\n";
	next;
    }

    my $overall_output = $working_directory.$uid."_".$rbp."_ReadsByLoc_combined.csv";
			
    if (-e $overall_output) { } else {
	print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/combine_ReadsByLoc_files.pl";
	for my $rep (@rep_listing) {
	    print SH " ".$output_fis{$rep};
	}
	print SH " > $overall_output\n";
	$needs_to_submit_flag = 1;
    }
    
    my $l2fc_output = $overall_output.".l2fcwithpval_enr.csv";
    if (-e $l2fc_output) { } else {
	my $mapped_read_fi = $filist_fi.".mapped_read_num";
	print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/convert_ReadsByLoc_combined_significancecalls.pl $overall_output $mapped_read_fi\n";
	$needs_to_submit_flag = 1;
    }

    my $go_output = $l2fc_output.".cutoff_counts.2";
    if (-e $go_output) { } else {
	print SH "perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/regionlevelanalysis_GOanalysisonl2fc.pl $l2fc_output $species\n";
	$needs_to_submit_flag = 1;
    }
    close(SH);
    
    system("qsub $sh_out") if ($needs_to_submit_flag == 1);
    
#    exit;
}
close(F);
