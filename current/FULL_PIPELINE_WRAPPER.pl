use warnings;
use strict;

my $filelist = $ARGV[0];
my $output_folder = $ARGV[1];
my $species = "hg19";
if (exists $ARGV[2]) {
    $species = $ARGV[2];
}


unless ($filelist && $output_folder) {
    print "Usage: perl FULL_PIPELINE_WRAPPER.pl filelist output_folder species\n";
    exit;
}

unless ($output_folder =~ /\/$/) {
    $output_folder = $output_folder."/";
}


system("perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/compare_vs_input_peakbased_wrapper2_OneOrTwoRepVersion_PEbamfileversion_submit.pl $filelist $output_folder $species");
system("perl /home/elvannostrand/data/clip/CLIPseq_analysis/scripts/LABshare/count_reads_OneOrTwoRepVersion.wrapper_organized_frombamfi_PE.pl $filelist $output_folder $species");


	    

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

