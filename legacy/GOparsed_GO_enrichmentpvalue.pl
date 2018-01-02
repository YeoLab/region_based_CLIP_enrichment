use warnings;
use strict;
use Statistics::Basic qw(:all);
use Statistics::Distributions;

my %precalculated_fisher;

my $species = $ARGV[3];
my $go_file = "/home/elvannostrand/data/gene_ontology/GO_parsedannotations_inclchildandparent.20160325.txt";
if ($species eq "hg19") {
} elsif ($species eq "mm9") {
    print STDERR "mouse annotations don't exist yet - error\n";
    exit;
#    $go_file = "/home/elvannostrand/data/gene_ontology/GO_annotations_fromEnsembl_mm9_111215.txt.gz";
} else {
    print STDERR "species mismatch $species\n";
}

my %go_hash;
open(GO,$go_file);
#open(GO,"gunzip -c $go_file |");
for my $line (<GO>) {
    chomp($line);
    next if ($line =~ /^Ensembl/);

#    my ($ensg,$enst,$goid,$goname,$godefn,$godomain,$gotermevidence) = split(/\t/,$line);
    my ($ensg,$gname,$goid,$go_flag,$goname) = split(/\t/,$line);
    next unless ($goid);
    next unless ($ensg);
    $go_hash{$goid}{name} = $goname;
    $go_hash{$goid}{ensgs}{$ensg} = 1;

}
close(GO);

my %bgd_hash;
my %desired_hash;
my $bgd_fi = $ARGV[0];
open(BGD,$bgd_fi);
for my $line (<BGD>) {
    chomp($line);
    my $ensg = $line;
    if ($line =~ m/^(ENSG\d+)\./) {
	$ensg = $1;
    }
    $bgd_hash{$ensg} = 1;
}
close(BGD);

my $n_bgd_genes = scalar(keys %bgd_hash);

my $input_fi = $ARGV[1];

my $fisher_tmp_fi = $input_fi.".tmp_fisher";
my $fisher_tmp_fi_out = $input_fi.".tmp_fisher.out";
my $fisher_tmp_fi_s = $fisher_tmp_fi.".s";


open(I,$input_fi);
for my $line (<I>) {
    chomp($line);
    my $ensg = $line;
    if ($line =~ m/^(ENSG\d+)\./) {
	$ensg =$1;
    }
    
    $desired_hash{$ensg} = 1;
    unless (exists $bgd_hash{$ensg}) {
	print STDERR "weird - differential expressed gene isn't in background list $line $ensg\n";
    }
}
close(I);
my $n_input_genes = scalar(keys %desired_hash);

#my $num_of_go_terms = scalar(keys %go_hash);
my $num_of_go_terms = 0;

for my $goid (keys %go_hash) {
    my $goid_count = 0;
    for my $ensg (keys %{$go_hash{$goid}{ensgs}}) {
        next unless (exists $bgd_hash{$ensg});
        $goid_count++;
    }
    if ($goid_count >= 5) {
	$num_of_go_terms++;
    }
}

my $out_file = $input_fi.".go_output.csv";
open(OUT,">$out_file");

my %output;
my @goidlist;
for my $goid (keys %go_hash) {


    my %count;    
    my $input_match = 0;
    my $goid_count = 0;
    for my $ensg (keys %{$go_hash{$goid}{ensgs}}) {
	next unless (exists $bgd_hash{$ensg});
	$goid_count++;
	$input_match++ if (exists $desired_hash{$ensg});
    }
    next unless ($goid_count >= 5);

    my $chisq_11 = $input_match;
    my $chisq_10 = $goid_count-$input_match;
    my $chisq_01 = $n_input_genes-$input_match;
    my $chisq_00 = $n_bgd_genes-$n_input_genes-$goid_count+$input_match;
 
    my ($chipval,$chival,$chitype,$chienrdepl) = &fisher_or_chisq($chisq_11,$chisq_10,$chisq_01,$chisq_00);
    my $log10pval = $chipval > 0 ? -1 * log($chipval)/log(10) : 400 ;
    
    my $log10pval_bfcorrect = $chipval > 0 ? -1 * log($chipval * $num_of_go_terms)/log(10) : 400 ;

    my $out_line = "$goid\t$chisq_11\t$chisq_10\t$chisq_01\t$chisq_00\t$chipval\t$chival\t$chitype\t$chienrdepl\t$log10pval\t$log10pval_bfcorrect\t".$go_hash{$goid}{name};
#    print "$goid\t$chisq\t$pval\n";
    $output{$out_line} = $log10pval_bfcorrect;
}
my @sorted = sort {$output{$b} <=> $output{$a}}  keys %output;
for my $s (@sorted) {
    print OUT "$s\n";
}
close(OUT);

system("rm $fisher_tmp_fi_out");
system("rm $fisher_tmp_fi");
system("rm $fisher_tmp_fi_s");





sub abs {
    my $x = shift;
    if ($x > 0) {
	return($x);
    } else {
	return(-1*$x);
    }
}


sub square {
    my $x = shift;
    return($x * $x);
}


sub fisher_or_chisq {
    my ($a,$b,$c,$d) = @_;
    unless ($a && $b && $c && $d) {
#        return("1","NA","NA");
    }

    my $tot = $a + $b + $c + $d;
    my $expa = ($a+$c)*($a+$b)/$tot;
    my $expb = ($b+$d)*($a+$b)/$tot;
    my $expc = ($a+$c)*($c+$d)/$tot;
    my $expd = ($b+$d)*($c+$d)/$tot;

    my $direction = "enriched";
    if ($a<$expa) {
        $direction = "depleted";
        return(1,"DEPL","N",$direction);
    }



    if ($expa < 5 || $expb < 5 || $expc < 5 || $expd < 5 || $a < 5 || $b < 5 || $c < 5 || $d < 5) {
        if (exists $precalculated_fisher{$a."|".$b}) {
            return($precalculated_fisher{$a."|".$b}{p},$precalculated_fisher{$a."|".$b}{v},"F",$direction);
        } else {
            my ($pval,$val) = &fisher_exact($a,$b,$c,$d);
            $precalculated_fisher{$a."|".$b}{p} = $pval;
            $precalculated_fisher{$a."|".$b}{v} = $val;
            return($pval,$val,"F",$direction);
        }
    } else {
        my ($pval,$val) = &chi_square($a,$b,$c,$d);
        return($pval,$val,"C",$direction);
    }
}


sub chi_square {
    my ($a,$b,$c,$d) = @_;
    #    print "$a\t$b\t$c\t$d\t";
    return(0) unless ($a && $b && $c && $d);
    #    $b = $b-$a;
    #    $c = $c-$a;
    #    $d = $d-$c-$b-$a;
    #    $d = $d - $c;

#        print "$a\t$b\t$c\t$d\t";
#    if ($a >= 5 && $b >= 5 && $c >= 5 && $d >= 5 ){
    my $tot = $a + $b + $c + $d;
    my $expa = ($a+$c)*($a+$b)/$tot;
    my $expb = ($b+$d)*($a+$b)/$tot;
    my $expc = ($a+$c)*($c+$d)/$tot;
    my $expd = ($b+$d)*($c+$d)/$tot;
    
    if ($expa >= 5 || $expb >= 5 || $expc >= 5 || $expd >= 5) {
        my $chival = &square(&abs($a-$expa)-0.5)/$expa +  &square(&abs($b-$expb)-0.5)/$expb + &square(&abs($c-$expc)-0.5)/$expc + &square(&abs($d-$expd)-0.5)/$expd;
        
        my $pval = Statistics::Distributions::chisqrprob(1,&abs($chival));

        if ($a<$expa) {
            $chival = $chival * -1;
        }
        return ($pval,$chival);
    } else {
         #       print "\n";
        print STDERR "shouldn't get to this - should have been shunted into fisher exact test\n";
        return(1);
    }
}

sub fisher_exact {
    my ($x1,$x2,$y1,$y2) = @_;

    system("rm $fisher_tmp_fi_out");
    system("rm $fisher_tmp_fi");
    system("rm $fisher_tmp_fi_s");

    open(LISTA,">$fisher_tmp_fi");
    print LISTA "".$x1."\t".$x2."\t".$y1."\t".$y2."\n";
    close(LISTA);

    my $p_value_vs_bgd = "";

    while (length($p_value_vs_bgd) < 1) {
        open(RFI,">$fisher_tmp_fi_s");
        printf RFI "A <- scan\(\"".$fisher_tmp_fi."\"\)\n";
        printf RFI "blah <- matrix(c(A),nrow=2)\n";
        printf RFI "sink(\"".$fisher_tmp_fi_out."\")\n";
        printf RFI "fisher.test(blah)\n";
        printf RFI "sink()\n";
        system("  R --save --quiet < $fisher_tmp_fi_s");


        open(FISH,$fisher_tmp_fi_out);
        for my $line (<FISH>) {
            chomp($line);
            if ($line =~ /p\-value\s\=\s(\S+)$/ || $line =~ /p\-value\s\<\s(\S+)$/) {
                $p_value_vs_bgd = $1;
            }
	}
        close(FISH);
    }
    my $val = "F";
    return($p_value_vs_bgd,$val);

}





