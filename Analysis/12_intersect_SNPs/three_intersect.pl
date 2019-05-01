#!/usr/bin/perl
##by Li Lei, 20190403, in St.Paul;
#this is to intersect the markers with the three datasets.
#usage: ./three_intersect.pl intersected_by_two_SNPsid third_sets_SNPid |grep "YES" >intersect_by_three_datasets
use strict;
use warnings;
use Data::Dumper;
my ($file1, $file2) = @ARGV;

my %gidhash;


open(SNPID,  "$file1") or die "Could not open $file1";

foreach my $row (<SNPID>){
        chomp $row;
        my @rtemp = split(/\t/,$row);
        my $g_id = $rtemp[0];
           $gidhash{$g_id}=$row;
        #print "$g_id\n";
}
close (SNPID);
#print Dumper(\%gidhash);
open(OUT,  "$file2") or die "Could not open $file2";

foreach my $row (<OUT>){
        chomp $row;
        my @rtemp1 = split(/\t/,$row);
           #print "$rtemp1[0]\n";
           my $SNP = $rtemp1[0]."_".$rtemp1[1];
           #print "$SNP\n";
        if(exists $gidhash{$SNP}){
           print "$SNP\tYES\n";
        }
        else{
           #
           print "$SNP\tNO\n";
        }
}
close (OUT);