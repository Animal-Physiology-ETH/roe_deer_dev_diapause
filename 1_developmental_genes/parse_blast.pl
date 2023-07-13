#!/usr/bin/perl

use strict;
use warnings;

my$blsFile1="transcripts_vs_dev-genes.bls.tabular";
my$blsFile2="dev-genes_vs_transcripts.bls.tabular";

my$outFile="annotation.tab";

my%results;
my%final;

open(IN1,'<',$blsFile1) or die $!;
while(<IN1>)
{
    chomp;
    my@line=split(/\s+/,$_);

    my$trID=$line[0];
    (my$gene)=$line[1]=~/^(.+)_.+_.+/;
    my$eval=$line[10];

    if($eval == 0)
    {
        $results{$trID}{$gene}{"high"}++;
    }
    elsif($eval <= 1e-3)
    {
        $results{$trID}{$gene}{"medium"}++;
    }
    else
    {
        $results{$trID}{$gene}{"low"}++;
    }

    $results{$trID}{"totalCount"}++;
}
close IN1 or die $!;

open(IN2,'<',$blsFile2) or die $!;
while(<IN2>)
{
    chomp;
    my@line=split(/\s+/,$_);

    (my$gene)=$line[0]=~/^(.+)_.+_.+/;
    my$trID=$line[1];
    my$eval=$line[10];

    if($eval <= 1e-3)
    {
        $results{$trID}{$gene}{"confirmed"}++;
    }
}
close IN2 or die $!;


foreach my$tr (keys %results)
{
    foreach my$gen (keys %{$results{$tr}})
    {

        next if $gen eq "totalCount";

        my$high = $results{$tr}{$gen}{"high"} || 0;
        my$medium = $results{$tr}{$gen}{"medium"} || 0;
        my$low = $results{$tr}{$gen}{"low"} || 0;
        my$total = $results{$tr}{"totalCount"}; 
        my$conf;

        ### TODO ###
        # I need a way to make sure, only one (the best) is printed
        # currently, it is the first
        # maybe use a final hash?


       if($high/$total >= 0.8 && exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=3;
        }
        elsif($high/$total >= 0.8 && !exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=2;
        }
        elsif(($high+$medium)/$total >= 0.8 && exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=3;
        }
        elsif(($high+$medium)/$total >= 0.8 && !exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=2;
        }
        elsif(($high+$medium)/$total >= 0.5 && exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=2;
        }
        elsif(($high+$medium+$low)/$total >= 0.8 && exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=2;
        }
        elsif(($high+$medium+$low)/$total >= 0.8 && !exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=1;
        }
        elsif(($high+$medium+$low)/$total >= 0.5 && exists $results{$tr}{$gen}{"confirmed"})
        {
            $conf=1;
        }
        else
        {
            $conf=0;
        }

        if(!exists $final{$tr} || $final{$tr}{"confidence"} < $conf)
        {
            $final{$tr}{"gene"}=$gen;
            $final{$tr}{"confidence"}=$conf;
        }
    }
}

open(OUT,'>',$outFile) or die $!;
foreach my$tr (keys %final)
{
    my$conf="NA";
    if($final{$tr}{"confidence"} == 3)
    {
        $conf="high"
    }
    elsif($final{$tr}{"confidence"} == 2)
    {
        $conf="medium"
    }
    elsif($final{$tr}{"confidence"} == 1)
    {
        $conf="low"
    }

    my$gene="undecided";
    $gene=$final{$tr}{"gene"} unless $conf eq "NA";

    print OUT join("\t",$tr,$gene,$conf),"\n";
}

close OUT or die $!;
