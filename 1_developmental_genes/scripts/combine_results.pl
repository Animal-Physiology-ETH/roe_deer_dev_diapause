#!/usr/bin/perl

use strict;
use warnings;

my$annFile="annotation.tab";
my@libNames=qw(E102 E103 E1041 E106 E1081 E1082 E1091 E1101 E1141 E1142 E1181 E1182 E1191 E1192 E1231 E128 E1281
E131 E1311 E1321 E1322 E1371 E1401 E1402 E1421 E1422 E1591 E1681 E1721 E1722 E181 E1821 E1822 E1831
E1832 E1861 E1871 E1872 E1901 E1902 E1911 E1912 E1921 E1922 E1931 E1932 E1941 E2011 E2012 E2031 E2032
E2091 E2092 E212 E2131 E2181 E2182 E2191 E2221 E2222 E2251 E231 E281 E291 E301 E321 E381 E382 E401
E472 E491 E501 E511 E591 E631 E632 E691 E721 E722 E741 E781 E782 E801 E811 E881 E98 E991);

my$outFile="combined.tab";

my%results;

open(ANN,'<',$annFile) or die $!;
while(<ANN>)
{
    chomp;
    my@line=split(/\t/,$_);

    $results{$line[0]}{"Ann"}=$line[1];
    $results{$line[0]}{"Conf"}=$line[2];

}
close ANN or die $!;

foreach my$lib (@libNames)
{
    open(IN,'<',"./Salmon/".$lib.".fastq.tabular") or die $!;
    <IN>;
    while(<IN>)
    {
        chomp;
        my@line=split(/\s+/,$_);
        
        if(!exists $results{$line[0]})   # there are some transcripts without hits in the blast search - they are not in the anno file, this here is to set a default for them
        {
            $results{$line[0]}{"Ann"}="unknown";
            $results{$line[0]}{"Conf"}="NA";
        }

        push(@{$results{$line[0]}{"TPM"}},$line[3]);
    }
    close IN or die $;
}


open(OUT,'>',$outFile) or die $!;
print OUT join("\t","ID","Annotation","Confidence",@libNames),"\n";
foreach my$tr (keys %results)
{
    print OUT join("\t",$tr,$results{$tr}{"Ann"},$results{$tr}{"Conf"},@{$results{$tr}{"TPM"}}),"\n";
}
close OUT or die $!;
