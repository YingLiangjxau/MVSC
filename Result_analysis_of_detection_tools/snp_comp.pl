#!/usr/bin/perl -w
use strict;
my $ref=$ARGV[0];  # the ground-truth SNPs
my $file=$ARGV[1];       # the detected SNPs
open IN,$ref;
open OUT,'>',$ARGV[2];      #the output file
my $count=0;
my $i=0;
my @array;
while(<IN>)
{
	$array[$i++]=(split)[0];

}
close IN;
open FIN,$file;
my $temp;
while(<FIN>)
{
	$temp=(split)[0];
	for(my $j=0;$j<$i;$j++)
	{
       if(($temp=~/^$array[$j]$/)&&($array[$j]=~/^$temp$/))
	   {
		   $count++;
		   last;
	   }
	}
}
print OUT"$count\n";
close FIN;
close OUT;
