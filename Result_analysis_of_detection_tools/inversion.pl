#£¡/usr/bin/perl
use strict;
my $standard=$ARGV[0];   #the ground-true inversions
my $file=$ARGV[1];       #the detected inversions
open SIN,$standard;
open FIN,$file;
open OUT,"> $ARGV[2]";   #the output file
my @array;
my $i=0;
my $count;
while(<SIN>)
{
	my @line=split;
    push @array,[@line];
	$i++;
}
while(<FIN>)
{
	my @line=split;
	my $j=0;
    for($j;$j<$i;$j++)
	{
		if($array[$j][0]==$line[0])
		{
			if((($array[$j][1]<$line[1])&&($line[1]<$array[$j][2]))||(($array[$j][1]<$line[2])&&($line[2]<$array[$j][2])))
			{
				$line[4]=1;
			    $count++;
				last;
			}
		}
	}
	if($j==$i)
	{
		$line[4]=0;
	}
    print OUT"@line\n";
}
print OUT"$count\n";
close SIN;
close FIN;
close OUT;
