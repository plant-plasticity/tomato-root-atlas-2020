#!/usr/bin/perl

if($#ARGV < 0){
  print "usage: Bam2PerBase.pl <bam>\n";
  exit;
}

my $i=1;
my $sam = $ARGV[0];
$sam =~ s/\.bam/\.sam/;
my $perbase = $sam;
$perbase =~ s/\.sam/\.perbase\.bed/;

system "samtools view $ARGV[0] > $sam";

system "/share/brady/people/AlexMason/Alex_scripts/Sam2PerBase_AMedits.sh $sam";
##system "mv big.perbase $perbase";

open(OUT,">$perbase");
open(IN,"<big.perbase");
while(<IN>){
  chomp($_);
  ($chr,$start,$end,$id,$score) = split(/\t/,$_);
  print OUT "$chr\t$start\t$end\tid-$i\t$score\n";
  $i++;
}
close(IN);
close(OUT);

system "rm big.perbase";
system "rm $sam";
