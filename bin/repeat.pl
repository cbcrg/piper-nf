#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

my $f = $ARGV[0];
my $inputGtf   = $ARGV[1];
my $repCov   = $ARGV[2];

# But if it is a file set the repeatCoverage per sequence!
# Or if no threshold given read query coverage

if (scalar @ARGV < 3){
  die "USAGE:
        perl repeatCoverage.pl inputFasta inputGtf repeatCoverageThreshold
";
}

# Based on...
# my $result = 0;
# $result++ while($string =~ m/\p{Uppercase}/g);    

#Take all repeat coverage
  my %bh;
  my $name = basename($f);
  chomp $name;
  $name=~s/\.fa$//;

  open (O1, ">repeatCoverage_${name}.id") or die "Error [repeatCoverage.pl]! Cannot create repeatCoverage_${name}.id  $! \n";
  open(F,"<$f") or die "Error [repeatCoverage.pl]! Cannot open $f $!\n";
  my ($fraction , $id , @seq);
  my $size    = 0;
  my $rep     = 0;
  my $starter = 0;
  while(my $line = <F>){
    chomp $line;
    #check
    if (($line=~/>/) && ($starter == 1)){
        $fraction = ($rep/$size)*100;
        print O1 "$fraction\n";
        if ($fraction <= $repCov){
            #print all
            printFunc($id , $name , \@seq , $inputGtf);
        }
        @seq = ();
        $id = '';
        $rep = 0;
        $size= 0;
    }
    $starter = 1;
    if($line=~/>(.+)/){
        $id = $1;
        print O1 "$id\t";
    }
    else {
      push (@seq,"$line\n");
      $rep++ while($line =~ m/\p{Lowercase}/g);
      $size += map $_, $line =~ /(\S)/gs;
    }
  }

#dump
  $fraction = ($rep/$size)*100;
  print O1 "$fraction\n";
  if ($fraction <= $repCov){
      printFunc($id , $name , \@seq ,$inputGtf);
  }

  close F;
  close O1;

  exit;

sub printFunc {
    my ($id , $name , $refseq , $gtf) = @_;
    my @seq = @{$refseq};
        system "grep -P \"$id\\S;\" ${gtf} >> rep${repCov}.ex.gtf";
        open (O2,">>rep${repCov}.fa") or die "Error [repeatCoverage.pl]! Cannot create rep${repCov}.fa $!\n";
        print O2 ">$id\n";
        foreach my $s (@seq){print O2 "$s";}
        close O2;
}
