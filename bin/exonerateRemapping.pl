#!/usr/bin/env perl
use strict;
use warnings;
#Author: Giovanni Bussotti

#
#
# Problems:
# - remove dependencies to absolute path
# - get the 'exonerateExtension' value as parameter instead of using file 'extensionFile' -- see function extensionStrategy()
# - get the 'info_querySizes' value as parameter instead of read it from the 'extensionFile' file -- see function readQuerySizes()
# - ?? function orthologFinderCommandline()
# - 'clusterConfigName' useless
# - use the current folder for '$exonerateOutDir'
# - check '$reportNerId_cmd'
# - pass parameters for 'orthologFinderCommandline()' function

#HELP
my $help  = 0;
foreach my $field (0..$#ARGV){
  $help = 2 if (($ARGV[$field] eq '-h') or ($ARGV[$field] eq '-help') or ($ARGV[$field] eq '--help'));
}
if ($help > 0){
 help_message();
}


#TAKE OPTIONS
acceptedVariableSpace();
my ($exonerateExtension , %info_querySizes , %chrLength , %allQueryNames  , $succesfull , %allExoneratedCoordinates , $querySizes  , $ID_of_mf2 , @clusterHeaderLines , $clusterConfigLine , $currentQuery_fa , $currentTarget_fa , $ortholog , $currentOrtholog_fa , $cmd_orthologFinder );
my ($ner , $queryFile , $targetGenomeFolder , $mf2 , $chr_subseq , $exonerateMinPercentLength , $exonerate_lines_mode  , $exonerate_success_mode  , $clusterName , $experimentName , $pipelineDirName) = options();
my $extensionFile      = "<deprecated>";
my $clusterConfigName  = "<deprecated>";
my $exonerateOutDir    = "./";
## orthologFinderCommandline();

$ID_of_mf2 = `basename $mf2`;
chomp $ID_of_mf2;
$ID_of_mf2 =~ s/.mf2$//;

#CLUSTER
if ($clusterName eq 'on'){
    readClusterConfig();
    doItOnTheCluster();
    exit;
}


#TAKE THE EXTENSION OVER THE BLAST HIT
extensionStrategy ();

#OUT FILES
makeOutFiles();

#READING TRANSCRIPTS
open (I,"<$queryFile") or die "cannot read the input queryFile $!\n";
my @allTranscriptFasta = <I>;
my %allTranscripts = loadMultifastaIntoHash(@allTranscriptFasta);
close I;




print "#exonerating ${mf2}...\n";
open (M,"<$mf2") or die "Error[exonerateRemapping.pl]! Cannot read the $mf2 file $!\n";
my $line = <M>;
my $query;
while ($line) {
  chomp $line;
  my @info_mf2;
  if ($line=~/^(\S+)\s+/){
    $query = $1;
  }
  my %currentHit = ABblastMformat2parser($line);
  push (@info_mf2 , \%currentHit);

  while (defined ($line = <M>) && ($line=~/^$query\s+/)){
    chomp $line;
    my %currentHit = ABblastMformat2parser($line);
    push (@info_mf2 , \%currentHit);
  }

  if (($exonerate_lines_mode eq 'exhaustive') && ($exonerate_success_mode eq 'exhaustive') ){
    @info_mf2 = sort {$a->{"targetName"} cmp $b->{"targetName"}  || $a->{"targetStart"} <=> $b->{"targetStart"} || $a->{"evalue"} <=> $b->{"evalue"}  || $b->{"bitscore"} <=> $a->{"bitscore"} || $b->{"alignlen"} <=> $a->{"alignlen"}} @info_mf2;
  }
  else {
    @info_mf2 = sort {$a->{"evalue"} <=> $b->{"evalue"}  || $b->{"bitscore"} <=> $a->{"bitscore"} || $b->{"alignlen"} <=> $a->{"alignlen"}} @info_mf2;
  }

  $succesfull =  0;
  $ortholog   =  0;
  my ($exonerated_chr , $exonerated_start , $exonerated_end);
  foreach my $hit (0..$#info_mf2){
    last if (($exonerate_success_mode eq 'ortholog')   &&   ($ortholog == 1) );
    last if (($exonerate_lines_mode ne 'exhaustive')   &&  ($hit == $exonerate_lines_mode) );
    last if (($exonerate_success_mode ne 'exhaustive') && ($exonerate_success_mode ne 'ortholog') && ($succesfull == $exonerate_success_mode) );


    #take blast coordinates
    my $b_target_chr   = $info_mf2[$hit]{"targetName"};
    my $b_target_start = $info_mf2[$hit]{"targetStart"};
    my $b_target_end   = $info_mf2[$hit]{"targetEnd"};

    #verify inclusion
    if (defined $exonerated_start){
      if (($b_target_chr eq $exonerated_chr) && ($b_target_start >= $exonerated_start) && ($b_target_end <= $exonerated_end)){next;}
    }

    #exonerate
    my $mf2_currentLine = $info_mf2[$hit]{"fullLine"};
    ($exonerated_chr , $exonerated_start , $exonerated_end) = exonerate($mf2_currentLine);
    if ($exonerated_chr eq '___skip___'){next;}
    next if ((-e 'extended_blastHit.gtf') && (skipEmptyOutput() == 1));

  }
}
#close LOG;
close OUT_GTF;
close OUT_FASTA;
close M;
(system "rm $currentQuery_fa $currentTarget_fa") == 0 or die "Error[exonerateRemapping.pl]! Cannot remove the temporary $currentQuery_fa $currentTarget_fa files $!\n"  if ((defined $currentQuery_fa)&&(-e "$currentQuery_fa")); ##DEBUG
(system "rm $currentOrtholog_fa") == 0 or die "Error[exonerateRemapping.pl]! Cannot remove the temporary$currentOrtholog_fa file $!\n"  if ((defined $currentOrtholog_fa) && (-e "$currentOrtholog_fa"));











sub exonerate {
  my ($line) = @_;
  my %currentHit = ABblastMformat2parser($line);
  my (%exoResult , $minStart , $maxEnd);
  system "rm $currentQuery_fa " if ((defined $currentQuery_fa) && (-e "$currentQuery_fa"));                                ###DEBUG
  system "rm $currentTarget_fa" if ((defined $currentTarget_fa) && (-e "$currentTarget_fa"));                              ###DEBUG
  system "rm $currentOrtholog_fa" if ((defined $currentOrtholog_fa) && (-e "$currentOrtholog_fa"));

  #PRINTING THE QUERY SEQUENCE
  $currentQuery_fa = fileNameGenerator("${exonerateOutDir}/currentQuery_$ID_of_mf2");
  open (Q,">$currentQuery_fa") or die "cannot create the current query fasta file $!\n";
  my $queryName = '>' . $currentHit{'queryName'};
  die "Error! in $query the query  $currentHit{'queryName'} is missing\n" if (! defined $allTranscripts{$queryName});
  print Q "$queryName\n$allTranscripts{$queryName}";
  close Q;
  if (defined $allQueryNames{$currentHit{'queryName'}} ){
    $allQueryNames{$currentHit{'queryName'}}++;
  }
  else{
    $allQueryNames{$currentHit{'queryName'}} = 1;
  }

  #GET TARGET SEQUENCE
  if (defined $info_querySizes{$currentHit{'queryName'}}){$exonerateExtension = $info_querySizes{$currentHit{'queryName'}};}
  my $targetStart   = min2 ($currentHit{'targetStart'} , $currentHit{'targetEnd'});
  my $targetEnd     = max2 ($currentHit{'targetStart'} , $currentHit{'targetEnd'});
  my $extendedStart = $targetStart - $exonerateExtension;
  my $extendedEnd   = $targetEnd + $exonerateExtension;

  my $chr_file = $targetGenomeFolder . "/" . $currentHit{'targetName'};
  $chr_file .= '.fa' unless (-e $chr_file);
  die "Error[exonerateRemapping.pl]! \'$chr_file\' The chromosome name on the BlastOutput and on the targetGenomeFolder must be the same!\n" unless (-e $chr_file);
  ($extendedStart , $extendedEnd) = chr_subseq_sanity($chr_file , $extendedStart , $extendedEnd);
  my $extendedTargetSequence = `$chr_subseq '$chr_file' $extendedStart $extendedEnd`;
  if ($?){
    #print LOG "Error while running; chr_subseq ../../screenings/mouse/assembly_mm9/chr1.fa $extendedStart $extendedEnd\nit returned $extendedTargetSequence\n";
    print "Error[exonerateRemapping.pl]! cannot run $chr_subseq $chr_file $extendedStart $extendedEnd  $!\n";
    #system "rm ${exonerateOutDir}/*";
    exit 1;
  }

  #PRINTING EXTENDED TARGET SEQUENCE
  my $hitName = "$currentHit{'queryName'}_hit$allQueryNames{$currentHit{'queryName'}}";
  $currentTarget_fa = fileNameGenerator("${exonerateOutDir}/currentTarget_$ID_of_mf2");
  open (T,">$currentTarget_fa") or die "cannot open the currentTarget file $!\n";
  print T ">$hitName\n$extendedTargetSequence\n";
  close T;

  #EXONERATE
  #print "exonerate --verbose 0 --showvulgar yes --showcigar no --showsugar no --showalignment no --model est2genome --softmaskquery no --softmasktarget no --bestn 1 $currentQuery_fa $currentTarget_fa\n";  ##DEBUG
  my $exoOutput = `exonerate --verbose 0 --showvulgar yes --showcigar no --showsugar no --showalignment no --model est2genome --softmaskquery no --softmasktarget no --bestn 1 $currentQuery_fa $currentTarget_fa`;
  if ($?){
    print "Error[exonerateRemapping.pl]! cannot run $exoOutput \n $!\n";
    #system "rm ${exonerateOutDir}/*";
    exit 1;
    #print LOG "exonerate didn't manage to align the query $currentHit{'queryName'}\n";
  }

  if ($exoOutput=~/vulgar:\s*(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\S+)\s+(.*)$/){
    my $q         = $1;
    my $q_start   = $2;
    my $q_end     = $3;
    my $q_strand  = $4;
    my $t         = $5;
    my $t_start   = min2($6,$7);
    my $t_end     = max2($7,$6);
    my $t_strand  = $8;
    my $score     = $9;
    my $rest      = $10;

    if ($t_strand eq '-'){
      $extendedTargetSequence = getReverseComplement($extendedTargetSequence);
      ($t_start , $t_end) = reversePosition($t_start , $t_end , length($extendedTargetSequence));
    }

    my $nucleotides;
    my $exonEnd   = $t_start;
    my $exonStart = $t_start;
    my $jump = 0;
    my $totAlignedQuery = 0;

    while ($rest ne "") {
      if ($rest =~ /^(\S+)\s+(\d+)\s+(\d+)(.*)$/) {
	my $label          = $1;
	my $sequenceLength = $2;
	my $databaseLength = $3;
	$rest              = trim($4);

	if (($label eq 'M') and ($databaseLength > 0)) {
	    if ($jump != 0){
		 my %exon;
		 ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
		 $exon{"start"} = $exonStart;
		 $exon{"end"}   = $exonEnd -1;
		 push @{$exoResult{"$hitName"}{"exons"}}, \%exon;
		 ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
		 $exonEnd      += $jump;
		 $exonStart     = $exonEnd;
		 $jump = 0;
	    }
	  $nucleotides .= substr ($extendedTargetSequence ,$exonEnd, $databaseLength);
	  $exonEnd   += $databaseLength;
	  $totAlignedQuery += $sequenceLength;
	}
	if ($label eq 'G'){
	  if ($databaseLength > 0){
	    $nucleotides .= substr ($extendedTargetSequence ,$exonEnd, $databaseLength);
	    $exonEnd     += $databaseLength;
	  }
	  if ($sequenceLength > 0){
	  }
	}
	if  (($label eq '5') or ($label eq '3') or ($label eq 'I')){
	  $jump   += $databaseLength;
	}
      }
      else{die"Error[exonerateRemapping.pl]! parsing error\n";}
    }

    #printing the last exon
    my %exon;
    ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
    $exon{"start"} = $exonStart;
    $exon{"end"}   = $exonEnd -1;
    push @{$exoResult{"$hitName"}{"exons"}}, \%exon;

    #length check and NER
    my $switch = checkLength(length ($allTranscripts{$queryName}) , $totAlignedQuery);
    if (($switch > 0) && ($ner eq 'yes')) {
      ($switch , my $newEexoResult_p  , $nucleotides , $q_strand , $t_strand) = ner_mode($hitName , $extendedTargetSequence , $queryName );
      %exoResult = %{$newEexoResult_p};
    }
    if ($switch > 0){
	#print LOG "$hitName has a too small length. Discarted\n";
	#next;
	return ("___skip___");
    }

    #check that the result it is not identical to some previous one
    my %currentExonCoordinates = exonCoordinates($hitName,$extendedStart,$currentHit{'targetName'},\%exoResult);
    #foreach my $ff(keys %currentExonCoordinates){print "$ff\n";}exit;
    if (checkIdenticalExonerateResult($currentHit{'queryName'},\%currentExonCoordinates) == 1){
	#next;
	return ("___skip___");
    }
    push (@{$allExoneratedCoordinates{$currentHit{'queryName'}}} , \%currentExonCoordinates);


    #print the output
    if ($exonerate_success_mode eq 'ortholog'){
      $currentOrtholog_fa = fileNameGenerator("${exonerateOutDir}/ortholog_$ID_of_mf2");
      open (ORT,">$currentOrtholog_fa") or die "Error[exonerateRemapping.pl]! connot create the current ortholog check file: $currentOrtholog_fa\n$!\n";
      my $nucleotides4ortholog = string2FASTA ($nucleotides);
      print ORT ">$currentHit{'queryName'}";
      print ORT "\n$nucleotides4ortholog\n";
      close ORT;
      $cmd_orthologFinder .= " -query $currentOrtholog_fa ";
      $ortholog = `$cmd_orthologFinder`; if ($?){die "Error[exonerateRemapping.pl] Error While calling orthologFinder.pl with commandline:\n$cmd_orthologFinder\n$!\n";}
      if ($ortholog == 0) {return ("___skip___");}
      else {print STDERR "#ortholog found for $currentHit{'queryName'}\n";}
    }
    $succesfull++;
    ($minStart , $maxEnd) = printExonerateSuccessfulOut(\$hitName , \$nucleotides , \%currentHit , \$q_strand , \%exoResult , \$t_strand , \$extendedStart);


  }
  else {
    #print LOG "$queryName returned an exonerate vulgar output impossible to parse...probable error\n";
      return ("___skip___");
  }
  return ($currentHit{'targetName'} , $minStart , $maxEnd);
}



sub ner_mode {
    my ($hitName , $extendedTargetSequence , $queryName ) = @_;
    my %exoResult;
    my $exoOutput = `exonerate --verbose 0 --showvulgar yes --showcigar no --showsugar no --showalignment no --model ner --softmaskquery no --softmasktarget no --bestn 1 $currentQuery_fa $currentTarget_fa`;
    if ($?){
      print "Error[exonerateRemapping.pl]! cannot run $exoOutput \n $!\n";
      exit 1;
    }
    if ($exoOutput=~/vulgar:\s*(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\S+)\s+(.*)$/){
      my $q         = $1;
      my $q_start   = $2;
      my $q_end     = $3;
      my $q_strand  = $4;
      my $t         = $5;
      my $t_start   = min2($6,$7);
      my $t_end     = max2($7,$6);
      my $t_strand  = $8;
      my $score     = $9;
      my $rest      = $10;
      if ($t_strand eq '-'){
        $extendedTargetSequence = getReverseComplement($extendedTargetSequence);
        ($t_start , $t_end) = reversePosition($t_start , $t_end , length($extendedTargetSequence));
      }
      my $nucleotides;
      my $exonEnd   = $t_start;
      my $exonStart = $t_start;
      my $jump = 0;
      my $totAlignedQuery = 0;

     while ($rest ne "") {
       if ($rest =~ /^(\S+)\s+(\d+)\s+(\d+)(.*)$/) {
      	 my $label          = $1;
      	 my $sequenceLength = $2;
      	 my $databaseLength = $3;
      	 $rest              = trim($4);


      	 if (($label eq 'M') and ($databaseLength > 0)) {
      	   if ($jump != 0){
      	     my %exon;
      	     ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
      	     $exon{"start"} = $exonStart;
      	     $exon{"end"}   = $exonEnd -1;
      	     push (@{$exoResult{"$hitName"}{"exons"}}, \%exon);
      	     ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
      	     $exonEnd      += $jump;
      	     $exonStart     = $exonEnd;
      	     $jump = 0;
      	   }
      	   $nucleotides .= substr ($extendedTargetSequence ,$exonEnd, $databaseLength);
      	   $exonEnd   += $databaseLength;
	   $totAlignedQuery += $sequenceLength;
      	 }
      	 if ($label eq 'G'){
      	   if ($databaseLength > 0){
      	     $nucleotides .= substr ($extendedTargetSequence ,$exonEnd, $databaseLength);
      	     $exonEnd     += $databaseLength;
      	   }
      	   if ($sequenceLength > 0){
      	   }
      	 }
      	 if  (($label eq '5') or ($label eq '3') or ($label eq 'N')){
      	   $jump   += $databaseLength;
      	 }
        }
        else{die"Error[exonerateRemapping.pl]! parsing error\n";}
      }

      #printing the last exon
      my %exon;
      ($exonStart , $exonEnd) = reversePosition($exonStart , $exonEnd , length($extendedTargetSequence))if ($t_strand eq '-');
      $exon{"start"} = $exonStart;
      $exon{"end"}   = $exonEnd -1;
      push @{$exoResult{"$hitName"}{"exons"}}, \%exon;

      #length check
      my $switch = checkLength(length ($allTranscripts{$queryName}) , $totAlignedQuery);
      if ($switch == 0){
	$queryName =~s/^>//;
	my $reportNerId_cmd = "echo \"$queryName\" >> ${pipelineDirName}/experiments/${experimentName}/STDERR/ner_${ID_of_mf2}_transcripts.id";
	(system "$reportNerId_cmd") == 0 or die "Error[exonerateRemapping.pl]! cannot run\n$reportNerId_cmd\n$!\n";
      }
      return ($switch , \%exoResult , $nucleotides , $q_strand , $t_strand );
    }
}
















####FUNCTIONS
sub printExonerateSuccessfulOut {
  my ($hitName_r , $nucleotides_r , $currentHit_r , $q_strand_r , $exoResult_r , $t_strand_r , $extendedStart_r) = @_;
  my $hitName     = ${$hitName_r};
  my $nucleotides = ${$nucleotides_r};
  my %currentHit  = %{$currentHit_r};
  my $q_strand    = ${$q_strand_r};
  my %exoResult   = %{$exoResult_r};
  my $t_strand    = ${$t_strand_r};
  my $extendedStart    = ${$extendedStart_r};


  my $headerName = $hitName;
  $headerName=~s/_hit\d+$/_hit/;
  $headerName .= $succesfull;
  $nucleotides = string2FASTA ($nucleotides);
  print OUT_FASTA ">$headerName";
  print OUT_FASTA "\n$nucleotides\n";
  my $annotations = "query_id \"$currentHit{'queryName'}\"; hitName \"$headerName\"; gene_id \"FAKE__$currentHit{'queryName'}\"; transcript_id \"$currentHit{'queryName'}\";";
  $annotations .= " revCom \" \";" if ($q_strand eq '-');
  my $minStart = 10000000000000000000000000000000000000000000000000000000000000000000000000;
  my $maxEnd   = 0;
  foreach my $exonNumber (0..$#{$exoResult{"$hitName"}{"exons"}}){
    my $start = $exoResult{"$hitName"}{"exons"}[$exonNumber]{"start"} + $extendedStart;
    my $end   = $exoResult{"$hitName"}{"exons"}[$exonNumber]{"end"}   + $extendedStart;
    $minStart = min2 ($start , $minStart);
    $maxEnd   = max2 ($end   , $maxEnd);
    if    ($t_strand eq '+'){ print OUT_GTF "$currentHit{'targetName'}\tBLAST\texon\t$start\t$end\t.\t+\t.\t$annotations\n";}
    elsif ($t_strand eq '-'){ print OUT_GTF "$currentHit{'targetName'}\tBLAST\texon\t$start\t$end\t.\t-\t.\t$annotations\n";}
    else{die "Error[exonerateRemapping.pl]! $t_strand can be either + o -\n";}
  }
  return ($minStart , $maxEnd);
}

sub chr_subseq_sanity{
  my ($chrFile , $start , $end) = @_;
  my ($allSeq , $allLength);

  if (! defined $chrLength{$chrFile}){
    open (CHR,"<$chrFile") or die "cannot open $chrFile \n";
    foreach my $line (<CHR>){
      chomp $line;
      next if (($line=~/>/) or ($line=~/^\s*$/));
      $line =~ s/ //g;
      $allSeq .= $line;
    }
    close CHR;
    $allLength = length ($allSeq);
    $chrLength{$chrFile} = $allLength;
  }

  $start = 1          if ($start < 0);
  $end   = $chrLength{$chrFile} if ($end > $chrLength{$chrFile});
  return ($start , $end);
}
sub reversePosition {
  my ($start,$end,$chunkLength) = @_;
  my $newStart = $chunkLength - $end;
  my $newEnd   = $chunkLength - $start;
  return ($newStart ,$newEnd);
}
sub checkLength{
  my ($q_length , $totAlignedQuery) = @_;
  my $minLength = sprintf("%.0f", ($q_length / 100) * $exonerateMinPercentLength); 
  my $switch = 0;
  $switch++ if ($totAlignedQuery < $minLength);                                                            #print "$q_length $totAlignedQuery  MIN LENGTH!! $minLength\n"; ####DEBUG
  return $switch;
}

sub loadMultifastaIntoHash{
  my (@multifasta) = @_;
  my (%allTheSequences , $sequence , $header);
  my $spy = 1;
  foreach my $line (@multifasta){
    chomp $line;
    next if ($line=~/^\s*$/);
    if ($line=~/(>\S+)/){
      if ($spy == 2){
        $spy = 1;
        $allTheSequences{$header} = $sequence;
        $header   = '';
        $sequence = '';
      }
      $header = $1;
      next;
    }
    if (($line!~/>/) and ($line=~/\w+/)){
      $sequence .= $line;
      $spy = 2;
      next;
    }
  }
  $allTheSequences{$header} = $sequence;
  return %allTheSequences;
}

sub trim {
  my ($string) = @_;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

sub skipEmptyOutput {
  my $filesize = -s "extended_blastHit.gtf";
  if ($filesize == 0){
    (system "rm extended_blastHit.fa extended_blastHit.gtf")==0 or die "Error[eonerateRemapping.pl]! Cannot rm temporary extended_blastHit.gtf and extended_blastHit.fa files\n$!\n";
    return 1;
  }
  else{
    return 0;
  }
}

sub string2FASTA {
    my ($string) = @_;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    my @fake = split (//, "$string");
    my $count = 0;
    my @copy;
    foreach my $letter (@fake){
        $count++;
        push (@copy,$letter);
        if ($count == 100){
            push (@copy, "\n");
            $count = 0;
        }
    }
    chomp $string if ($count == 0);
    $string       = join( "", @copy );
    return ($string);
}

sub exonCoordinates{
  my ($hitName , $extendedStart , $chr , $exoResult_ref) = @_;
  my %exoResult = %{$exoResult_ref};
  my %out;
  foreach my $exonNumber (0..$#{$exoResult{"$hitName"}{"exons"}}){
    my $start = $exoResult{"$hitName"}{"exons"}[$exonNumber]{"start"} + $extendedStart;
                                         my $end   = $exoResult{"$hitName"}{"exons"}[$exonNumber]{"end"}   + $extendedStart;
    #push(@out,"${chr}_${start}_${end}");
    $out{"${chr}_${start}_${end}"}=1;
  }
  return (%out);
}

sub checkIdenticalExonerateResult {
  my ($queryName,$currentExonCoordinates_ref) = @_;
  my %currentExonCoordinates = %{$currentExonCoordinates_ref};
  my $c = 0;
  foreach my $hit (0..$#{$allExoneratedCoordinates{$queryName}}){
      foreach my $currentExonCoordinate (keys %currentExonCoordinates){
	  $c++ if (defined $allExoneratedCoordinates{$queryName}[$hit]{$currentExonCoordinate});
      }
      if ($c == scalar keys (%currentExonCoordinates)){
	  return 1;
      }
      $c=0;
  }
  return 0;
}

sub min2 {
  my ($a, $b) = @_;
  return $b unless (defined $a);
  return $a unless (defined $b);
  return (($a < $b)? $a: $b);
}

sub max2 {
  my ($a, $b) = @_;
  return $b unless (defined $a);
  return $a unless (defined $b);
  return (($a > $b)? $a: $b);
}

sub readQuerySizes{
  open (S,"<$querySizes") or die "Error[exonerateRemapping.pl]!cannot open the gtf file $querySizes \n$!\n";
  foreach my $line (<S>){
    if ($line=~/^(\S+)\s+(\S+)\s*$/){
      my $id   = $1;
      my $size = $2;
      $info_querySizes{$id} = $2;
    }
    else{die "Error[exonerateRemapping.pl]!cannot parse line:\n$line\n";}
  }
  close S;
}


sub getReverseComplement {
  my ($sequence) = @_;
  $sequence = reverse $sequence;
  my $complementedSequence = "";
  for (my $key = 0; $key < length($sequence); $key++) {
    my $char = substr $sequence, $key, 1;
    if    ($char eq "A") {$complementedSequence .= "T";}
    elsif ($char eq "C") {$complementedSequence .= "G";}
    elsif ($char eq "G") {$complementedSequence .= "C";}
    elsif ($char eq "T") {$complementedSequence .= "A";}
    elsif ($char eq "U") {$complementedSequence .= "A";}
    elsif ($char eq "a") {$complementedSequence .= "t";}
    elsif ($char eq "c") {$complementedSequence .= "g";}
    elsif ($char eq "g") {$complementedSequence .= "c";}
    elsif ($char eq "t") {$complementedSequence .= "a";}
    elsif ($char eq "u") {$complementedSequence .= "a";}
    elsif ($char eq "(") {$complementedSequence .= ")";}
    elsif ($char eq ")") {$complementedSequence .= "(";}
    elsif ($char eq "<") {$complementedSequence .= ">";}
    elsif ($char eq ">") {$complementedSequence .= "<";}
    else {
        $complementedSequence .= $char;
    }
  }
  return $complementedSequence;
}


sub help_message {
my $helpMessage = "\nNAME
exonerateRemapping.pl - It accepts the mformat out of blast and it exonerates the hits\n
SYNOPSIS
exonerateRemapping.pl -pipeline_dir -experiment -mf2 -targetGenomeFolder -query [-cluster -exonerate_lines_mode -exonerate_success_mode -exonerateMinPercentLength -chr_subseq -ner]

DESCRIPTION
   * exonerateRemapping.pl takes as input the multi-FASTA query file, the \"chr\" directory contained in \"allGenomeInfo\" relative to the target genome, and the mformat=2 output of blast.
   * It uses exonerate to remap the queries on the target genomes using the blast hit as seed.
   * It is possible to extend over the blast hit editing the -extensionFile
   * It uses chr_subseq to extract the nucleotidic sequences from the genome assemblies.
   * It returns in experiment folder \"EXONERATE_OUT\"  a .gtf annotation file and a multi-FASTA out file

OPTIONS
   * By default chr_subseq program name is \"chr_subseq\". If this is not the case you must provide the field -chr_subseq with the chr_subseq path.
   * The user can filter out an exonerate hit if this does not cover the query for at least a certain threshold. By default it is 70%. You can edit it setting the parameter -exonerateMinPercentLength
   * By default exonerateRemapping.pl will sort the blast output putting on the top of the list the best candidates in terms of e-value, bitscore and coverage.
     The user can set how many blast output lines exonerete will screen before giving up and passing to the next query hits.
     You can edit this by the parameter -exonerate_lines_mode.
     You can choose either an integer [Default:1000] either \"exhaustive\"
     If you choose an integer the exonerateRemapping.pl will screen just the specified number of lines.
     If you choose \"exhaustive\" it will screen all the blast output lines.
   * The user can set how many successful exonerate extension to keep.
     You can edit this by the parameter -exonerate_success_mode
     You can choose either an integer [Default:1], either \"exhaustive\" either \"ortholog\".
     If you choose \"exhaustive\" it will keep all the successful extension
     If you choose \"ortholog\" it will keep just the successful extension that prooved to be the ortholog
   * exonerateRemapping.pl will stop screening the blast output if -exonerate_lines_mode or -exonerate_success_mode is fulfilled
   * The ortholog relationship is established by a reciprocal blast test performed by orthologFinder.pl
   * The user can run exonerateRamapping.pl on the cluster (if available). ExonerateRamapping.pl will generate bash scripts in the CLUSTER_FILES experiment folder.
     The cluster will be called using the command line specified under \"#COMMAND\" in CONFIG/clusterConfig, while the bash scripts will include the header lines specified under \"#COMMAND\" in CONFIG/clusterConfig.
     ExonerateRamapping.pl will wait untill all the jobs return an exit code from the cluster.
   * The -ner option is not used by default. If set to \"yes\" the script will try to use the ner (non-equivalent regions) exonerate model after that est2genome failed. This is not the default, but it is recommended to use as it can handle target genome rich in gaps (marked by Ns). A normal est2genome would fail in aligning regions having N breaking the alignment. If the alignment is too short, the coverage threshold of 70 would discard the hit.
The problem in using the ner option is that the splicing events are not modelled. This means that ner alignment is just an assembly of HSP blast like alignment, rather than a real, refined transcript assembly.

TROUBLESHOOTING
   * Exonerate will extend the blast hit on both direction either using an arbitrary extension, either using the genomic size (intron + exons) of each query as extension value.
     Edit the experiment file \"CONFIG/exonerateExtensionFile\" to include for each query \"queryName genomicSize\", either just state the arbitrary extension with the following syntax:  __arbitrary_extension__XXX
   * To make the annotation output in  a standard .gtf format both gene_id and transcript_id fields are provided. However the gene_id is a fake name, this is because if the imput is just a multiFASTA file, the information about the query's gene_id is not provided. The transcript_id field will be just the name of the query
";
print "$helpMessage\n\n\n";
exit;
}


sub options {
  my ($ner , $queryFile , $targetGenomeFolder , $mf2 , $chr_subseq , $exonerateMinPercentLength , $exonerate_lines_mode  , $exonerate_success_mode  , $clusterName , $experimentName , $pipelineDirName);
  my $spyQuery                     = 1;
  my $spyTargetGenomeFolder        = 1;
  my $spyMf2                       = 1;
  my $spyChr_subseq                = 1;
  my $spyExonerateMinPercentLength = 1;
  my $spyCluster                   = 1;
  my $spyExperiment                = 1;
  my $spyPipelineDir               = 1;
  my $spyExonerate_lines_mode      = 1;
  my $spyExonerate_success_mode    = 1;
  my $spyNer                       = 1;

  foreach my $field (0..$#ARGV){
    if ($ARGV[$field] eq '-ner'){
	$ner = $ARGV[1+$field];
	$spyNer = 2;
	next;
    }
    if ($spyNer == 2){
	$spyNer = 3;
	next;
    }
    if ($ARGV[$field] eq '-experiment'){
	$experimentName = $ARGV[1+$field];
	$spyExperiment = 2;
	next;
    }
    if ($spyExperiment == 2){
	$spyExperiment = 3;
	next;
    }
    if ($ARGV[$field] eq '-pipeline_dir'){
      $pipelineDirName = $ARGV[1+$field];
      $spyPipelineDir = 2;
      next;
    }
    if ($spyPipelineDir == 2){
      $spyPipelineDir = 3;
      next;
    }
    if ($ARGV[$field] eq '-cluster'){
	$clusterName = $ARGV[1+$field];
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
	next;
    }
    if ($ARGV[$field] eq '-exonerate_success_mode'){
	$exonerate_success_mode = $ARGV[1+$field];
	$spyExonerate_success_mode = 2;
	next;
    }
    if ($spyExonerate_success_mode == 2){
      $spyExonerate_success_mode = 3;
      next;
    }
    if ($ARGV[$field] eq '-exonerate_lines_mode'){
	$exonerate_lines_mode = $ARGV[1+$field];
	$spyExonerate_lines_mode = 2;
	next;
    }
    if ($spyExonerate_lines_mode == 2){
      $spyExonerate_lines_mode = 3;
      next;
    }
     if ($ARGV[$field] eq '-exonerateMinPercentLength'){
	$exonerateMinPercentLength = $ARGV[1+$field];
	$spyExonerateMinPercentLength = 2;
	next;
    }
    if ($spyExonerateMinPercentLength == 2){
      $spyExonerateMinPercentLength = 3;
      next;
    }
    if ($ARGV[$field] eq '-chr_subseq'){
	$chr_subseq = $ARGV[1+$field];
	$spyChr_subseq = 2;
	next;
    }
    if ($spyChr_subseq == 2){
      $spyChr_subseq = 3;
      next;
    }
    if ($ARGV[$field] eq '-mf2'){
      $mf2 = $ARGV[1+$field];
      $spyMf2 = 2;
      next;
    }
    if ($spyMf2 == 2){
      $spyMf2 = 3;
      next;
    }
    if ($ARGV[$field] eq '-targetGenomeFolder'){
      $targetGenomeFolder = $ARGV[1+$field];
      $spyTargetGenomeFolder = 2;
      next;
    }
    if ($spyTargetGenomeFolder == 2){
      $spyTargetGenomeFolder = 3;
      next;
    }
   if ($ARGV[$field] eq '-query'){
      $queryFile = $ARGV[1+$field];
      $spyQuery = 2;
      next;
    }
    if ($spyQuery == 2){
      $spyQuery = 3;
      next;
    }
  }

  die "Error[exonerateRemapping.pl]! Provide the -mf2 parameter. Please indicate the blast output\n"                  if ($spyMf2 != 3);
  die "Error[exonerateRemapping.pl]! Provide the -query parameter. Please indicate the multi-FASTA query file\n"      if ($spyQuery != 3);
  die "Error[exonerateRemapping.pl]! Provide the -targetGenomeFolder parameter.Please indicate the \"chr\" folder\n"  if ($spyTargetGenomeFolder != 3);
#  die "Error[exonerateRemapping.pl]! Please provide the -experiment parameter\n"   if (! defined $experimentName);
#  die "Error[exonerateRemapping.pl]! Please provide the -pipeline_dir parameter\n" if (! defined $pipelineDirName);
  die "Error[exonerateRemapping.pl]! -ner parameter can be either yes or no\n"     if ((defined $ner)&&($ner ne 'no')&&($ner ne 'yes'));
  $ner = 'no'                     if (! defined $ner);
  $chr_subseq = 'chr_subseq'      if (! defined $chr_subseq);
  $clusterName = 'off'            if (! defined $clusterName);
  $exonerateMinPercentLength = 70 if (! defined $exonerateMinPercentLength);
  if ((defined $exonerate_lines_mode) && ($exonerate_lines_mode ne "exhaustive") && ($exonerate_lines_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_lines_mode modes are either \"exhaustive\" either an integer number of iteration\n";
  }
  if ((defined $exonerate_success_mode) && ($exonerate_success_mode ne "exhaustive") && ($exonerate_success_mode ne "ortholog") && ($exonerate_success_mode =~/\D/)){
      die "Error[exonerateRemapping.pl]! The supported exonerate_success_mode modes are either \"exhaustive\" either \"ortholog\" either  an integer number of iteration\n";
  }
  $exonerate_lines_mode   = 1000  if (! defined $exonerate_lines_mode);
  $exonerate_success_mode = 1     if (! defined $exonerate_success_mode);


  return ($ner , $queryFile , $targetGenomeFolder , $mf2 , $chr_subseq , $exonerateMinPercentLength , $exonerate_lines_mode  , $exonerate_success_mode   , $clusterName , $experimentName , $pipelineDirName);
}

sub acceptedVariableSpace {
  my %space = ('-ner' => 1 , '-pipeline_dir' => 1 , '-experiment' => 1  , '-cluster' => 1  , '-exonerate_lines_mode' => 1  , '-exonerate_success_mode' => 1  , '-exonerateMinPercentLength' => 1  , '-chr_subseq' => 1  , '-mf2' => 1  , '-targetGenomeFolder' => 1  , '-query' => 1 );
  foreach my $field (0..$#ARGV){
    if (($ARGV[$field] =~/^-/) && (! defined $space{"$ARGV[$field]"})){
      print "Error[exonerateRemapping.pl]! $ARGV[$field] it is not a valid parameter\n";
      help_message();
    }
  }
}

sub extensionStrategy {
#    open (S,"<$extensionFile") or die "Error! cannot read $extensionFile $!\n";
#    my $lineCounter = 0;
#    foreach my $line (<S>){
#	next if ($line=~/^\s*$/);
#	last if ($lineCounter == 1);
#	if ($line=~/__arbitrary_extension__(\d+)/){
#	    $exonerateExtension = $1;
#	}
#	$lineCounter = 1;
#    }
#    close S;

    $exonerateExtension = 20000;

#    if (! defined $exonerateExtension) {
#	$querySizes = $extensionFile;
#	readQuerySizes();
#    }
}

sub readClusterConfig {
    open (CC,"<$clusterConfigName") or die "Error! Cannot read $clusterConfigName \n$!\n";
    my $spyHeader = 0;
    foreach my $line (<CC>){
	chomp $line;
	next if ($line=~/^\s*$/);
	$clusterConfigLine = $line if (($line=~/##SCRIPT##/) and ($spyHeader == 0));

	if ($line=~/^#HEADER/){$spyHeader = 1;}
	push (@clusterHeaderLines , $line) if ($spyHeader == 1);
    }
    close CC;
}

sub makeOutFiles {
    #open (LOG      ,">log") or die "cannot create the log$!\n";
    open (OUT_FASTA,">${exonerateOutDir}/${ID_of_mf2}.fa")  or die "Error[exonerateRemapping]! Cannot create the output fasta file$!\n";
    open (OUT_GTF  ,">${exonerateOutDir}/${ID_of_mf2}.ex.gtf") or die "Error[exonerateRemapping]! Cannot create the output gtf file$!\n";
}

sub doItOnTheCluster {
    my $clusterFileDir = "${pipelineDirName}/experiments/${experimentName}/CLUSTER_FILES/";
    my $spyClusterConfig = 1;
    my $spyCluster      = 1;

    my $shellCmdLine = "perl " . $0 . " ";
    foreach my $f (@ARGV){
    if ($f eq '-cluster_config'){
	$spyClusterConfig = 2;
	next;
    }
    if ($spyClusterConfig == 2){
	$spyClusterConfig = 3;
	next;
    }

    if ($f eq '-cluster'){
	$spyCluster = 2;
	next;
    }
    if ($spyCluster == 2){
	$spyCluster = 3;
	next;
    }
    $shellCmdLine .= " $f";
    }

    #create the script
    my $scriptName = "${clusterFileDir}/${ID_of_mf2}.ex.sh";
    open (S,">$scriptName") or die "Error[exonerateRemapping.pl]! Cannot create the cluster script for ${ID_of_mf2}  $!\n";
    foreach my $h (@clusterHeaderLines){
	print S "$h\n";
    }
    print S $shellCmdLine . "\n\n";
    print S "touch ${clusterFileDir}/${ID_of_mf2}.ex."."___exonerate_done___" . "\n";
    close S;

    #do the cmd
    my $cmd = $clusterConfigLine;
    $cmd =~s/##SCRIPT##/$scriptName/;
    system "$cmd";  if ($?){die "Error[exonerateRemapping.pl]! Error message returned with:\n$cmd\n$!";}
}

sub fileNameGenerator{
  my ($nameRoot) = @_;
  my $tmp_name_counter = 0;
  my $tmp_name;
  while (!$tmp_name || -f $tmp_name) {
    $tmp_name_counter++;
    $tmp_name = "${nameRoot}_$$".".${tmp_name_counter}.fa";
  }
  return $tmp_name;
}

sub ABblastMformat2parser{
  my ($line) = @_;
  chomp $line;
  my %mf2Line;
  my @fields = split (/\s+/,$line);
  $mf2Line{'queryName'}   = $fields[0];
  $mf2Line{'targetName'}  = $fields[1];
  $mf2Line{'evalue'}      = $fields[2];
  $mf2Line{'bitscore'}    = $fields[4];
  $mf2Line{'score'}       = $fields[5];
  $mf2Line{'alignlen'}    = $fields[6];
  $mf2Line{'pcident'}     = $fields[10];
  $mf2Line{'pcpos'}       = $fields[11];
  $mf2Line{'strand'}      = $fields[19];
  $mf2Line{'targetStart'} = $fields[20];
  $mf2Line{'targetEnd'}   = $fields[21];
  $mf2Line{'fullLine'} = $line;
  return %mf2Line;
}

#
# Bad dependency to: '${pipelineDirName}/experiments/${experimentName}/RNAmapping_pipeline.txt'
#
sub orthologFinderCommandline {
  open (P,"<${pipelineDirName}/experiments/${experimentName}/RNAmapping_pipeline.txt") or die "Error[exonerateRemapping.pl]! cannot read the RNAmapping_pipeline.txt\n$!\n";
  my ($blast_ort , $strategy_ort , $xdformat_ort , $referencegenome_ort , $querygtf_ort);
  foreach my $line (<P>){
    chomp $line;
    if ($line=~/^BLAST=(\S+)$/)          {$blast_ort = $1;}
    if ($line=~/^STRATEGY=(\S+)$/)       {$strategy_ort = $1;}
    if ($line=~/^XDFORMAT=(\S+)$/)       {$xdformat_ort = $1;}
    if ($line=~/^REFERENCEGENOME=(\S+)$/){$referencegenome_ort = $1;}
    if ($line=~/^QUERYGTF=(\S+)$/)       {$querygtf_ort = $1;}
    if ($line=~/^############$/)         {last;}
  }
  close P;
  $cmd_orthologFinder  = "${pipelineDirName}/scripts/orthologFinder.pl -experiment $experimentName -pipeline_dir $pipelineDirName -blast $blast_ort -blast_strategy $strategy_ort  -xdformat $xdformat_ort -reference_genome $referencegenome_ort -query_gtf $querygtf_ort";
}
