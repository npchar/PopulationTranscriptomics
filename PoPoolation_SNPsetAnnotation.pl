#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Estimate if SNP is synonymous or not
use Bio::Seq;
use Bio::SeqIO;
use Bio::Index::Fasta;
use Bio::Tools::CodonTable;

# From a mapping on a reference (multiple)
#=========================================
# Estimate
#Flanking region :
# FOr every SNP from PoPoolation2 :
# - convert position
# - find if UTR or CDS (and SYN or NONSYN)
# - flanking region on reference

# Not implemented :
#==================
# In order to produce a solid SNP set I should check as much as possible if there is an INDEL close to the SNP (flanking distance ?) 
# So, I have to determine the ancestral alelic state ! 
# option I have : if allele is > 90% in IPERS / just consider Major allele
# 	- determine state diretly from the Sync file
# 	- determine state from BAM file(s)


#Description :
my $usage ;
$usage= "$0 -s <rc file from PoPoolation2> -i <correspondance> -a <Annotation> -f <Fasta Reference>\n\n" ;
$usage .= "     Description: \n";
$usage .= "             -s <rc file from PoPoolation2> stores the SNP abundance called by population (should be filtered: Indel, low coverage, ...).\n";
$usage .= "             -i <correspondance> Link between transcripts and artifical chromosomes (from concatenation of transcript to speed up previous computionnal steps)\n";
$usage .= "             -a <Annotation> Trinotate report (as a tsv file)\n";
$usage .= "             -f <Fasta Reference> Original file containing individual transcripts\n ";
$usage .= "     Option: \n";
$usage .= "             -c codon table [default:1] (standard)\n";  
$usage .= "Note that this script is working for transcript in strand, if a CDS is predicted by transdecoder as oposite, UTR/CDS nature as well as the Synonymous/Non-Synonymous nature of the variation will be noted as NA\n\n";

# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

my $SNPfile ;
my $help ;
my $ANNOTfile ;
my $CORRESfile ;
my $FASTAfile ;
my $codontable = 1 ;

GetOptions(
           's=s'	  => \$SNPfile,
           'c=s'	  => \$codontable,
	   'a=s'	  => \$ANNOTfile,
	   'f=s'	  => \$FASTAfile,
	   'i=s'	  => \$CORRESfile,
           'h'        => \$help,
          );


#my $OrthogroupFile = "Orthogroups.csv";
if($help){
    print $usage;
    exit 0;
}
unless(defined $SNPfile or defined $ANNOTfile or defined $CORRESfile or defined $FASTAfile){
	print $usage;
	exit 0;
}

# Open Correspondance and store data
#===================================
open IN, "<$CORRESfile" or die "Unable to open \"$CORRESfile\" !\n";
my %DictConcat ;
# Wrapping row by row ($i)
my $count = 0 ;
while(my $i = <IN>){
  chomp $i;
  $count++;
  my @words = split "\t", $i ;
  my $transcript = $words[0] ;
  my $chr = $words[1] ;
  my $start = $words[2] ;
  my $stop = $words[3] ;
  
  $DictConcat{$chr}{$count}{"start"}= $start ;
  $DictConcat{$chr}{$count}{"stop"}= $stop ;
  $DictConcat{$chr}{$count}{"transcript"}= $transcript ;
}
close IN;
# Open Annotation and Store data
#===============================
open IN, "<$ANNOTfile" or die "Unable to open \"$ANNOTfile\" !\n";
my %DictAnnot ;
while(my $i = <IN>){
  chomp $i;
  my @words = split "\t", $i ;
  my $transcript = $words[0] ;
  if($words[5] =~ /^([0-9]+)\-([0-9]+)\[([+-])\]/){
	#my @positionProt = split "-", $words[5] ;
	my $protstart = $1 ;
	my $protstop = $2 ;
	my $strand = $3 ;
# 	print "$transcript $protstart $protstop $strand\n" ;
	$DictAnnot{$transcript}{"CDS"} = "TRUE" ;
	$DictAnnot{$transcript}{"strand"} = $strand ;
	$DictAnnot{$transcript}{"start"} = $protstart ;
	$DictAnnot{$transcript}{"stop"} = $protstop ;
  }
  else{$DictAnnot{$transcript}{"CDS"} = "FALSE" ;}
}
close IN;

# Index Fasta
#===============================
my $IndexFasta ;
unless (-e "${FASTAfile}.inx"){
  print "Indexing ${FASTAfile}\n";
  $IndexFasta = Bio::Index::Fasta->new(-filename => ${FASTAfile}.'.inx',-write_flag => 1);
  $IndexFasta->make_index(${FASTAfile});
}else{ $IndexFasta = Bio::Index::Fasta->new(-filename => ${FASTAfile}.'.inx',-write_flag => 0)}

my $myCodonTable   = Bio::Tools::CodonTable->new( -id => $codontable );

# Open SNP file and evaluate
#===========================
open IN, "<$SNPfile" or die "Unable to open \"$SNPfile\" !\n";
print "chr\tpos\ttrans\tpost\tflan\ttype\tsyn\tAAref\tAAMa\tAAMi\n" ;
while(my $i = <IN>){
	chomp $i;
	
	my @words = split "\t", $i ;
	my $chr = shift @words ;
	my $position = shift @words ;
	if($position eq 'pos' ){ next ;}
	my $embRef = shift @words ;
	my $NumberOfAllele = shift @words ;
	my $AllelicStates = shift @words ;
	my $UnknownCol = shift @words ;
	my $UnknownCol2 = shift @words ;
	my $MajorAlleleStateVector = shift @words ;
	my $MinorAlleleStateVector = shift @words ;
	
	# check if allele is biallelic
	#=============================
	my $Maa ;
	my $Mia ;
	if($NumberOfAllele > 2){
		next ;
	}
	else{
		my @allele = split '/', $AllelicStates ;
		$Maa = $allele[0] ;
		$Mia = $allele[1] ;
	}
	
	# Is coverage sufficient ?
	#=========================
	#
	#Not Yet implemented
	#
	
	#Find transcript Name and Position on transcript :
	#=================================================
	my $transcript ;
	my $TranscriptIndex ;
	my $TranscriptIndexStop ;
	foreach my $p (sort keys %{$DictConcat{$chr}}){
	      if($position <= $DictConcat{$chr}{$p}{"stop"} and  $position >= $DictConcat{$chr}{$p}{"start"}){
		      $transcript = $DictConcat{$chr}{$p}{"transcript"} ;
		      $TranscriptIndex = $DictConcat{$chr}{$p}{"start"} ;
		      $TranscriptIndexStop  = $DictConcat{$chr}{$p}{"stop"} ;
		      last ;
	      }
	}
	my $positionOnTranscript = $position - $TranscriptIndex + 1; # +1 because we do not start at 0 ! ## double checked !
	
	my $flanckingRegion ;
	# retrieving Reference transcript :
	# =================================
	my $SeqObj ;
	if(defined $IndexFasta->fetch($transcript)){ 
		$SeqObj = $IndexFasta->fetch($transcript) ;
# 		$flanckingRegion = "In progress!" ;
	}else{
		$flanckingRegion = "Unable To retrieve $transcript !";
	}

# 	print "$chr $position $transcript $positionOnTranscript $flanckingRegion ";
	
	#Is UTR or CDS ?
	# ==============
	my $codingInformations;
	my $codonPosition ;
	my $synonymous = 'NA' ;
	my $AAMaa = 'NA' ;
	my $AAMia = 'NA' ;
	my $AAref = 'NA' ;
	if($DictAnnot{$transcript}{"CDS"} eq "TRUE" ){
		if($DictAnnot{$transcript}{"strand"} eq "+" and $positionOnTranscript >= $DictAnnot{$transcript}{"start"} and $positionOnTranscript <= $DictAnnot{$transcript}{"stop"}){
# 			print "strand+ ";
			$codingInformations="CDS" ;
			#Determine position on the CDS
			my $positionOnCDS = $positionOnTranscript - $DictAnnot{$transcript}{"start"} + 1 ; # +1 because we do not start at 0 ## double checked !
			
			#Determine if SNP is in 1st, 2nd or 3rd position (Working Only for + strand)
			# if position modulo 3 equal 0 -> SNP is in third position
			my $testpos = 0 ; 
			while( ($positionOnCDS + $testpos) % 3 ne 0){
			      $testpos++ ;
			}
			if($testpos eq 1){ 
				$codonPosition = 2 ;
			} elsif ($testpos eq 2 ) { 
				$codonPosition = 1 ;
			} elsif ($testpos eq 0) { 
				$codonPosition = 3 ;
			}
# 			print "$positionOnCDS $codonPosition ";
			
			# find Fasta and check if synonymous
			if(defined $IndexFasta->fetch($transcript)){ 
# 				$flanckingRegion = "Unable To retrieve $transcript !";
# 			} else {
				my $startCodonPosition ;
				my $stopCodonPosition;
				my $codonRef; 
				my $codonMaa ; 
				my $codonMia ; 
				if($codonPosition eq 1){
					$startCodonPosition = $positionOnTranscript ;
					$stopCodonPosition = $positionOnTranscript + 2 ;
					$codonRef = $SeqObj->trunc($startCodonPosition, $stopCodonPosition) ;
					$codonMaa = $codonRef->seq() ;
					$codonMaa =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/${Maa}$2$3/ ;
					$codonMia = $codonRef->seq() ;
					$codonMia =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/${Mia}$2$3/ ;
				} elsif ($codonPosition eq 2) {
					$startCodonPosition = $positionOnTranscript - 1 ;
					$stopCodonPosition = $positionOnTranscript + 1 ;
					$codonRef = $SeqObj->trunc($startCodonPosition, $stopCodonPosition) ;
					$codonMaa = $codonRef->seq() ;
					$codonMaa =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/$1${Maa}$3/ ;
					$codonMia = $codonRef->seq() ;
					$codonMia =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/$1${Mia}$3/ ;
				} elsif ($codonPosition eq 3) {
					$startCodonPosition = $positionOnTranscript - 2 ;
					$stopCodonPosition = $positionOnTranscript ;
					$codonRef = $SeqObj->trunc($startCodonPosition, $stopCodonPosition) ;
					$codonMaa = $codonRef->seq() ;
					$codonMaa =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/$1$2${Maa}/ ;
					$codonMia = $codonRef->seq() ;
					$codonMia =~ s/^([ATCGatcg])([atcgATCG])([ATCGatcg])$/$1$2${Mia}/ ;
				}
				$AAMaa = $myCodonTable->translate($codonMaa);
				$AAMia = $myCodonTable->translate($codonMia);
				$AAref = $myCodonTable->translate($codonRef);
				
# 				print "$AAMia $AAMaa " ;
				if($AAMia ne $AAMaa){
					$synonymous = "FALSE" ;
				} else {
					$synonymous = "TRUE" ;
				}
			}
		      
		}else{
		      $codingInformations="UTR" ; # It could also be a part of the CDS missed (frameshift producing a stop codon).....
		}
	}else{
		$codingInformations="UTR" ; # It could also be an undected CDS by Transdecoder, or whatever....
	}

	# Retrieve Reference FlanckingRegion
	# ==================================
	unless(defined $IndexFasta->fetch($transcript)){ 
		$flanckingRegion = "NA";
	}elsif ($positionOnTranscript <= 1 or $position >= $TranscriptIndexStop){
		$flanckingRegion = "NA-init";
	}else{
	      $flanckingRegion = $SeqObj->subseq($positionOnTranscript-1, $positionOnTranscript+1)
	}

	# report informations
	print "$chr\t$position\t$transcript\t$positionOnTranscript\t$flanckingRegion\t$codingInformations\t$synonymous\t$AAref\t$AAMaa\t$AAMia\n";
}
close IN;


#  INFOS Codontable ID :
#  =====================
#      'Standard',        #1
#      'Vertebrate Mitochondrial',#2
#      'Yeast Mitochondrial',# 3
#      'Mold, Protozoan, and CoelenterateMitochondrial and Mycoplasma/Spiroplasma',#4
#      'Invertebrate Mitochondrial',#5
#      'Ciliate, Dasycladacean and Hexamita Nuclear',# 6
#      '', '',
#      'Echinoderm Mitochondrial',#9
#      'Euplotid Nuclear',#10
#      '"Bacterial"',# 11
#      'Alternative Yeast Nuclear',# 12
#      'Ascidian Mitochondrial',# 13
#      'Flatworm Mitochondrial',# 14
#      'Blepharisma Nuclear',# 15
#      'Chlorophycean Mitochondrial',# 16
#      '', '',  '', '',
#      'Trematode Mitochondrial',# 21
#      'Scenedesmus obliquus Mitochondrial', #22
#      'Thraustochytrium Mitochondrial', #23
#      'Strict', #24, option for only ATG start
