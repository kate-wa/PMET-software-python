#!/bin/bash
set -e

# 22.1.18 Charlotte Rich
# last edit: 7.2.18 - removed the make 1 big fimohits files

# cl_index_wrapper.sh
# mac -> server Version differences
# ggrep = GNU grep, chnage if required
#pmetroot=/Users/lfuler/Documents/PhD/PMET/PMET_index/scripts
# pmetroot=/home/charlotterich/Apps/PMET_index

function usage () {
    cat >&2 <<EOF
USAGE: PMETindexgenome [options] <genome> <gff3> <memefile>

Creates PMET index for Paired Motif Enrichment Test using genome files.
Required arguments:
-r <PMETindex_path>	: Full path of PMET_index. Required.
-i <gff3_identifier> : gene identifier in gff3 file e.g. gene_id=

Optional arguments:
-o <output_directory> : Output directory for results
-n <topn>	: How many top promoter hits to take per motif. Default=5000
-k <max_k>	: Maximum motif hits allowed within each promoter.  Default: 5
-p <promoter_length>	: Length of promoter in bp used to detect motif hits default: 1000
-j <max_jobs>	: Max number of jobs to run at once.  Default: 1
-v <include_overlaps> :  Remove promoter overlaps with gene sequences. AllowOverlap or NoOverlap, Default : AllowOverlap
-u <include_UTR> : Include 5' UTR sequence? Yes or No, default : No
-f <fimo_threshold> : Specify a minimum quality for hits matched by fimo. Default: 0.05

EOF
}

function error_exit() {
    echo "ERROR: $1" >&2
    usage
    exit 1
}

function check_file() {
  file=$1
    name=$2
    flag=$3

    if [ ! -f "$file" ]
    then error_exit "$name is not a file.  $flag flag is required"
    fi
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    # -z string  = True if the length of string is zero.
    # exits if a string has not been specified
    if [ -z "$value" ]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi

}

# set up defaults
topn=5000
maxk=5
maxjobs=1
promlength=1000
fimothresh=0.05
overlap="AllowOverlap"
utr="No"
gff3id='gene_id'
# set up empty variables
pmetroot=
outputdir=
genomefile=
gff3file=
memefile=

# deal with arguments
# if none, exit
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"  >&2
    usage
    exit 1
fi

while getopts ":r:i:h:o:n:k:n:p:f:j:v:u:" options; do
 case $options in
   r) echo "Full path of PMET_index:  $OPTARG" >&2
   pmetroot=$OPTARG;;
   i) echo "GFF3 feature identifier: $OPTARG" >&2
   gff3id=$OPTARG;;
   o) echo "Output directory for results: $OPTARG" >&2
   outputdir=$OPTARG;;
   n) echo "Top n promoter hits to take per motif: $OPTARG" >&2
   topn=$OPTARG;;
   k) echo "Top k motif hits within each promoter: $OPTARG" >&2
   maxk=$OPTARG;;
   p) echo "Promoter length: $OPTARG" >&2
   promlength=$OPTARG;;
   f) echo "Fimo threshold: $OPTARG" >&2
   fimothresh=$OPTARG;;
   j) echo "Max no. parallel jobs: $OPTARG" >&2
   maxjobs=$OPTARG;;
   v) echo "Remove promoter overlaps with gene sequences: $OPTARG" >&2
   overlap=$OPTARG;;
   u) echo "Include 5' UTR sequence?: $OPTARG" >&2
   utr=$OPTARG;;
   \?) echo "Invalid option: -$OPTARG" >&2
   exit 1;;
   :)  echo "Option -$OPTARG requires an argument." >&2
   exit 1;;
 esac
done

shift $((OPTIND - 1))
genomefile=$1
gff3file=$2
memefile=$3


# check things exist
check_set "$pmetroot" "PMET index scripts path" "-r"
check_file "$genomefile" "fasta input file" ""
check_file "$gff3file" "gff3 input file" ""
check_file "$memefile" "meme input file" ""

# check if output directory is specified, if it is check it already exists, if not create it
# and enter that directory
# if a motif_found_final.txt file is already in the folder exit to avoid overwriting data.
if [ -n $outputdir ]; then
echo "Output directory:" $outputdir
  if [ -d $outputdir ]; then
    cd $outputdir
    if [ -d fimo ]; then
      error_exit "Output directory already exists and contains PMET_index data, exiting to avoid overwriting"
    fi
  else
    echo "Output directory does not exist, creating directory"
    mkdir $outputdir
    cp $genomefile $outputdir
    cp $gff3file $outputdir
    cp $memefile $outputdir
    cd $outputdir
  fi
fi


# Now we're ready to start
date

#start off by filtering the .gff3 to gene lines only
if [ ! -f annot.gff3 ]
then
	cp $gff3file annot.gff3
fi

grep -P '\tgene\t' annot.gff3 > genelines.gff3
echo "filtering the .gff3 to gene lines only COMPLETE"

#strip the potential FASTA line breaks. creates genome_stripped.fa
if [ ! -f genome.fa ]
then
	cp $genomefile genome.fa
fi
python3 $pmetroot/strip_newlines.py
echo "strip line breaks COMPLETE"

#create the .genome file which contains coordinates for each chromosome start
samtools faidx genome_stripped.fa
cut -f 1-2 genome_stripped.fa.fai > bedgenome.genome
echo "chromosome coordinates COMPLETE"


#parse up the .bed for promoter extraction
python3 $pmetroot/parse_genelines.py $gff3id
#the python script takes the genelines.gff3 file and makes a genelines.bed out of it

cut -f 4 genelines.bed > universe.txt

#prepare promoter region information
bedtools flank -l $promlength -r 0 -s -i genelines.bed -g bedgenome.genome > promoters.bed
#remove overlapping promoter chunks
if [ $overlap == 'NoOverlap' ]
then
	bedtools subtract -a promoters.bed -b genelines.bed > promoters2.bed
	mv promoters2.bed promoters.bed
fi
#check that we have no split promoters. if so, keep the bit closer to the TSS
python3 $pmetroot/assess_integrity.py
#possibly add 5' UTR
if [ $utr == 'Yes' ]
then
  echo adding UTRs
  python3 $pmetroot/parse_utrs.py $promlength
fi
python3 $pmetroot/parse_promoter_lengths.py promoters.bed
echo "promoter preparation COMPLETE"

#get promoters
bedtools getfasta -fi genome_stripped.fa -bed promoters.bed -s -fo promoters_rough.fa
python3 $pmetroot/parse_promoters.py promoters_rough.fa promoters.bed
echo "promoters COMPLETE"


#now we can actually FIMO our way to victory
fasta-get-markov promoters.fa > promoters.bg
#FIMO barfs ALL the output. that's not good. time for individual FIMOs
#on individual MEME-friendly motif files too
mkdir memefiles
echo "FIMO bg COMPLETE"


python3 $pmetroot/parse_memefile.py $memefile
echo "parse memfile COMPLETE"

python3 $pmetroot/calculateICfrommeme.py
echo "Calc IC from meme COMPLETE"

echo "Start indexing"
runIndexing () {
  fid=$1
  #echo $fid
  bfid=`basename $fid`
  #echo $bfid
  progpath=$4
  fimothresh=$5
  # get all the possible motif hits using fimo
  fimo --text --thresh $fimothresh --verbosity 1 --bgfile promoters.bg $fid promoters.fa > fimo_$bfid
  echo "Within runIndexing loop a fimo COMPLETE"

  # # parse the fimo output to get top n promoters containing top n hits
   python3 $progpath/parse_matrix_n.py fimo_$bfid $2 $3
   echo "parse matrix n.py COMPLETE within runindexing loop"
   rm fimo_$bfid
}
export -f runIndexing
echo "export runIndexing function COMPLETE"

mkdir fimohits

find memefiles -name \*.txt  | parallel  --jobs=$maxjobs "runIndexing {} $maxk $topn $pmetroot $fimothresh"
wait

#there's a lot of intermediate files that need blanking

if [ ! "$genomefile" == "genome.fa" ]
then
rm genome.fa
fi

if [ ! "$gff3file" == "annot.gff3" ]
then
rm annot.gff3
fi


rm bedgenome.genome
rm genelines.bed
rm genelines.gff3
rm -r memefiles
rm genome_stripped.fa
rm genome_stripped.fa.fai
rm promoters.bed
rm promoters.bg
rm promoters.fa
rm promoters_rough.fa
