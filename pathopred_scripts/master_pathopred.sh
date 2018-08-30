
#!/bin/bash



exec_cmd(){ #{{{
    case $VERBOSE in
        yes|1)
        echo -e "\n$*\n"
    esac
    eval "$*"
}
#}}}

usage="
USAGE: $0 SEQ MUT OUTDIR ID
OPTIONS:
  SEQ - FASTA file containing sequence
  MUT - File containing variants in format reference AA, position, altered AA, e.g. A509G
  OUTDIR - Path to save annotated variant file
  ID - Sequence identifier (e.g. Uniprot ID)
"

VERBOSE=0
DEBUG=0

echo "Starting pathopred script"


if [ "$DISPLAY" == "" ];then
    export DISPLAY=localhost:0.0
fi

#argument parser
if [ $# -lt 4 ];then
	echo "$usage"
	exit 1
fi


positionalArgList=()

isNonOptionArg=0
while [ "$1" != "" ]; do
    if [ $isNonOptionArg -eq 1 ]; then
        positionalArgList+=("$1")
        isNonOptionArg=0
    elif [ "$1" == "--" ]; then
        isNonOptionArg=true
    elif [ "${1:0:1}" == "-" ]; then
        case $1 in
            -h | --help) echo "$usage"; exit;;
            -debug|--debug) DEBUG=1;;
            -verbose|--verbose) VERBOSE=1;;
            -*) echo Error! Wrong argument: $1 >&2; exit;;
        esac
    else 
        positionalArgList+=("$1")
    fi
    shift
done

numPositionalArgs=${#positionalArgList[@]}

if [ $numPositionalArgs -ne 4 ];then
  echo "Wrong number of positional arguments, must be 4"
  echo "$usage"
  exit 1
fi

FASTA=${positionalArgList[0]}
VARFILE=${positionalArgList[1]}
OUTDIR=${positionalArgList[2]}
IDNAME=${positionalArgList[3]}

#FASTA=$(readlink -f $FASTA)
#VARFILE=$(readlink -f $VARFILE)
#OUTDIR=$(readlink -f $OUTDIR)

TMPDIR="${OUTDIR}"/tmp/
#create the directories if they do not exist

if [ ! -d "$OUTDIR" ]; then
  mkdir -p $OUTDIR
fi

if [ ! -d "$TMPDIR" ]; then
  mkdir -p $TMPDIR
fi


echo "$FASTA"
echo "$VARFILE"
echo "$OUTDIR"

#let's make sure outdir actually exists 
if [ ! -d "$OUTDIR" ]; then
    echo "$0: $OUTDIR does not exist"
    exit 1	
fi

#Where is the script located?
#rundir=$(dirname $0)
#rundir=$(readlink -f $rundir)
rundir="/Users/alex/pathopred_scripts" #fix for OSX

#Start time
res1=$(/bin/date +%s)

#Write to log file in the out directory..
echo "Attempting to print the script working directory" >> "$OUTDIR"/pathopred_log.txt
echo "$rundir" >> "$OUTDIR"/pathopred_log.txt

#Now we need to run a blast and save this as identifier.xml
echo $(cat $FASTA) >> "$OUTDIR"/pathopred_log.txt
base_name=$(basename "$FASTA" | cut -d. -f1)
echo "$base_name"
echo "Printing base name" >> "$OUTDIR"/pathopred_log.txt
echo "$base_name" >> "$OUTDIR"/pathopred_log.txt

nr_path="$rundir/blastdb_soft/nr"
#blastpath="/usr/bin/blastp"
blastpath="/usr/local/ncbi/blast/bin/blastp"
xml_path="${OUTDIR}/${base_name}.xml"

echo "Running BLAST"

echo "Starting BLAST with xml path $xml_path" >> "$OUTDIR"/pathopred_log.txt
blast_output=$("$blastpath" -query "$FASTA" -db "$nr_path" -evalue=0.001 -out "$xml_path" -outfmt 5 2>&1)
echo "$blast_output" >> "$OUTDIR"/pathopred_log.txt
#For testing, actually just copy the dummy xml in the script folder to this path
#echo "Copying dummy xml" >> "$OUTDIR"/pathopred_log.txt
#echo $(cp $rundir/P26439.xml $xml_path) >> "$OUTDIR"/pathopred_log.txt
#echo "BLAST attempt finished" >> "$OUTDIR"/pathopred_log.txt


#Feed xml into python script
#Store intermediate files in tmpdir
echo $(python $rundir/prepare_alignments.py $FASTA $xml_path $TMPDIR $nr_path)

#Check output to see if successful
if [[ ! -f "$TMPDIR"/"$base_name"_100-90.a3m || ! -f "$TMPDIR"/"$base_name"_90-0.a3m ]]; then
    echo "Alignments do not exist"
    exit 1
fi

echo "Attempting to generate prediction" >> "$OUTDIR"/pathopred_log.txt
#Then feed these into predictor
echo $(python $rundir/predictor.py $FASTA ${TMPDIR}/${base_name}_90-0.a3m ${TMPDIR}/${base_name}_100-90.a3m $VARFILE $OUTDIR $IDNAME 2>&1) >> "$OUTDIR"/pathopred_log.txt

#Make sure the output prediction file exists
if [ ! -f "$OUTDIR"/output_predictions ]; then
    echo "Predictions were not generated"
    exit 1
fi

echo "Exiting pathopred script" >> "$OUTDIR"/pathopred_log.txt
echo "Exiting pathopred script"


#End time
res2=$(/bin/date +%s.%N)
timefile=$OUTDIR/time.txt
runtime=$(echo "$res2 - $res1"|/usr/bin/bc)
echo "0;$runtime" > $timefile
exit 0
