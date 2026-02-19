# Lough_Catherine_Shotgun_Analysis


We hereafter parsed all QC reads through the Holi pipeline for taxonomic profiling and DNA
damage estimation (for details about version of tools and options set see SI). To increase
resolution and sensitivity of our taxonomic assignments, we supplemented the RefSeq (92
excluding bacteria) and the nucleotide database (NCBI) with a recently published Arctic-boreal
plant database (PhyloNorway). All alignments were hereafter merged using samtools and sorted
by coordinate and parsed through filterBAM reassign and filter functions to refine alignments
and generate reference-wide statistics. Cytosine deamination frequencies were then estimated
using metaDMG, by first finding the lowest common ancestor across all possible alignments for
The copyright holder for this preprint this version posted October 8, 2025. ; https://doi.org/10.1101/2025.10.08.681091 doi: bioRxiv preprint
(which was not certified by peer review) is the author/funder, who has granted bioRxiv a license to display the preprint in perpetuity. It is made
. CC-BY-NC-ND 4.0 International license available under a each read and then calculating the nucleotide misincorporations (deamination frequencies) due to
DNA damage at each taxonomic level (See SI for details on damage filtering). In parallel, we
computed the mean read length as well as number of reads per taxonomic nodes. Using the DNA
damage for plant species with >500 reads allowed us to create a DNA damage model which we
used to filter all eukaryotic taxa as described in the SI

###############################################################################
# Holi pipeline 
###############################################################################

############################
# 0) INSTALLATION (TOOLS)  #
############################
# Step 0.1 Install with conda/mamba (recommended). Run manually:
#
#   conda create -y -n holi-pipe -c conda-forge -c bioconda \
#     adapterremoval fastq-tools sga bowtie2 samtools pigz gawk coreutils sed
#   conda activate holi-pipe
#
# Tools used:
#   AdapterRemoval, fastq-grep, sga, bowtie2, samtools, awk/sed/coreutils, gzip/pigz

########################
# 1) CONFIG (EDIT ME)  #
########################

# Step 1.1 Threading
THREADS_MAP=80
THREADS_SGA=30
THREADS_SAMTOOLS=30

# Step 1.2 Temporary directory (portable)
TMPDIR="${TMPDIR:-./tmp}"
mkdir -p "$TMPDIR"

# Step 1.3 Database root (EDIT THIS)
# Use a relative path for portability, e.g. DB_ROOT="../databases"
DB_ROOT="${DB_ROOT:-./databases}"

# Step 1.4 Database globs/prefixes (EDIT THESE TO MATCH YOUR SETUP)
# Each entry must expand to Bowtie2 index prefixes (what you pass to bowtie2 -x).
DB_NORPLANT_GLOB="${DB_NORPLANT_GLOB:-$DB_ROOT/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.?}"
DB_NT_GLOB="${DB_NT_GLOB:-$DB_ROOT/ncbi_nt/nt.?}"
DB_VERT_OTHER_GLOB="${DB_VERT_OTHER_GLOB:-$DB_ROOT/refseq/vert_other/vert_other.?}"
DB_VERT_MAM_GLOB_1="${DB_VERT_MAM_GLOB_1:-$DB_ROOT/refseq/vert_mam/vert_mam.?}"
DB_VERT_MAM_GLOB_2="${DB_VERT_MAM_GLOB_2:-$DB_ROOT/refseq/vert_mam/vert_mam.??}"
DB_INVERT_GLOB="${DB_INVERT_GLOB:-$DB_ROOT/refseq/invert/invert.?}"

# If viral_fungi_archaea is a bowtie2 index prefix at your site, keep as-is.
# If it's a fasta (not indexed), you must build a bowtie2 index and point to that prefix instead.
DB_VIRAL_FUNGI_ARCHAEA="${DB_VIRAL_FUNGI_ARCHAEA:-$DB_ROOT/refseq/microbe/viral_fungi_archaea.fa}"

# "Second round" DBs, integrated here
DB_GTDB="${DB_GTDB:-$DB_ROOT/microbial_dbs/GTDB/r89/bowtie2/gtdb_r89}"
DB_PLANT_GLOB="${DB_PLANT_GLOB:-$DB_ROOT/refseq_23dec2020/plant/plant.?}"

########################
# 2) SANITY CHECKS     #
########################
# Step 2.1 Ensure required commands exist
req_cmds=(AdapterRemoval fastq-grep sga bowtie2 samtools awk sed gzip)
for c in "${req_cmds[@]}"; do
  command -v "$c" >/dev/null 2>&1 || {
    echo "ERROR: Missing required command: $c"
    echo "Install tools (see installation section) and try again."
    exit 1
  }
done

#############################################
# 3) ADAPTERREMOVAL ON ALL *.fastq.gz       #
#############################################

# Step 3.1 AdapterRemoval on each input file (parallel)
for file in *.fastq.gz; do
  AdapterRemoval \
    --file1 "$file" \
    --mm 3 \
    --minlength 30 \
    --basename "$file" \
    --trimns \
    --trimqualities \
    --minquality 30 &
done

# Step 3.2 Wait for AdapterRemoval jobs
wait

########################################################
# 4) RENAME *.truncated OUTPUTS TO *.fq                 #
########################################################

# Step 4.1 Convert "something.fastq.gz.truncated..." to "something.fq"
for infile in *.truncated; do
  bname="$(basename "$infile")"
  bname2="$(echo "$bname" | sed 's/.fastq.gz.truncated*/.fq/')"
  echo "Renaming: $bname -> $bname2"
  mv "$bname" "$bname2"
done

########################################################
# 5) PER-SAMPLE PROCESSING + INTEGRATED MAPPING         #
########################################################

for infile in ./*.fq; do
  # Step 5.1 Define sample names
  bname="$(basename "$infile")"                # e.g. sample.fq
  bprefix="$(echo "$bname" | sed 's/.fq*//')"  # e.g. sample
  sample_dir="${bprefix}_holi"

  echo "============================================================"
  echo "Processing sample: $bname"
  echo "Output folder:     $sample_dir"
  echo "============================================================"

  # Step 5.2 Create and enter sample directory
  mkdir -p "./$sample_dir"
  cd "./$sample_dir"

  ############################################
  # Step 5.3 Filtering with fastq-grep        #
  ############################################

  echo "Step 5.3.1 Removing poly A tails"
  fastq-grep -v "AAAAA$" "../$bname" > "kmer_$bname"

  echo "Step 5.3.2 Removing reverse complemented A tails"
  fastq-grep -v "^TTTTT" "kmer_$bname" > "kmer2_$bname"

  echo "Step 5.3.3 Removing remnants adapter sequence 1"
  fastq-grep -v "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" "kmer2_$bname" > "adap1_kmer2_$bname"

  echo "Step 5.3.4 Removing remnants adapter sequence 2"
  fastq-grep -v "ATCTCGTATGCCGTCTTCTGCTTG" "adap1_kmer2_$bname" > "adap2_kmer2_$bname"

  ############################################
  # Step 5.4 SGA preprocess + dedup           #
  ############################################

  echo "Step 5.4.1 sga preprocess"
  sga preprocess --dust-threshold=1 -m 30 "adap2_kmer2_$bname" -o "adap2_kmer2_$bname.pp.fq"

  echo "Step 5.4.2 sga index"
  sga index --algorithm=ropebwt --threads="$THREADS_SGA" "adap2_kmer2_$bname.pp.fq"

  echo "Step 5.4.3 sga filter (rmdup)"
  sga filter --threads="$THREADS_SGA" --no-kmer-check "adap2_kmer2_$bname.pp.fq" -o "adap2_kmer2_$bname.pp.rmdup.fq"

  echo "Step 5.4.4 Read length distribution"
  awk 'NR%4==2 {print length($0)}' "adap2_kmer2_$bname.pp.rmdup.fq" \
    | sort -n | uniq -c > "adap2_kmer2_$bname.pp.rmdup.fq.read_length.txt"

  ############################################
  # Step 5.5 Integrated mapping to ALL DBs    #
  ############################################

  # Step 5.5.1 Define a mapping helper
  map_db () {
    local db_prefix="$1"
    echo "Mapping against: $db_prefix"
    bowtie2 --threads "$THREADS_MAP" -k 1000 -x "$db_prefix" \
      -U "adap2_kmer2_$bname.pp.rmdup.fq" --no-unal \
      | samtools view -bS - > "${bprefix}.$(basename "$db_prefix").bam"
  }

  # Step 5.5.2 Map against every DB group
  for DB in $DB_NORPLANT_GLOB;       do map_db "$DB"; done
  for DB in $DB_NT_GLOB;            do map_db "$DB"; done
  for DB in $DB_VERT_OTHER_GLOB;    do map_db "$DB"; done
  for DB in $DB_VERT_MAM_GLOB_1;    do map_db "$DB"; done
  for DB in $DB_VERT_MAM_GLOB_2;    do map_db "$DB"; done
  for DB in $DB_INVERT_GLOB;        do map_db "$DB"; done
  for DB in $DB_VIRAL_FUNGI_ARCHAEA;do map_db "$DB"; done
  for DB in $DB_GTDB;               do map_db "$DB"; done
  for DB in $DB_PLANT_GLOB;         do map_db "$DB"; done

  ############################################
  # Step 5.6 Merge + Sort             #
  ############################################

  echo "Step 5.6.1 Merging all BAMs into one"
  samtools merge --verbosity 5 -@ "$THREADS_SAMTOOLS" "${bprefix}.allDBs.merged.bam" ./*.bam

  echo "Step 5.6.2 Sorting merged BAM (coordinate sort)"
  samtools sort -@ "$THREADS_SAMTOOLS" \
    -T "$TMPDIR/${bprefix}_sorttmp" \
    -o "${bprefix}.allDBs.merged.sorted.bam" \
    "${bprefix}.allDBs.merged.bam"

  echo "Step 5.6.3 Indexing sorted BAM"
  samtools index "${bprefix}.allDBs.merged.sorted.bam"

  ############################################
  # Step 5.7 Cleanup (optional)               #
  ############################################

  echo "Step 5.7 Cleaning intermediates (optional)"
  rm -f kmer* adap1* "adap2_kmer2_$bname" "adap2_kmer2_$bname.pp.fq" \
        *wt *sai *rmdup.discard.fa || true
  rm -f ./*.bam "${bprefix}.allDBs.merged.bam" || true

  ############################################
  # Step 5.8 metaDMG note                     #
  ############################################

  echo "Step 5.8 metaDMG input ready:"
  echo "  $(pwd)/${bprefix}.allDBs.merged.sorted.bam"
  echo "These merged, coordinate-sorted BAM files were hereafter parsed to metaDMG."

  # Step 5.9 Return to project root
  cd ..
done

echo "DONE"
