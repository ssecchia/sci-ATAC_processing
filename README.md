#  sci-ATAC-seq processing pipeline

Please note that most of the code here is taken directly from Darren Cusanovich's fly-atac repository: https://github.com/shendurelab/fly-atac, all credits go to the respective author(s). I made minor changes to adapt the scripts to our cluster environment.

Also please note that while this repository is public, it is intended for internal use by the Furlong lab members.

## Step 0. Preparations

The starting point of this pipeline are the sequencing FASTQ files for READ1 and READ2. Sequencing is performed with four custom primers, two for READ1 and READ2 and two dedicated primers for reading INDEX1 and INDEX2. Our sequencing facility normally performs the bcl2toFASTQ convertion, and the indexes are added to the name of each read. Before starting, confirm that the two indexes (18 bp each) are correctly reported in the read name and separated by the symbol "+".

This pipeline was written for Python2, and it won't work if you are using Python3. Moreover several libraries are required, therefore I suggest creating a dedicated conda environment. I provide a .yml file from which you should be able to recreate the conda environment I use. Make sure to activate the conda environment before starting.

```bash
conda env create --prefix /path/yourenvironmentname/ --file sciATAC.yml

source activate /path/yourenvironmentname/
```

## Step 1. Barcode fixation and trimming

The indexes of each read are checked against the whole space of possible barcode combinations, and they are flagged as exact, edited (if they can be corrected according to the distance rules) or failed matches. Note that no reads are discarded at this stage. Moreover reads are trimmed for the Nextera transposome sequence with Trimmomatic. New FASTQ files are generated for READ1 and READ2, which serve as input for mapping.

```bash
python sc_atac_fastq_fix_trim.py -R1 READ1 -R2 READ2 -O OUTDIR -P PREFIX
```

## Step 2. Mapping

Standard mapping with bowtie2. "GENOME" points to your favorite bowtie genome index. Standard error is redirected to a LOG file of your choosing. The output is a combined BAM file.

```bash
bowtie2 -p 8 -X 2000 -3 1 -x GENOME -1 READ1 -2 READ2 2> LOG | samtools view -bS - > OUTBAM
```

## Step 3. Quality filtering

After mapping, the BAM file is filtered for low quality mapping reads, unwanted chromosomes and reads with failed or ambigous barcodes. The output is a sorted BAM file.

### Drosophila

```bash
samtools view -h -f3 -F12 -q10 INBAM | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrU | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 8 - -T temp -o OUTBAM; samtools index OUTBAM
```

### Human

```bash
samtools view -h -f3 -F12 -q10 INBAM | grep -v MT | grep -v CHR_ | grep -v KI | grep -v GL | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 8 - -T temp -o OUTBAM; samtools index OUTBAM
```

## Step 4. Deduplication

Duplicate reads are removed for each cell barcode.

```bash
python sc_atac_true_dedup.py INBAM OUTBAM; samtools index OUTBAM
```

## Step 5. Cell calling

In this step, the barcodes that exist in the experiment are retrieved and used together with a read count based cutoff to call cells.

### Generate index table of experiment barcodes

```bash
TBC
```

### Count reads per barcode present in the experiment

Using the deduplicated BAM generated above, reads are counted for each barcode. The output is a tab-delimited text file which reports, for each barcode, the read count and the condition / sample assignment (barcodes that don't exist in the experiment are labelled as 'bkgd').

```bash
python sc_atac_barcode_read_counter.py INBAM INDEXTABLE OUTFILE
```

### Read count cutoff

Barcodes that exist in the experiment and that are above the specified cutoff are determined to be cells. The script requires the same PREFIX used for the output report above. The cutoff parameter can be either a numeric value (ex. 500, it retains cells > 500 reads) or "auto", in which case a cutoff based on mixture modelling is applied (using package mclust). The cutoff is applied equally to any condition / sample present in the experiment. The output is a new index table containing only the barcodes surviving the cutoff and a .pdf file showing the read count distribution per condition / sample.

```bash
Rscript sc_atac_cell_cutoff.R PREFIX CUTOFF
```

## Step 6. Library deconvolution

New BAM files, retaining only the barcodes identified as cells above, are generated per condition / sample.

```bash
python sc_atac_library_deconvoluter.py INBAM INDEXTABLE OUTPREFIX .bam
```

## Step 7. Count matrix generation

Reads are counted by cell and by region specified in the bed file. The ouput is a tab-delimited count matrix of regions (rows) by cells (columns).

```bash
python sc_atac_window_counter.py INBAM INDEXTABLE INBED OUTFILE [Include sites with no reads (True/False)]
```
