#  sci-ATAC-seq processing pipeline

## Step 1. Barcode fixation and trimming

```bash
python sc_atac_fastq2bam_SS.py -R1 READ1 -R2 READ2 -O OUTDIR -P PREFIX
```

## Step 2. Mapping

```bash
bowtie2 -p 8 -X 2000 -3 1 -x GENOME -1 READ1 -2 READ2 -S OUTILE 2> LOG
```

## Step 3. Quality filtering

```bash
samtools view -h -f3 -F12 -q10 INBAM | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrU | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 8 - -T temp -o OUTBAM
```

## Step 4. Deduplication

```bash
python sc_atac_true_dedup.py INBAM OUTBAM
```

## Step 5. Count matrix

```bash
python sc_atac_bam2matrix_SS.py -B INBAM -I INDEX -O OUTDIR -P PREFIX -C "auto" -W BED
```
