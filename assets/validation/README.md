# Overview
The workflow below describes how to validate VAPER assemblies.

## Step 1. Create the samplesheet
Create a samplesheet in the manner shown below. Rows with the same sample name will be compared. Rows with the source listed as "truth" will be treated as the _truth_ assembly for accuracy calculations. All other replicate rows per sample will be used for the precision calculations.

`samplesheet.csv`:
```csv
sample,segment,source,assembly
sample01,wg,test,path/to/test/assembly/sample01.fa
sample01,wg,truth,path/to/truth/assembly/sample01.fa
sample02,wg,test-1,path/to/truth/assembly/sample02-1.fa
sample02,wg,test-2,path/to/truth/assembly/sample02-2.fa
sample02,wg,test-3,path/to/truth/assembly/sample02-3.fa
```

In the example above, the **accuracy** of the sample01 test sample would be determined using the truth sequence and the **precision** of each sample02 sample would be determined in a pairwise manner.

## Step 2. Run the validation script
Run the validation script using the command below. Terminal indels are ignored, as these are normally considered artifacts of read mapping or reference selection.
```
./validate.py \
    --input samplesheet.csv \
    --outdir ./ \
    --ignore-termini
```

## Step 3. Examine the output
A summary of the accuracy and precision should be reported to screen. This summary is based on the thresholds defined using the `--min-completeness`, `--min-accuracy`, and `--min-precision` paramaters. Samples that do not meet these thresholds will be reported. This summary, along with the full list of comparisons, can be found in the `outdir`.