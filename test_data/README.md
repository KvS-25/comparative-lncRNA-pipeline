# Test Data

This directory contains minimal synthetic data to verify the pipeline works
and to illustrate the expected input file formats.

## File Formats

### Reference Transcriptome (FASTA)
- Standard FASTA format
- Headers should be sequence identifiers
- Used as the alignment target for minimap2

### Candidate lncRNA Transcripts (FASTA)
- Standard FASTA format
- Headers should follow StringTie MSTRG naming convention: `>MSTRG.X.Y`
  - X = gene number
  - Y = transcript number
- One file per sample

## Running the Test

1. Update config/config.yaml to point to test data:
```yaml
genome:
  reference: "test_data/reference/test_reference.fasta"

samples:
  directory: "test_data/samples"
  names:
    - TEST1
    - TEST2
```

2. Run the dry run:
```bash
snakemake --dry-run --cores 2
```

3. Run the pipeline:
```bash
snakemake --cores 2 --use-conda
```

## Expected Output Structure

After running, you should see:
```
results/
├── paf/
│   ├── TEST1.paf
│   └── TEST2.paf
├── bed/
│   ├── TEST1.bed
│   └── TEST2.bed
└── multiinter_output.bed
```

## Note
These sequences are synthetic and not biologically meaningful.
They are only intended to verify file formats and pipeline structure.
