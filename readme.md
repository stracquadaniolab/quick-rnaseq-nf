# quick-rnaseq-nf

![](https://github.com/stracquadaniolab/quick-rnaseq-nf/workflows/build/badge.svg)

A basic and quick workflow for differential expression analysis

## Configuration

- param1: this is the parameter description (default: "hello")
- param2: this is the parameter description (default: "world")
- ...
- paramN: this is the parameter description (default: "flow")

## Running the workflow

### Install or update the workflow

```bash
nextflow pull stracquadaniolab/quick-rnaseq-nf
```

### Run the analysis

```bash
nextflow run stracquadaniolab/quick-rnaseq-nf
```

## Results

- `results/analysis.txt`: file with the analysis results.
- `results/tuning.txt`: file with the parameter tuning.
- ...

## Authors

- Giovanni Stracquadanio
