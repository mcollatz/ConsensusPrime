# ConsensusPrime
Pipeline to identify ideal consensus regions for primer design.

## Installation & Dependencies
python3.8

pandas

mafft

primer3

clustalx (or another alignment visualisation programm)

## Usage
** Example**

```bash
consensus_prime.py -i /path_to/multifasta.fna --primer3 /path_to/primer3_parameters.txt
```

**Options:**

command | what it does
  ------------- | -------------
-i, --infile          |Multi- or Singe- Fasta file with protein sequences.  [required]
-o, --outdir          |Specifies output directory. Default = .
--delim               |Delimiter char for fasta header. Default = White space
--idpos               |Position of gene ID in fasta header. Zero based. Default = 0
-x, --primer3         | Primer3 parameter file. [required]
-p, --processes       |Number of processes used for predictions. Default = #CPU-cores
-h, --help            |show this message and exit



## We also provide a Docker image for ConsensusPrime
Simply pull and run a ready-to-use image from Dockerhub:  
```bash
docker run -t --rm -v /path/to/dir/with/your/input/files/:/in \
-v /path/to/dir/for/results/:/out \
-u `id -u $USER`:`id -g $USER` \
mcollatz/consensusprime:1.0 \
/consensus_prime.py --infile /in/multifasta.fas --primer3 /in/primer3_parameters.txt --outdir /out
```
