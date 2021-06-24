# ConsensusPrime
Pipeline to identify ideal consensus regions for primer design.

## Installation & Dependencies
Python3.8 in miniconda.
https://docs.conda.io/en/latest/miniconda.html

Install MAFFT, Primer3 and ClustalX (optional for alignment visualzation)
```bash
sudo apt-get update -y
sudo apt-get install -y mafft
sudo apt-get install -y primer3
sudo apt-get install -y clustalx
```

Create and activate new Python environment for ConsensusPrime
```bash
conda create -n consensus_prime
source activate consensus_prime
```
Install pandas
```bash
conda install pandas
```



## Usage
**Example**

```bash
python3.8 /path_to/consensus_prime.py -infile /path_to/multifasta.fna --primer3 /path_to/primer3_parameters.txt
```

**Options:**

command | what it does
  ------------- | -------------
-i, --infile              |Multi-Fasta file with gene sequences.  [required]
-x, --primer3             |Primer3 parameter file. [required]
-o, --outdir              |Specifies output directory. Default = .
-t, --threads             |Number of threads used by MAFFT. Default = -1 (all)
-k, --keepduplicates      |Keep duplicate sequences. Default = False
-c, --consensusthreshold  |Removes sequences with higher gap to sequence ratio. Default = 0.2
--primers                 |Known primers for visualisation in the final alignment in multifasta format.
--negativesequences       |File with sequences that get their consensus sequence added to the final alignment in multifasta format.
-h, --help                |show this message and exit



## We also provide a Docker image for ConsensusPrime
Simply pull and run a ready-to-use image from Dockerhub:

**Show --help**
```bash
docker run mcollatz/consensusprime:1.0
```

**Example**
```bash
docker run -t --rm -v /path/to/dir/with/your/input/files/:/in \
-v /path/to/dir/for/results/:/out \
-u `id -u $USER`:`id -g $USER` \
mcollatz/consensusprime:1.0 \
/consensus_prime.py --infile /in/multifasta.fas --primer3 /in/primer3_parameters.txt --outdir /out
```
