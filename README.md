# ConsensusPrime
Pipeline to identify ideal consensus regions from homologue sequences for primer design.

## System Requirements

The pipeline is developed and testet for:

Software | Version
  ------------- | -------------
Linux   | 20.04
Mafft   | 7.453
Primer3 | 2.5.0

On the hardware side, the alignments with mafft are the bottleneck. The length and number of sequences play a decisive role. Alignments with a few hundred sequences of moderate length can be calculated in a few seconds to minutes even on simple laptops. Larger alignments require a system with more RAM.

The duration of the installation is a few minutes and depends among other things on whether various requirements are already installed, such as conda.

The processing time of the sample dataset is about 15 sec on a laptop with Intel(R) Core(TM) i7-10750H CPU @ 2.60GHz and 16 GB RAM.

## Installation & Dependencies
Install Python3.8 in miniconda. https://docs.conda.io/en/latest/miniconda.html

Download consensus_prime.py and primer3_parameters.txt from GitHub.

Adept the primer3_parameters.txt to your needs. For details see https://primer3.org/manual.html

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
conda activate consensus_prime
```
Install pandas
```bash
conda install pandas
```


## Usage

Download the example_data.fna and primer3_parameters.txt then run the following command (remember to adjust the file paths accordingly).
```bash
python3.8 /path_to/consensus_prime.py -infile /path_to/example_data.fna --primer3 /path_to/primer3_parameters.txt
```

**Options:**

command | what it does
  ------------- | -------------
-i, --infile              |Multi-Fasta file with gene sequences.  [required]
-x, --primer3             |Primer3 parameter file. [required]
-o, --outdir              |Specifies output directory. Default = .
-t, --threads             |Number of threads used by MAFFT. Default = -1 (all)
-k, --keepduplicates      |Keep duplicate sequences. Default = False
-c, --consensusthreshold  |Consensus threshold bitween 0 and 1 with 1 beeing a perfect consensus. Default = 0.95
-g, --gapthreshold        |Removes sequences with higher gap to sequence ratio. Default = 0.2
--primers                 |Known primers for visualisation in the final alignment in multifasta format.
--negativesequences       |File with sequences that get their consensus sequence added to the final alignment in multifasta format.
-h, --help                |show this message and exit


Exit the conda environment when you are done
```bash
conda deactivate
```

## We also provide a Docker image for ConsensusPrime
Simply pull and run a ready-to-use image from Dockerhub:

**Show --help**
```bash
docker run mcollatz/consensusprime:1.0
```

**Example**

Download the example_data.fna and primer3_parameters.txt then run the following command (remember to adjust the file paths accordingly).

```bash
docker run -t --rm -v /path/to/dir/with/your/input/files/:/in \
-v /path/to/dir/for/results/:/out \
-u `id -u $USER`:`id -g $USER` \
mcollatz/consensusprime:1.0 \
/consensus_prime.py --infile /in/example_data.fna --primer3 /in/primer3_parameters.txt --outdir /out
```

## Results

The results are stored in the specified /out directory or in the current directory under "/results". The intermediate results of the individual filter steps are also located in this directory. The predicted primers are listed in the "consensus_prime_summary.html" with all related details. The predicted primers can be viewed in "final_alignment.fna" using an alignment visualization program such as ClustalX.

For an example output check the example_results.zip.

