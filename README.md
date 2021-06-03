# ConsensusPrime
Pipeline to identify ideal consensus regions for primer design.

## Dependencies
python3.8

pandas

mafft

primer3

clustalx (or another alignment visualisation programm)


## We also provide a Docker image for ConsensusPrime
Simply pull and run a ready-to-use image from Dockerhub:  
```bash
docker run -t --rm -v /path/to/dir/with/your/input/files/:/in \
-v /path/to/dir/for/results/:/out \
-u `id -u $USER`:`id -g $USER` \
mcollatz/consensusprime:1.0 \
/consensus_prime.py --infile /in/multifasta.fas --primer3 /in/primer3_parameters.txt --outdir /out
```
