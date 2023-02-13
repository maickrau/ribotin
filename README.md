## rdnaConsensus

rDNA consensus sequence builder. Input hifi. Extracts rDNA-specific reads based on k-mer matches to a reference rDNA sequence or based on a [verkko](https://github.com/marbl/verkko) assembly, builds a DBG out of them, extracts the most covered path as a consensus and bubbles as variants.

#### Compilation

- `git clone https://github.com/maickrau/rdnaConsensus.git`
- `git submodule update --init --recursive`
- `make all`

Also needs [MBG](https://github.com/maickrau/MBG) version 1.0.13 or more recent.

#### Usage

Reference based:
```
bin/rdnaConsensus-ref -r reference.fa -o output_folder --mbg /path/to/MBG -i hifi_reads1.fa -i hifi_reads2.fq.gz
```

This extracts rDNA-specific reads based on k-mer matches to `reference.fa`, builds a graph and a consensus, and finds variants supported by at least 3 reads. Results are written to `output_folder`.

Verkko based:

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then pick the nodes in each rDNA cluster manually, and save them to files eg `node_cluster1.txt`, `node_cluster2.txt`, `node_cluster3.txt`. Then run:

```
bin/rdnaConsensus-verkko -i /path/to/verkko/assembly -o output_folder_prefix -c node_cluster1.txt -c node_cluster2.txt -c node_cluster3.txt
```

This extracts HiFi reads uniquely assigned to each node cluster, and for each cluster builds a graph and a consensus and finds variants supported by at least 3 reads. Results are written per cluster to `output_folder_prefix[x]` where `[x]` is the cluster number.

#### Todo

- flip and rotate the consensus to match reference orientation
- useful variant output
