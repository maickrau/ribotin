## rdnaConsensus

rDNA consensus sequence builder. Input hifi. Extracts rDNA-specific reads based on k-mer matches to a reference rDNA sequence or based on a [verkko](https://github.com/marbl/verkko) assembly, builds a DBG out of them, extracts the most covered path as a consensus and bubbles as variants.

#### Compilation

- `git clone https://github.com/maickrau/rdnaConsensus.git`
- `git submodule update --init --recursive`
- `make all`

Also needs [MBG](https://github.com/maickrau/MBG) version 1.0.13 or more recent.

#### Usage

###### Reference based:
```
bin/rdnaConsensus-ref -r reference.fa -o output_folder --mbg /path/to/MBG -i hifi_reads1.fa -i hifi_reads2.fq.gz
```

This extracts rDNA-specific reads based on k-mer matches to `reference.fa`, builds a graph and a consensus, and finds variants supported by at least 3 reads. Results are written to `output_folder`.

###### Verkko based:

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then pick the nodes in each rDNA cluster manually, and save them to files eg `node_cluster1.txt`, `node_cluster2.txt`, `node_cluster3.txt`. Then run:

```
bin/rdnaConsensus-verkko -i /path/to/verkko/assembly --mbg /path/to/MBG -o output_folder_prefix -c node_cluster1.txt -c node_cluster2.txt -c node_cluster3.txt
```

This extracts HiFi reads uniquely assigned to each node cluster, and for each cluster builds a graph and a consensus and finds variants supported by at least 3 reads. Results are written per cluster to `output_folder_prefix[x]` where `[x]` is the cluster number.

#### Output

The output folder will contain several files:

- `nodes.txt`: List of nodes used in this cluster. Only in verkko based mode. Same as the input node cluster.
- `reads.fa`: Reads used in this cluster.
- `graph.gfa`: de Bruijn graph of the reads.
- `paths.gaf`: Paths of the hifi reads in `graph.gfa`.
- `consensus.fa`: Consensus sequence of the rDNA.
- `consensus_path.gaf`: Path of the consensus sequence in `graph.gfa`.
- `variants.txt`: A list of variants supported by at least 3 reads. Format is: variant path, reference path, variant read support, reference read support.
- `variant-graph.gfa`: The graph filtered only to the consensus path and variants.

#### Todo

- flip and rotate the consensus to match reference orientation
- useful variant output
