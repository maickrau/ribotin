## ribotin

rDNA consensus sequence builder. Input hifi or duplex. Extracts rDNA-specific reads based on k-mer matches to a reference rDNA sequence or based on a [verkko](https://github.com/marbl/verkko) assembly, builds a DBG out of them, extracts the most covered path as a consensus and bubbles as variants.

#### Compilation

- `git clone https://github.com/maickrau/ribotin.git`
- `git submodule update --init --recursive`
- `make all`

Also needs [MBG](https://github.com/maickrau/MBG) version 1.0.13 or more recent.

#### Usage

##### Reference based:

```
bin/ribotin-ref -r template_seqs/chm13_rDNAs.fa -o output_folder --mbg /path/to/MBG -i hifi_reads1.fa -i hifi_reads2.fq.gz --orient-by-reference template_seqs/rDNA_one_unit.fasta
```

This extracts rDNA-specific reads based on k-mer matches to `template_seqs/chm13_rDNAs.fa`, builds a graph and a consensus, and finds variants supported by at least 3 reads. Results are written to `output_folder`.

##### Verkko based (automatic):

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then run:

```
bin/ribotin-verkko -i /path/to/verkko/assembly --mbg /path/to/MBG -o output_folder_prefix --guess-clusters-using-reference template_seqs/chm13_rDNAs.fa --orient-by-reference template_seqs/rDNA_one_unit.fasta
```

This finds the rDNA clusters based on k-mer matches and assembly graph topology, extracts HiFi reads uniquely assigned to each cluster, and for each cluster builds a graph and a consensus and finds variants supported by at least 3 reads. Results are written per cluster to `output_folder_prefix[x]` where `[x]` is the cluster number. You can also add the extra parameter `--guess-clusters-using-reference template_seqs/chm13_mito.fa` to detect the mitochondria and create a consensus sequence of it.

##### Verkko based (manual):

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then manually pick the nodes in each rDNA cluster from `assembly.homopolymer-compressed.noseq.gfa`, and save them to files with one cluster per file eg `node_cluster1.txt`, `node_cluster2.txt`, `node_cluster3.txt`. Format of the node cluster files should be eg `utig4-1 utig4-2 utig4-3...` or `utig4-1, utig4-2, utig4-3...` or each node in its own line. Then run:

```
bin/ribotin-verkko -i /path/to/verkko/assembly --mbg /path/to/MBG -o output_folder_prefix -c node_cluster1.txt -c node_cluster2.txt -c node_cluster3.txt --orient-by-reference template_seqs/rDNA_one_unit.fasta
```

This extracts HiFi reads uniquely assigned to each node cluster, and for each cluster builds a graph and a consensus and finds variants supported by at least 3 reads. Results are written per cluster to `output_folder_prefix[x]` where `[x]` is the cluster number.

##### Annotations

You can lift over annotations with the optional parameters `--annotation-reference-fasta` and `--annotation-gff3`, for example `bin/ribotin-verkko ... --annotation-reference-fasta template_seqs/rDNA_one_unit.fasta --annotation-gff3 template_seqs/rDNA_annotation.gff3`. This requires [liftoff](https://github.com/agshumate/Liftoff) to be installed.

#### Output

The output folder will contain several files:

- `nodes.txt`: List of nodes used in this cluster. Only in verkko based mode.
- `hifi_reads.fa`: HiFi or duplex reads used in this cluster.
- `graph.gfa`: de Bruijn graph of the reads.
- `paths.gaf`: Paths of the hifi reads in `graph.gfa`.
- `consensus.fa`: Consensus sequence.
- `consensus_path.gaf`: Path of `consensus.fa` in `graph.gfa`.
- `variants.txt`: A list of variants supported by at least 3 reads. Format is: variant ID, variant path, reference path, variant read support, reference read support, variant sequence, reference sequence.
- `variant-graph.gfa`: `graph.gfa` filtered only to the consensus path and the variant paths in `variants.txt`.
- `variants.vcf`: A list of variants supported by at least 3 reads. Variant IDs match `variants.txt`
- `annotation.gff3`: Annotations lifted over from a previous reference. Only if using parameters `--annotation-reference-fasta` and `--annotation-gff3`
