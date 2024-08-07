## ribotin

rDNA consensus sequence builder. Input hifi or duplex, and optionally ultralong ONT. Extracts rDNA-specific reads based on k-mer matches to a reference rDNA sequence or based on a [verkko](https://github.com/marbl/verkko) or [hifiasm](https://github.com/chhylp123/hifiasm) assembly, builds a DBG out of them, extracts the most covered path as a consensus and bubbles as variants. Optionally assembles highly abundant rDNA morphs using the ultralong ONT reads.

#### Installation

Installation with [conda](https://docs.conda.io/projects/miniconda/en/latest/) is recommended. Use `conda install ribotin`.

#### Compilation

- `git clone https://github.com/maickrau/ribotin.git`
- `git submodule update --init --recursive`
- `make all`

Also needs [MBG](https://github.com/maickrau/MBG) version 1.0.13 or more recent, [GraphAligner](https://github.com/maickrau/GraphAligner) and [Liftoff](https://github.com/agshumate/Liftoff).

#### Usage

##### Quick start with CHM13 test dataset:

Download the [test dataset](https://zenodo.org/records/10468773/files/ribotin_testdata.tar.gz?download=1). Unzip with `tar -xzf ribotin_testdata.tar.gz` and then run `ribotin-ref -x human -i data/hifi_reads.fa --nano data/ont_reads.fa -o output`. This will run ribotin-ref on the test dataset and store the results in the folder `output`. The assembled morphs will be in `output/morphs.fa` which describes the morph sequences as well as their ONT coverages.

The dataset also has example results with ribotin version 1.2 in `result_v1.2` and instructions for replicating them in `README`.

##### Reference based:

```
bin/ribotin-ref -x human -i hifi_reads1.fa -i hifi_reads2.fq.gz --nano ont_reads.fa -o output_folder
```

This extracts rDNA-specific reads based on k-mer matches to human rDNA, builds a graph and a consensus, and finds variants supported by at least 3 reads. `--nano` is optional, if it is present then ribotin also builds morph consensuses. Results are written to `output_folder`.

##### Verkko based (automatic):

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then run:

```
bin/ribotin-verkko -x human -i /path/to/verkko/assembly -o output_folder_prefix
```

The folder in parameter `-i` should be the same folder as verkko's parameter `-d`. This finds the rDNA tangles based on k-mer matches to human rDNA and assembly graph topology, extracts HiFi reads uniquely assigned to each tangle, and for each tangle builds a graph and a consensus and finds variants supported by at least 3 reads and builds morph consensuses. Results are written per tangle to `output_folder_prefix[x]` where `[x]` is the tangle number.

##### Verkko based (manual):

First you must run a whole genome assembly with [verkko](https://github.com/marbl/verkko). Then manually pick the nodes in each rDNA tangle from `assembly.homopolymer-compressed.noseq.gfa`. Each tangle should only have the rDNA tangle nodes, and not the surrounding nodes. Save the tangle nodes to files with one tangle per file eg `node_tangle1.txt`, `node_tangle2.txt`, `node_tangle3.txt`. Format of the node tangle files should be eg `utig4-1 utig4-2 utig4-3...` or `utig4-1, utig4-2, utig4-3...` or each node in its own line. Then run:

```
bin/ribotin-verkko -x human -i /path/to/verkko/assembly -o output_folder_prefix -c node_tangle1.txt -c node_tangle2.txt -c node_tangle3.txt
```

This extracts HiFi reads uniquely assigned to each node tangle, and for each tangle builds a graph and a consensus and finds variants supported by at least 3 reads and build morph consensuses. Results are written per tangle to `output_folder_prefix[x]` where `[x]` is the tangle number. 

##### Hifiasm based (automatic):

First you must run a whole genome assembly with [hifiasm](https://github.com/chhylp123/hifiasm). Then run:

```
bin/ribotin-hifiasm -x human -a /path/to/hifiasm/assembly_prefix -o output_folder_prefix -i hifi_reads.fa --nano ont_reads.fa
```

The prefix in parameter `-a` should be the same prefix as hifiasm's parameter `-o`, and the reads in `-i` and `--nano` should be the same reads that were given to hifiasm. This finds the rDNA tangles based on k-mer matches to human rDNA and assembly graph topology, extracts HiFi reads uniquely assigned to each tangle, and for each tangle builds a graph and a consensus and finds variants supported by at least 3 reads and builds morph consensuses. Results are written per tangle to `output_folder_prefix[x]` where `[x]` is the tangle number.

##### Hifiasm based (manual):

First you must run a whole genome assembly with [hifiasm](https://github.com/chhylp123/hifiasm). Then manually pick the nodes in each rDNA tangle from `<assembly>.bp.r_utg.noseq.gfa`. Each tangle should only have the rDNA tangle nodes, and not the surrounding nodes. Save the tangle nodes to files with one tangle per file eg `node_tangle1.txt`, `node_tangle2.txt`, `node_tangle3.txt`. Format of the node tangle files should be eg `utig4-1 utig4-2 utig4-3...` or `utig4-1, utig4-2, utig4-3...` or each node in its own line. Then run:

```
bin/ribotin-hifiasm -x human -a /path/to/hifiasm/assembly_prefix -o output_folder_prefix -i hifi_reads.fa --nano ont_reads.fa -c node_tangle1.txt -c node_tangle2.txt -c node_tangle3.txt
```

This extracts HiFi reads uniquely assigned to each node tangle, and for each tangle builds a graph and a consensus and finds variants supported by at least 3 reads and build morph consensuses. Results are written per tangle to `output_folder_prefix[x]` where `[x]` is the tangle number. 

##### Nonhumans

For running `ribotin-ref` on nonhumans replace `-x human` with `--approx-morphsize <morphsize> -r path_to_example_morphs.fa` where `<morphsize>` is the estimated size of a single morph (45000 for human) and `path_to_example_morphs.fa` is a fasta/fastq file which contains (partial or complete) example morphs from the same sample or species. It does not matter if the example morphs are complete or fragments as long as most rDNA k-mers are present. You can get partial reference morphs by doing a whole genome assembly with hifi reads using MBG or a similar hifi based assembly tool, and extracting the sequences of the rDNA tangle from the assembly. 

For `ribotin-verkko` and `ribotin-hifiasm`, replace `-x human` with either `--approx-morphsize <morphsize> --guess-tangles-using-reference path_to_example_morphs.fa` if you have an example morphs file, or with `--approx-morphsize <morphsize> -c tangle1.txt -c tangle2.txt` if you have manually located rDNA tangles from a verkko or hifiasm assembly, where `tangle1.txt tangle2.txt` etc. are files with manually selected rDNA tangle nodes from the verkko or hifiasm assembly with every tangle in a separate file.

If you additionally have one complete morph from the same or related species, you can also include `--orient-by-reference previous_reference_single_morph.fa` to have the results in the same orientation (forward / reverse complement) and offset (rotation) as the previous reference. This file should contain exactly one complete morph.

##### Species with short rDNA morphs

If your species has short rDNA morph size (about 10000bp per morph) you can use the HiFi reads to get very accurate morphs. In this case input the HiFi reads as the nanopore reads as well and adjust the clustering parameters: `-i hifi_reads.fa --nano hifi_reads.fa --morph-cluster-maxedit 10 --morph-recluster-minedit 1 --approx-morphsize <morphsize> -r path_to_reference_kmers.fa`. This will resolve morphs with very small differences.

##### Clustering morphs with ultralong ONT reads

If you have ultralong ONT reads, you can include them to produce consensuses of highly abundant rDNA morphs similar to the CHM13 assembly. For `ribotin-ref` and `ribotin-hifiasm`, add the parameter `--nano /path/to/ont/reads.fa` (multiple files may be added with `--nano file1.fa --nano file2.fa` etc). `ribotin-verkko` will automatically check if ONT reads were used in the assembly and use them, and can be overrode with `--do-ul=no`. This will error correct the ultralong ONT reads by aligning them to the allele graph, extract rDNA morphs from the corrected reads, cluster them based on sequence similarity, and compute a consensus for each cluster. This requires [GraphAligner](https://github.com/maickrau/GraphAligner) to be installed.

##### Annotations

You can lift over annotations with the optional parameters `--annotation-reference-fasta` and `--annotation-gff3`, for example `bin/ribotin-verkko ... --annotation-reference-fasta template_seqs/rDNA_one_unit.fasta --annotation-gff3 template_seqs/rDNA_annotation.gff3`. This requires [liftoff](https://github.com/agshumate/Liftoff) to be installed. The parameter preset `-x human` includes annotations from the KY962518.1 rDNA reference.

#### Output

The output folder will contain several files:

- `nodes.txt`: List of nodes used in this cluster. Only in verkko and hifiasm based modes.
- `hifi_reads.fa`: HiFi or duplex reads used in this cluster.
- `graph.gfa`: de Bruijn graph of the reads.
- `paths.gaf`: Paths of the hifi reads in `graph.gfa`.
- `consensus.fa`: Consensus sequence.
- `consensus_path.gaf`: Path of `consensus.fa` in `graph.gfa`.
- `variants.txt`: A list of variants supported by at least 3 reads. Format is: variant ID, variant path, reference path, variant read support, reference read support, variant sequence, reference sequence.
- `variant-graph.gfa`: `graph.gfa` filtered only to the consensus path and the variant paths in `variants.txt`.
- `variants.vcf`: A list of variants supported by at least 3 reads. Variant IDs match `variants.txt`
- `allele-graph.gfa`: `graph.gfa` filtered only to the consensus path and the variant paths in `variants.txt`, with nodes duplicated so that a path cannot recombine in the middle of a variant.
- `annotation.gff3`: Annotations lifted over from a previous reference. Only if using parameters `--annotation-reference-fasta` and `--annotation-gff3`

The following files are created when ultralong ONT reads are included:

- `ont-alns.gaf`: Aligned paths of ultralong ONT reads to `allele-graph.gfa`.
- `loops.fa`: A list of individual rDNA morphs found in the ultralong ONT reads.
- `morphs.fa`: A list of rDNA morph consensuses and their abundances.
- `morphs.gaf`: The paths of the rDNA morph consensuses in `allele-graph.gfa`.
- `morphgraph.gfa`: A graph describing how the morph consensuses connect to each others.
- `readpaths-morphgraph.gaf`: Paths of the ONT reads in `morphgraph.gfa`. Only shows reads which are assigned to complete morphs.
