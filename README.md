## rdnaConsensus

rDNA consensus sequence builder. Input hifi. Extracts rDNA-specific reads, builds a DBG out of them, extracts the most covered path as a consensus and bubbles as variants.

#### Compilation

- `git clone https://github.com/maickrau/rdnaConsensus.git`
- `git submodule update --init --recursive`
- `make all`

Also needs [MBG](https://github.com/maickrau/MBG) commit [bf5a22d](https://github.com/maickrau/MBG/commit/bf5a22dc9914e752dc807384e99d4b3e9c7f21a0) or more recent.

#### Usage

```
bin/seqPicker 201 2000 template_seqs/chm13_rDNAs.fa input_hifi_read_file.fa > rDNA_reads.fa
MBG -i rDNA_reads.fa -o graph.gfa -k 101 -w 70 -a 2 -u 3 -r 15000 -R 4000 --error-masking=msat --output-sequence-paths paths.gaf --only-local-resolve
scripts/get_heaviest_path_acyclic.py heavy_path < graph.gfa > consensus_seq.fa
scripts/get_variants.py graph.gfa paths.gaf consensus_seq.fa 10 > variants.txt
```

#### Todo

- flip and rotate the consensus to match reference orientation
- useful variant output
