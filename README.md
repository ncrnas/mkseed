mkseed
======

`mkseed` is a simple tool that aims to find microRNA seed sites.

Algorithms
----------

`mkseed` has been developed with C++ and [SeqAn](https://www.seqan.de). `mkseed` uses several algorithms offered by [SeqAn](https://www.seqan.de) to find microRNA target sites.

-   Suffix array
-   N-mer
-   Horspool's algorithm
-   Brute-force

Installation
------------

`mkseed` bundles all dependent header libraries and uses [CMake](https://cmake.org/) to compile the source code.

Programs
--------

The current version of `mkseed` only provides a single executable.

| Program              | Description                                                |
|:---------------------|:-----------------------------------------------------------|
| find_seed            | Find microRNA seed sites in mRNA sequences                 |


`find_seed` requires three command line arguments to process seed finding. 

| Input                 | Description                                                |
|:--------------------- |:-----------------------------------------------------------|
| miRNA fasta file name | micorRNA sequences in Fasta format                         |
| mRNA fasta file name  | messengerRNA sequences in Fasta foramt                     |
| Output file name      | Output results will be written in this file                |

Usage
-----

```
find_seed /path/to/mirna.fasta /path/to/mrna.fasta /path/to/output.txt
```
