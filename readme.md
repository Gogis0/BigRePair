# Big-Repair

*bigrepair* is a grammar compressor for huge files with many repetitions. *bigrepair* uses Context Triggered Piecewise Hashing on the input file to parse it into phrases which are later processed by RePair. See [1] for further details and experimental results. 

Copyrights 2019- by the authors. 
 

## Installation

* Download/Clone the repository
* `make` (create the C/C++ executables) 
* `bigrepair -h` (get usage instruction)

Note that `bigrepair` is a Python script so you need at least Python 3.4 installed.
 


## Sample usage

To build a grammar for file *yeast.fasta* just type

       bigrepair yeast.fasta

If no errors occur the files yeast.fasta.C and yeast.fasta.R are created.

To recover the original file, type

       bigrepair -d yeast.fasta

this command will read the yeast.fasta.C and yeast.fasta.R files and produce a yeast.fasta.out file identical to the original input yeast.fasta. 

The CTPH parsing step has limited support for multiple threads. Use `bigrepair` option `-t` to specify the number of helper threads: in our tests `-t 8` reduced the running time of the parsing step by roughly a factor 6. 

For very large input files (or not so large but without many repetitions), RePair may run out of memory and crash. In this case use option `-m` to limit RePair RAM usage
 

### Experimental feature (April 2020)

Using the option `-i` *bigrepair* will consider the input file as a sequence of 32bit integers (hence the input file size must be a multiple of 4). The output grammar will have integers as non-terminal. To get the original file use again *bigrepair* with option `-di`. Currently all the integers in the input file must be smaller than 2<sup>30</sup>.


## References

\[1\] Travis Gagie, Tomohiro I, Giovanni Manzini, Gonzalo Navarro, Hiroshi Sakamoto, Yoshimasa Takabatake: *Rpair: Rescaling RePair with Rsync*. [Proc. SPIRE '19](https://link.springer.com/chapter/10.1007%2F978-3-030-32686-9_3) [CoRR abs/1906.00809](https://arxiv.org/abs/1906.00809) (2019)

