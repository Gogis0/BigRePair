# Big-Repair

*bigrepair* is a grammar compressor for huge files many repetitions. *bigrepair* uses Context Triggered Piecewise Hashing on the input file to parse it into phrases which are later proceesed by RePair. See [1] for further details and experimental results. 

Copyrights 2019- by the authors. 
 

## Installation

* Download/Clone the repository
* `make` (create the C/C++ executables) 
* `bigrepair -h` (get usage instruction)

Note that `bigrepair` is a Python script so you need at least Python 3.4 installed.
 


## Sample usage

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02 which are used internally by the CTPH parsing. To build a grammar for file *yeast.fasta* just type

       bigrepair yeast.fasta

If no errors occur the files yeast.fasta.C and yeast.fasta.R are created.

To recover the original file, type

       bigrepair -d yeast.fasta

this command will read the yeast.fasta.C and yeast.fasta.R files and produce a yeast.fasta.out file identical to the original input yeast.fasta. 

The CTPH parsing step has limited support for multiple threads. Use `bigrepair` option `-t` to specify the number of helper threads: in our tests `-t 8` reduced the running time of the parsing by roughly a factor 6.

For very large input files (or not so large but without many repetitions), RePair may run out of memory and crash. In this case use option `-m` to limit RePair RAM. usage

## References

\[1\] Travis Gagie, Tomohiro I, Giovanni Manzini, Gonzalo Navarro, Hiroshi Sakamoto, Yoshimasa Takabatake: *Rpair: Rescaling RePair with Rsync*. [CoRR abs/1906.00809](https://arxiv.org/abs/1906.00809) (2019)
 
