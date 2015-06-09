# AssemblyComparison

AssemblyComparison provides insight into the differences of genome assemblers Velvet and Minia. By starting with a contig file, AssemblyComparison will use Velvet and Minia to assemble a genome then use bowtie2 to compare the differences between the coverages. More metrics will be added soon.

# Required Software and Libraries

In order to use AssemblyComparison, the following software and libraries must be installed and added to your PATH:

- Python 2.7
- Velvet
- Minia
- bowtie2
- pySAM
- numpy
- mathplotlib.py
