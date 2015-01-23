# CRiSPy
Computing Species Richness in 16S rRNA Pyrosequencing Datasets

## Run
- To run CRiSPy 
```
./crispy.sh <input_fasta_file> [k_value] [band_value]
```
- `k_value`: parameter for kmer distance computation. By default: 6
- `band_value`: parameter of banded global alignment for genetic distance computation. By default: 10 i.e. 1/10 banded alignment


## Require

- Hardware: 
	+ NVIDIA GPU card of compute capability 2.0 or above
- Software:
	+ Linux operating system, preferably Ubuntu
	+ CUDA toolkit installed for code compilation
	+ GCC installed for code compilation
	+ R with R `changepoint` package installed for code execution

## Compile

- To compile:
```
./compile.sh
```

- To compile with for an NVIDIA GPU with compute capability other than 5.0. For example for a Fermi card with compute capability 2.0:
```
./compile.sh 20
```

## License

The project is licensed under the GPL 3.0 license (http://www.gnu.org/licenses/gpl.html).
