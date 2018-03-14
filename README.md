# Tandem Repeat Finder
Finds exact repeats in a fasta file.

## Usage
python main.py --repeat 2 --min 16 < fastafile 2>log >result.csv

## Command line parameters
* *repeat*: length of short repeat. e.g. 1 will find mononucleotide repeats, 2 will find dinucleotide repeats
* *min*: only report repeats of at least this length

## Output
CSV file with the following fields:
* chrom: name of fasta sequence
* pos: start position of tandem repeat
* length: total length of tandem repeat
* kmer: repeated kmer

