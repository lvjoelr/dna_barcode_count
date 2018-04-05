# dna_barcode_count
Counting multiple DNA barcodes in gzipped fastq 

The python script intends to count sequencing reads with perfect matching with pre-designed DNA barcodes at specific locations.

Usage:
python3 barcode_count.py -b barcode_list.xlsx -f first1M.fastq.gz -o first1M_barcode_recognization_result -p 4

Dependencies:
numpy
pandas
openpyxl
matplotlib
multiprocessing

