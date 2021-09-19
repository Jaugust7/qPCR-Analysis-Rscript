# qPCR Analysis Rscript
 a simple command line tool for qPCR analysis

**Usage**

### First run the script with the help file to load the requisite packages


```console
Rscript <path>/qPCR.R -h

-f .XLSX, --file=.XLSX
  dataset file name

-o .CSV, --out=.CSV
  output file name for stats

-c CONTROL, --control=CONTROL
  name of control sample to be used

-p .PDF, --pdf=.PDF
  output file name for PDF

-w INTEGER, --width=INTEGER
  width of PDF

-t INTEGER, --height=INTEGER
  height of PDF

-g GENE, --gene=GENE
  gene to be used for comparison

-h, --help
  Show this help message and exit
```
