# EpiPred
Wrapper for the [IEDB](http://www.iedb.org/) MHC class II epitope predictor.

#Installation
\- Download and configure the IEDB MHC class II epitope predictor from `http://tools.iedb.org/mhcii/download/`
\- Checkout this source: `git clone https://github.com/kamanufred/EpiPred.git`  

#Usage
Run `python EpiPred.py -h` to see all the options.  

Usage: EpiPred.py [options]

    Options:
    -h, --help            show this help message and exit
    -p PEPTIDEPATH, --peptides=PEPTIDEPATH
                        A path to the directory holding the peptides
    -a ALLELEFILE, --allele=ALLELEFILE
                        File with list of alleles
    -m METHOD, --method=METHOD
                    Prediction method   
    -o REPORTFILE, --output=REPORTFILE
                        The final output file
