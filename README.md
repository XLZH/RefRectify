# RefRectify
Rectify the reference (germline Snp, Insert and Deletion) with the result of mpileup

# Background
Some of the SNP, Insert and Deletion have an alt-frequence more than 100%, we call this kind of mutation "Germline mutation"
This program is used to rectify the germline mutation with the input of mpileup

# Usage
## compile
gcc -std=c99 -o RefRectify RefRectify.3.c
## comand line
        RefRectify [options]
        -h        helpinfo
        -i        inputfile  [mpileup file]
        -o        outputfile [rectified reference]
        -r        original reference that will be rectified
        -c        cutoff for alt-frequence [default = 0.9]
