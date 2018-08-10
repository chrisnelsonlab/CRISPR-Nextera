# CRISPR-Nextera
Code for characterizing complex genomic alterations caused by CRISPR gene editing using Illumina Nextera. 


# Parameter file
Line 1 - Amplicon sequence if no editing takes place

Line 2 - Amplicon sequence if a targeted deletion occurs

Line 3 - Amplicon sequence if a targeted inversion occurs

Line 4 - DSB location in bp

Line 5 - BP window to check for indels on each side of DSB

Line 6 - Minimum score for aligning 

Line 7 - Minimum score for aligning short portion

Line 8 - Minimum alignment score for the AAV genome

Line 9 - Large portion of target gene with inversion for unexpected modifications (e.g. chewback)

Line 10 -  Large portion of target gene for unexpected modifications (e.g. chewback)

Line 11 - minimum amplicon size. Reject shorter amplicons

Line 12 - AAV Genome 1

Line 13 - AAV genome 2

Line 14 - Value for gap opening when aligning to the dystrophin gene

Line 15 - Value for gap opening when aligning to the AAV genome

Line 16 - Value for gap extending for AAV genome
