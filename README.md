Workflow for Tri-C analysis:

1. Removal of adapter sequences using Trim Galore!.
http://www. bioinformatics.babraham.ac.uk/projects/trim_galore/

2. Reconstruction of overlapping reads using FLASH.
Magoč, T., and Salzberg, S.L. (2011). FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics 27, 2957–2963.

3. In silicon restriction enzyme digestion using custom script.
fastq_digester.pl

4. Alignment of digested reads to reference genome using Bowtie.
Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 10, R25.

5. Extraction of unique reads containing restriction fragments interacting with viewpoints of interest using custom script.
tric_reporters.pl

6. Calculation of interaction frequencies for each reported interaction between fragments within a ≥3plet using custom script.
tric_frequencies.pl
