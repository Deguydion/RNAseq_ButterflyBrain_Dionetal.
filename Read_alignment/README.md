** Purpose:
Aligns reads to the genome, creates an assembly based on reads and generates abundance tables (gene and transcript counts) to be used in DEG, DSG analysis and annotation. 
Program used are hisat2 (http://daehwankimlab.github.io/hisat2/) and stringtie (https://ccb.jhu.edu/software/stringtie/)

**Steps:
1. Create .bam files for each library using hisat2.loop.sh
2. Create .gtf files from each .bam files using stringtie_gtf_to_bam.sh
3. Follow Assembly_building to get an assembly made from the libraries. We alos alos create extra assembly files that contain short gene ids, and only one isoform (the longest) for annotation purposes). We also use assembly_stat to get basic stats on the assembly.
4. Create abundance tables (gene and transcript count tables) using stringtie_abundance.sh and abundance.tables

