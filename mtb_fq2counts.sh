#Baljai Sundararaman
#11-14-22

#!/bin/bash

mkdir $1
cd $1

if [ ! -f ./$1.sub.cutad.R1.fq.gz ]
then
    seqtk sample -s100 $2 3000000 | gzip > ./$1.sub.R1.fq.gz
    seqtk sample -s100 $3 3000000 | gzip > ./$1.sub.R2.fq.gz
    wait
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --pair-filter=any -m25 -o ./$1.sub.cutad.R1.fq -p ./$1.sub.cutad.R2.fq ./$1.sub.R1.fq.gz ./$1.sub.R2.fq.gz > ./$1.sub.cutad.log
    wait
    gzip ./*.fq
fi

if [ ! -f ./$1.sub.rd.s.q20.bam ]
then
    bwa mem -t16 /path/to/Mtb_reference/GCF_000195955.2_ASM19595v2_genomic.fna.gz ./$1.sub.cutad.R1.fq.gz ./$1.sub.cutad.R2.fq.gz | samtools view -Sbh - | samtools sort - -o ./$1.sub.s.bam
    wait
    samtools index ./$1.sub.s.bam
    wait
    samtools rmdup -Ss ./$1.sub.s.bam ./$1.sub.rd.bam
    wait
    samtools sort ./$1.sub.rd.bam -o ./$1.sub.rd.s.bam
    wait
    samtools index ./$1.sub.rd.s.bam
    wait
    samtools view -q20 ./$1.sub.rd.s.bam -o ./$1.sub.rd.s.q20.bam
    wait
    #samtools stats ./$1.sub.rd.s.bam | head -50 > ./$1.sub.sam.stats1.txt
    #samtools stats ./$1.sub.rd.s.q20.bam | head -50 > ./$1.sub.sam.stats2.txt
    #wait
fi

wait
if [ ! -f ./$1.sub.stats.csv ]
then
    samtools stats ./$1.sub.s.bam | head -50 > ./$1.sub.sam.stats1.txt
    wait
    samtools stats ./$1.sub.rd.s.bam | head -50 > ./$1.sub.sam.stats2.txt
    wait
    samtools stats ./$1.sub.rd.s.q20.bam | head -50 > ./$1.sub.sam.stats3.txt
    wait
    echo -e 'Sample\n'$1'' > ./$1.sub.stats.csv
    sed -i -e "1s/$/,Raw_readPairs/; 2s/$/,$(zcat $2 | wc -l | awk '{print$1/4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Trimmed_readPairs/; 2s/$/,$(zcat ./$1.sub.cutad.R1.fq.gz | wc -l | awk '{print$1/4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_Reads/; 2s/$/,$(grep '^SN.seq' ./$1.sub.sam.stats1.txt | awk '{print$3}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_Bases/; 2s/$/,$(grep '^SN.total.length' ./$1.sub.sam.stats1.txt | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_UQ_reads/; 2s/$/,$(grep '^SN.seq' ./$1.sub.sam.stats2.txt | awk '{print$3}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_UQ_bases/; 2s/$/,$(grep '^SN.total.length' ./$1.sub.sam.stats2.txt | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_Mq20_reads/; 2s/$/,$(grep '^SN.seq' ./$1.sub.sam.stats3.txt | awk '{print$3}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,PF_Mq20_bases/; 2s/$/,$(grep '^SN.total.length' ./$1.sub.sam.stats3.txt | awk '{print$4}')/" ./$1.sub.stats.csv

    sed -i -e "1s/$/,All_mappedReads/; 2s/$/,$(grep '^SN.reads.mapped:' ./$1.sub.sam.stats1.txt | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,All_mappedBases/; 2s/$/,$(grep '^SN.bases.mapped.(cigar):' ./$1.sub.sam.stats1.txt | awk '{print$5}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,UQ_mappedReads/; 2s/$/,$(grep '^SN.reads.mapped:' ./$1.sub.sam.stats2.txt | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,UQ_mappedBases/; 2s/$/,$(grep '^SN.bases.mapped.(cigar):' ./$1.sub.sam.stats2.txt | awk '{print$5}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Mq20_Reads/; 2s/$/,$(grep '^SN.reads.mapped:' ./$1.sub.sam.stats3.txt | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Mq20_Bases/; 2s/$/,$(grep '^SN.bases.mapped.(cigar):' ./$1.sub.sam.stats3.txt | awk '{print$5}')/" ./$1.sub.stats.csv

    sed -i -e "1s/$/,Genome_Cov/; 2s/$/,$(bedtools coverage -a /redser4/projects/bsundara/refs/mtb/Mtb_H37Rv_genome.bed -b ./$1.sub.s.bam -mean | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_Cov/; 2s/$/,$(bedtools coverage -a /redser4/projects/bsundara/refs/mtb/Mtb_H37Rv_genome.bed -b ./$1.sub.rd.s.bam -mean | awk '{print$4}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_Mq20_Cov/; 2s/$/,$(bedtools coverage -a /redser4/projects/bsundara/refs/mtb/Mtb_H37Rv_genome.bed -b ./$1.sub.rd.s.q20.bam -mean | awk '{print$4}')/" ./$1.sub.stats.csv

    sed -i -e "1s/$/,Genome_all_1X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 1 | grep 'genome' | awk 'END {if($2==1) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_all_10X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 10 | grep 'genome' | awk 'END {if($2==10) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_all_25X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 25 | grep 'genome' | awk 'END {if($2==25) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_all_50X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 50 | grep 'genome' | awk 'END {if($2==50) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_all_100X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 100 | grep 'genome' | awk 'END {if($2==100) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_all_125X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.s.bam -max 125 | grep 'genome' | awk 'END {if($2==125) print$5; else print"0"}')/" ./$1.sub.stats.csv

    sed -i -e "1s/$/,Genome_UQ_1X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 1 | grep 'genome' | awk 'END {if($2==1) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_10X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 10 | grep 'genome' | awk 'END {if($2==10) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_25X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 25 | grep 'genome' | awk 'END {if($2==25) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_50X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 50 | grep 'genome' | awk 'END {if($2==50) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_100X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 100 | grep 'genome' | awk 'END {if($2==100) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_UQ_125X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.bam -max 125 | grep 'genome' | awk 'END {if($2==125) print$5; else print"0"}')/" ./$1.sub.stats.csv

    sed -i -e "1s/$/,Genome_q20_1X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 1 | grep 'genome' | awk 'END {if($2==1) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_q20_10X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 10 | grep 'genome' | awk 'END {if($2==10) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_q20_25X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 25 | grep 'genome' | awk 'END {if($2==25) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_q20_50X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 50 | grep 'genome' | awk 'END {if($2==50) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_q20_100X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 100 | grep 'genome' | awk 'END {if($2==100) print$5; else print"0"}')/" ./$1.sub.stats.csv
    sed -i -e "1s/$/,Genome_q20_125X/; 2s/$/,$(bedtools genomecov -ibam ./$1.sub.rd.s.q20.bam -max 125 | grep 'genome' | awk 'END {if($2==125) print$5; else print"0"}')/" ./$1.sub.stats.csv
fi

#bedtools multicov -bams ./$1.sub.rd.s.bam -bed /redser4/projects/bsundara/refs/mtb/Mtb_H37Rv_100bb_bins.bed > ./$1.sub.cov_100bp.tsv
#wait
#rm ./$1.sub.rd.bam
#rm ./$1.sub.sam.stats*.txt
#rm ./$1.sub.s.bam
#echo "Done counting $1.sub.for MTb capture!"
