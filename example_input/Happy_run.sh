pip install numpy
pip install pysam

python3 ../HapChIP/HapChIP.py -b example.bam -v example.vcf.gz -o ../example_output

samtools index ../example_output/example.bam.hap0.bam
samtools index ../example_output/example.bam.hap1.bam
samtools index ../example_output/example.bam.hap2.bam
samtools index ../example_output/example.bam.hap3.bam
