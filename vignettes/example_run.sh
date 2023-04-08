python3 HapChIP/HapChIP.py -b example_input/example.bam -v example_input/example.vcf.gz -o ./example_output/
samtools index example_output/example.bam.hap1.bam
samtools index example_output/example.bam.hap2.bam
python3 HapChIP/cistrans_parse.py -b example_output/example.bam.hap1.bam -o example_output/ --name hap1
python3 HapChIP/cistrans_parse.py -b example_output/example.bam.hap2.bam -o example_output/ --name hap2