#!/bin/bash
##This script will summarize the amino acid and nucleotide results into results folders that can then be processed an analyzed.
#example usage: sh summarize_aa.bash trimmed75
prefix_name=$1

mkdir ${prefix_name}
mkdir ${prefix_name}/aln
mkdir ${prefix_name}/aln_aa
cd scos/aligns/aligns_trim/

for file in $(find . -name '*.nogap.nt.fasta');
do
og=$(basename $file .nogap.nt.fasta)
cp $file ../../../${prefix_name}/aln/${og}.fasta
done

cd ../../../scos_aa/aligns/aligns_trim/
for file in $(find . -name '*.nogap.aa.fasta');
do
og=$(basename $file .nogap.aa.fasta)
cp $file ../../../${prefix_name}/aln_aa/${og}.fasta
done

cd ../../../${prefix_name}/aln
python fasta_to_catnex.py

sh part_converter.sh all_orthologs.nex DNA

cp all_orthologs.phy ../../${prefix_name}.dna.phy
cp all_orthologs.nex ../../${prefix_name}.dna.nex


#sed -i '' 's/\.\///g' COMBINED_hmm_aa.nex_raxml_partition.txt
cp all_orthologs.nex_raxml_partition.nt.txt ../../${prefix_name}.raxml_partition.dna.txt
cp all_orthologs_nexus_charsets.nt.nexus ../../${prefix_name}_nexus_charsets.dna.nexus
rm *.nex
rm *phy
cd ..

cd aln_aa

python fastaaa_to_catnex.py

sh part_converter.sh all_orthologs.nex AA

cp all_orthologs.phy ../../${prefix_name}.aa.phy
cp all_orthologs.nex ../../${prefix_name}.aa.nex


#sed -i '' 's/\.\///g' COMBINED_hmm_aa.nex_raxml_partition.txt
cp all_orthologs.nex_raxml_partition.aa.txt ../../${prefix_name}.raxml_partition.aa.txt
cp all_orthologs_nexus_charsets.aa.nexus ../../${prefix_name}_nexus_charsets.aa.nexus

echo done

rm *.nex
rm *.phy

cd ../..

##convert to codons; use only positions 1 and 2
cat ${prefix_name}.raxml_partition.dna.txt | while read l
do
text="$l"
part1=$(echo $text | awk -F" = " '{print $1}')
part2=$(echo $text | awk -F" = " '{print $2}' )
start=$(echo $part2 | awk -F"-" '{print $1}')
end=$(echo $part2 | awk -F"-" '{print $2}')
three="\3"
pos1=$(printf "${start}-${end}%s" "${three}")
start2=$((start+1))
pos2=$(printf "${start2}-${end}%s" "${three}")
start3=$((start+2))
pos3=$(printf "${start3}-${end}%s" "${three}")
p1=$(echo ${part1}_1 = $pos1)
p2=$(echo ${part1}_2 = $pos2)
p3=$(echo ${part1}_3 = $pos3)
echo $p1 >> ${prefix_name}.12.parts
echo $p2 >> ${prefix_name}.12.parts
echo $p2 >> ${prefix_name}.2.parts
echo $p1 >> ${prefix_name}.123.parts
echo $p2 >> ${prefix_name}.123.parts
echo $p3 >> ${prefix_name}.123.parts
done