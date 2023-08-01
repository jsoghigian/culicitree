#!/bin/bash
#Convert an input partition file to raxml style and outputs both raxml style and nexus style.

data_type=$2
if [[ $data_type == "AA" ]];
then
	file=$(echo $1 | sed 's/.nex//g')
	sed '/begin sets/,/charpartition/!d;//d' $1 | sed 's/charset/LG,/g' | sed 's/;//' | sed 's/\.\/scos_aa\/aligns\/aligns_trim\///g' | sed 's/\.nex//' | sed 's/\.\///' > ${1}_raxml_partition.aa.txt
	echo \#NEXUS > tmp.nex
	awk '/^begin sets;$/,0' all_orthologs.nex >> tmp.nex
	sed 's/\.\///g' tmp.nex >> tmp2.nex
	cat tmp2.nex > ${file}_nexus_charsets.aa.nexus
	rm tmp.nex
	rm tmp2.nex
else
	file=$(echo $1 | sed 's/.nex//g')
	sed '/begin sets/,/charpartition/!d;//d' $1 | sed 's/charset/DNA,/g' | sed 's/;//' | sed 's/\.\/scos_aa\/aligns\/aligns_trim\///g' | sed 's/\.nex//' | sed 's/\.\///' > ${1}_raxml_partition.nt.txt
	echo \#NEXUS > ${file}_nexus_charsets.nt.nexus
	awk '/^begin sets;$/,0' all_orthologs.nex >> tmp.nex
	sed 's/\.\///g' tmp.nex >> tmp2.nex
	cat tmp2.nex > ${file}_nexus_charsets.nt.nexus
	rm tmp.nex
	rm tmp2.nex
fi