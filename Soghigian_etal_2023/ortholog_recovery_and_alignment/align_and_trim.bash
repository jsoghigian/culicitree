#!/bin/bash
#Input var 1 should be the path to the ortholog set from orthograph; alternatively, a directory structure in which the aligned reference sequences are stored in the subfolder aln, and formatted as ortholog_name.aln.fa 
#Input var2 should be the set of directories to use as input.  Each directory is expected to contain an aa_summarized and nt_summarized folder.  These directories must be in a subfolder called input. E.g., if you want to align the sample cmmex, your var2 should list cmmex , which would also have a corresponding directory in input/cmmex
#Input var3 should be the set of orthologs to align from the set as a text file.
#This script epxects a python script, fasta_subset.py, to be available; change path below, as needed.
#example:
#sh align_and_trim.bash /work/soghigian_lab/databases/culicologs/sets/culicologs samples_to_align.txt threshold_scos.txt
mkdir scos_aa
mkdir scos_aa/folders
mkdir scos_aa/aligns
mkdir scos_aa/aligns/aligns_raw
mkdir scos_aa/aligns/aligns_trim
mkdir scos
mkdir scos/folders
mkdir scos/aligns
mkdir scos/aligns/aligns_trim
mkdir logs
i=1
tot=$(awk 'END{print NR}' $3 )
echo Aligning orthologs from $3 ... Total of $tot orthologs.
for usco in `cat $3`
    do
    mkdir scos_aa/folders/${usco}
    mkdir scos/folders/${usco}
    for sample in `cat $2`
        do
        base=$(basename $sample)
        if [ -f "input/${sample}/aa_summarized/${usco}.aa.summarized.fa" ]; then
               echo ">"${base} >> scos_aa/folders/${usco}/${base}_${usco} 
               tail -n 1 input/${sample}/aa_summarized/${usco}.aa.summarized.fa >> scos_aa/folders/${usco}/${base}_${usco} 2> /dev/null
               echo ">"${base} >> scos/folders/${usco}/${base}_${usco}
               tail -n 1 input/${sample}/nt_summarized/${usco}.nt.summarized.fa >> scos/folders/${usco}/${base}_${usco} 2> /dev/null
        else 
            echo "${sample} does not have ${usco}" >> logs/${usco}_align_log.txt
        fi
    done
#join together AA and NT files per ortholog so they are singles
    cat scos_aa/folders/${usco}/* > scos_aa/folders/${usco}/${usco}.fa
    cat scos/folders/${usco}/* > scos/folders/${usco}/${usco}.fa
#align amino acids with mafft --addfragments
    mafft --globalpair --addfragments scos_aa/folders/${usco}/${usco}.fa --maxiterate 1000 --6merpair --thread 8 ${1}/aln/${usco}.aln.fa > scos_aa/aligns/aligns_raw/${usco}.full.fa 2>> logs/${usco}_align_log.txt
    cp scos_aa/aligns/aligns_raw/${usco}.full.fa scos_aa/aligns/aligns_raw/${usco}.fix.fa
    rm scos_aa/aligns/aligns_raw/${usco}.full.fa
#retrieve reference IDs
    grep ">" ${1}/aln/${usco}.aln.fa | sed 's/>//g' > scos_aa/aligns/aligns_raw/ref_${usco}.ids.txt
#remove the references from the amino acid set, so it now has the same species as the NT set. This is used to report gaps to the NT set, then gets removed.  You could leave the reference set in, too, so long as you didn't pass the reference set through the pipeline as well.
    python fasta_subset.py scos_aa/aligns/aligns_raw/${usco}.fix.fa scos_aa/aligns/aligns_raw/ref_${usco}.ids.txt > scos_aa/aligns/aligns_raw/${usco}.full.noref.fa
#let's remove Ns and Xs from alignments:

#We'll use gappyout here to prune the alignments.  An alternative is to use something like a fixed threshold to remove columns with heavy missingness... Or not to trim at all.
    rm scos_aa/aligns/aligns_raw/${usco}.fix.fa
    trimal -gappyout -in scos_aa/aligns/aligns_raw/${usco}.full.noref.fa -backtrans scos/folders/${usco}/${usco}.fa -out scos/aligns/aligns_trim/${usco}.fix.nt.trimal.fasta  -ignorestopcodon 2>> logs/${usco}_align_log.txt
    trimal -gappyout -in scos_aa/aligns/aligns_raw/${usco}.full.noref.fa -out scos_aa/aligns/aligns_trim/${usco}.trimal.fasta 2>>  logs/${usco}_align_log.txt
#last step is just to remove the Ns, Xs, and columns that are only gaps (e.g. missing due to the reference alignment - probably removed during gappyout, but because of the Ns/Xs, may not have been).
    sed -e '/^[^>]/s/X/-/g' scos_aa/aligns/aligns_trim/${usco}.trimal.fasta > scos_aa/aligns/aligns_trim/${usco}.xtogap.aa.fasta
    sed -e '/^[^>]/s/N/-/g' scos/aligns/aligns_trim/${usco}.fix.nt.trimal.fasta > scos/aligns/aligns_trim/${usco}.ntogap.nt.fasta
    trimal -gappyout -in scos_aa/aligns/aligns_trim/${usco}.xtogap.aa.fasta -out scos_aa/aligns/aligns_trim/${usco}.nogap.aa.fasta 2>>  logs/${usco}_align_log.txt
    trimal -gappyout -in scos/aligns/aligns_trim/${usco}.ntogap.nt.fasta -out scos/aligns/aligns_trim/${usco}.nogap.nt.fasta 2>>  logs/${usco}_align_log.txt
    rm -r scos/folders/${usco}
    rm -r scos_aa/folders/${usco}
    echo $i / $tot : ${usco} complete.
    i=$((i + 1))
done