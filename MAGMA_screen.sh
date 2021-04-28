inFile=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/MAGMA/ukb2_sumstats/f.50.0.0_res.EUR.sumstats.MACfilt.txt

geneannot=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/MAGMA/f.50.0.0_res.EUR.sumstats.MACfilt/f.50.0.0_res.EUR.sumstats.MACfilt.genes.annot
generesults=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/MAGMA/f.50.0.0_res.EUR.sumstats.MACfilt/f.50.0.0_res.EUR.sumstats.MACfilt.genes.raw
geneset=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/MAGMA/genesets/signature_conserved.txt

outPath=/share/ScratchGeneral/scoyou/scott_projects/osteochondro_signature/data/humanchondro_bulk/project_results/MAGMA/f.50.0.0_res.EUR.sumstats.MACfilt
mkdir $outPath
outPref=$outPath/f.50.0.0_res.EUR.sumstats.MACfilt

geneloc=/share/ScratchGeneral/scoyou/local/bin/MAGMA-1.09/reference_files/NCBI37.3.gene.loc
bfile=/share/ScratchGeneral/scoyou/local/bin/MAGMA-1.09/reference_files/g1000_eur

gunzip $inFile'.gz'
/share/ScratchGeneral/scoyou/local/bin/MAGMA-1.09/magma --annotate window=1 --gene-loc $geneloc --out $outPref  --snp-loc $inFile
/share/ScratchGeneral/scoyou/local/bin/MAGMA-1.09/magma --bfile $bfile --gene-annot $geneannot --gene-model snp-wise=mean --seed 1 --pval $inFile use=SNPID_UKB,P duplicate=first ncol=NMISS --out $outPref
/share/ScratchGeneral/scoyou/local/bin/MAGMA-1.09/magma --gene-results $generesults --set-annot $geneset --out $outPref
gzip $inFile

