######### Scripts for LDSC analysis #########

 
for i in *.txt;do
    basefile_name=$(basename $i .txt)
    for((j=1; j<=22; j=j+1));do
    ./make_annot.py \
    --gene-set-file $i \
    --gene-coord-file /rds/general/user/rf1116/home/LDSC_2022/REF/ENSG_coord.txt \
    --windowsize 100000 \
    --bimfile /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j.bim \
    --annot-file $basefile_name.$j.annot.gz

    mv $basefile_name.$j.annot.gz /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_ALL
    done
done



for i in *.txt;do
    basefile_name=$(basename $i .txt)
    for((j=1; j<=22; j=j+1));do
    ./ldsc.py \
    --l2 \
    --bfile /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j \
    --ld-wind-cm 1 \
    --thin-annot \
    --annot /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_ALL/$basefile_name.$j.annot.gz \
    --out $basefile_name.$j \
    --print-snps /rds/general/user/rf1116/home/LDSC_2022/REF/hapmap3_snps/hm.$j.snp

    mv $basefile_name.$j.l2.M $basefile_name.$j.l2.M_5_50 $basefile_name.$j.l2.ldscore.gz -t /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_ALL

	done
done

for i in TB2022.All.E.vs.All.N.Human.AllRegulated.txt;do
basefile_name=$(basename $i .txt)
./ldsc.py \
--h2 /rds/general/user/rf1116/home/LDSC_2022/SumStat/SavageJansen_2018.sumstats.gz \
--w-ld-chr /rds/general/user/rf1116/home/LDSC_2022/REF/weights_hm3_no_hla/weights. \
--ref-ld-chr /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_ALL/$basefile_name.,/rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_baselineNEW/baseline. \
--overlap-annot \
--frqfile-chr /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_Phase3_frq/1000G.EUR.QC. \
--out $basefile_name \
--print-coefficients
mv $basefile_name.results /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Results_TB_ALL_IQ
done



library(data.table)
library(ggplot2)
library(cowplot)

FOLDERS=list.files()
count=0
for(FOLD in FOLDERS){
	All.files=list.files(path=FOLD,pattern=".results")
	for(file in All.files){
		data=fread(sprintf("%s/%s",FOLD,file))
		tmp=data[1,]
		tmp$FOLD=FOLD
		tmp$file=file
		if(count==0){
			R.Table=tmp
		}else{
			R.Table=rbind(R.Table,tmp)
		}
		count=count+1
	}
}
R.Table$ct=gsub("_\\w[^_]*.results","",R.Table$file,perl=TRUE)
R.Table$ct=gsub(".results","",R.Table$file,perl=TRUE)
R.Table$ct=gsub(paste(FOLDERS,collapse="|"),"",R.Table$ct)
R.Table$ct=gsub("_$","",R.Table$ct)
ggplot(R.Table) + geom_bar(aes(x=ct,y=Enrichment_p,fill=FOLD),stat="identity",position="dodge") + scale_y_log10() + coord_flip() + theme_cowplot() + geom_hline(yintercept=0.05/dim(R.Table)[1])
