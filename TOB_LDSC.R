### Scripts to run LDSC analysis #######



for i in TB2022.*.Human.DownRegulated.txt;do
    basefile_name=$(basename $i .txt)
    for((j=1; j<=22; j=j+1));do
    ./make_annot.py \
    --gene-set-file $i \
    --gene-coord-file /rds/general/user/rf1116/home/LDSC_2022/REF/ENSG_coord.txt \
    --windowsize 100000 \
    --bimfile /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j.bim \
    --annot-file $basefile_name.$j.annot.gz

    mv $basefile_name.$j.annot.gz /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_DOWN
    done
done



for i in TB2022.N.GAERS.vs.N.NEC.Human.DownRegulated.txt;do
    basefile_name=$(basename $i .txt)
    for((j=1; j<=22; j=j+1));do
    ./ldsc.py \
    --l2 \
    --bfile /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_plink/1000G.EUR.QC.$j \
    --ld-wind-cm 1 \
    --thin-annot \
    --annot /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_DOWN/$basefile_name.$j.annot.gz \
    --out $basefile_name.$j \
    --print-snps /rds/general/user/rf1116/home/LDSC_2022/REF/hapmap3_snps/hm.$j.snp

    mv $basefile_name.$j.l2.M $basefile_name.$j.l2.M_5_50 $basefile_name.$j.l2.ldscore.gz -t /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_DOWN

	done
done

for i in TB2022.*.Human.DownRegulated.txt;do
basefile_name=$(basename $i .txt)
./ldsc.py \
--h2 /rds/general/user/rf1116/home/LDSC_2022/SumStat/SavageJansen_2018.sumstats.gz \
--w-ld-chr /rds/general/user/rf1116/home/LDSC_2022/REF/weights_hm3_no_hla/weights. \
--ref-ld-chr /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Annot_TB_DOWN/$basefile_name.,/rds/general/user/rf1116/home/LDSC_2022/REF/1000G_EUR_Phase3_baselineNEW/baseline. \
--overlap-annot \
--frqfile-chr /rds/general/user/rf1116/home/LDSC_2022/REF/1000G_Phase3_frq/1000G.EUR.QC. \
--out $basefile_name \
--print-coefficients
mv $basefile_name.results /rds/general/user/rf1116/home/LDSC_2022/TB_DOWN/Results_TB_DOWN_IQ
done

###########
