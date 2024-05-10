#####宏基因组数据分析流程#####
###################################################################################################
###质控
fastp  --thread  16 -i  01_raw_data/A1_1.fq.gz  -I 01_raw_data/A1_2.fq.gz -o  02_clean_data/A1_1.fq.gz -O 02_clean_data/A1_2.fq.gz  -j  02_clean_data/A1.fastp.json -h  02_clean_data/A1.fastp.html
fastp  --thread  16 -i  01_raw_data/A2_1.fq.gz  -I 01_raw_data/A2_2.fq.gz -o  02_clean_data/A2_1.fq.gz -O 02_clean_data/A2_2.fq.gz  -j  02_clean_data/A2.fastp.json -h  02_clean_data/A2.fastp.html
fastp  --thread  16 -i  01_raw_data/A3_1.fq.gz  -I 01_raw_data/A3_2.fq.gz -o  02_clean_data/A3_1.fq.gz -O 02_clean_data/A3_2.fq.gz  -j  02_clean_data/A3.fastp.json -h  02_clean_data/A3.fastp.html
fastp  --thread  16 -i  01_raw_data/B1_1.fq.gz  -I 01_raw_data/B1_2.fq.gz -o  02_clean_data/B1_1.fq.gz -O 02_clean_data/B1_2.fq.gz  -j  02_clean_data/B1.fastp.json -h  02_clean_data/B1.fastp.html
fastp  --thread  16 -i  01_raw_data/B2_1.fq.gz  -I 01_raw_data/B2_2.fq.gz -o  02_clean_data/B2_1.fq.gz -O 02_clean_data/B2_2.fq.gz  -j  02_clean_data/B2.fastp.json -h  02_clean_data/B2.fastp.html
fastp  --thread  16 -i  01_raw_data/B3_1.fq.gz  -I 01_raw_data/B3_2.fq.gz -o  02_clean_data/B3_1.fq.gz -O 02_clean_data/B3_2.fq.gz  -j  02_clean_data/B3.fastp.json -h  02_clean_data/B3.fastp.html
fastp  --thread  16 -i  01_raw_data/C1_1.fq.gz  -I 01_raw_data/C1_2.fq.gz -o  02_clean_data/C1_1.fq.gz -O 02_clean_data/C1_2.fq.gz  -j  02_clean_data/C1.fastp.json -h  02_clean_data/C1.fastp.html
fastp  --thread  16 -i  01_raw_data/C2_1.fq.gz  -I 01_raw_data/C2_2.fq.gz -o  02_clean_data/C2_1.fq.gz -O 02_clean_data/C2_2.fq.gz  -j  02_clean_data/C2.fastp.json -h  02_clean_data/C2.fastp.html
fastp  --thread  16 -i  01_raw_data/C3_1.fq.gz  -I 01_raw_data/C3_2.fq.gz -o  02_clean_data/C3_1.fq.gz -O 02_clean_data/C3_2.fq.gz  -j  02_clean_data/C3.fastp.json -h  02_clean_data/C3.fastp.html
###################################################################################################
#汇总质控文件
conda install bioconda::multiqc
multiqc ./02_clean_data/ #生成报告
mv multiqc_report_1.html multiqc_data 
mv multiqc_data 03_multiqc_data 
###################################################################################################
###去宿主(根据自己的样本类型添加相应物种基因组)
bowtie2-build genome.fa genome.db #构建宿主基因组index
bowtie2 --threads 60 -x ./genome.db -1 02_clean_data/A1_1.fq.gz -2 02_clean_data/A1_2.fq.gz -S A1.sam #bowtie2比对
samtools view -f 12 A1.sam > A1.unmap.bam #去除宿主数据
samtools fastq  -1 A1_1.clean.fq.gz -2 A1_2.clean.fq.gz -s A1_singleton.clean.fq.gz A1.unmap.bam #bam转换回fq格式
#上述步骤可进行串联操作
#bowtie2 --threads 60 -x ./genome.db -1 02_clean_data/A1_1.fq.gz -2 02_clean_data/A1_2.fq.gz | samtools view -f 12 | samtools fastq -1  A2_1.clean.fq.gz -2 A2_2.clean.fq.gz -s A2_singleton.clean.fq.gz 
###################################################################################################
#基于reads物种注释
conda install bioconda::kraken2 #安装kraken2做物种注释
https://benlangmead.github.io/aws-indexes/k2 #网址
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20240112.tar.gz #下载数据库、上传服务器
tar -zxvf k2_pluspfp_16gb_20240112.tar.gz #解压缩
#物种注释（reads）
mkdir 04_annotation
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/A1.kreport --output 04_annotation/A1.kraken  02_clean_data/A1_1.fq.gz 02_clean_data/A1_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/A2.kreport --output 04_annotation/A2.kraken  02_clean_data/A2_1.fq.gz 02_clean_data/A2_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/A3.kreport --output 04_annotation/A3.kraken  02_clean_data/A3_1.fq.gz 02_clean_data/A3_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/B1.kreport --output 04_annotation/B1.kraken  02_clean_data/B1_1.fq.gz 02_clean_data/B1_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/B2.kreport --output 04_annotation/B2.kraken  02_clean_data/B2_1.fq.gz 02_clean_data/B2_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/B3.kreport --output 04_annotation/B3.kraken  02_clean_data/B3_1.fq.gz 02_clean_data/B3_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/C1.kreport --output 04_annotation/C1.kraken  02_clean_data/C1_1.fq.gz 02_clean_data/C1_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/C2.kreport --output 04_annotation/C2.kraken  02_clean_data/C2_1.fq.gz 02_clean_data/C2_2.fq.gz --confidence 0.2
kraken2 --threads 64 --paired --db database/K2/  --report 04_annotation/C3.kreport --output 04_annotation/C3.kraken  02_clean_data/C3_1.fq.gz 02_clean_data/C3_2.fq.gz --confidence 0.2
###################################################################################################
#物种丰度估计
conda install bioconda::bracken #安装bracken软件做丰度定量、安装后显示找不到软件
conda remove bracken #删除软件
https://github.com/jenniferlu717/Bracken #在github上下载、上传至服务器、在该文件目录下运行以下脚本安装软件
bash install_bracken.sh #在software文件夹下面
/icdc/Users/caozhijie/metagenomics/database/Bracken-master/bracken --help #显示安装成功
export PATH="/icdc/Users/caozhijie/metagenomics/database/Bracken-master:$PATH" #添加到环境变量
-----------------------------------------------------------------------------------------------------------------------
mkdir 05_abundance/
bracken -d database/K2/ -i annotation/A1.kreport -o 05_abundance/A1.bracken.S -w 05_abundance/A1.kreport -l S -t 10  #-l S 表示分类水平，如 D,P,C,O,F,G,S
bracken -d database/K2/ -i annotation/A2.kreport -o 05_abundance/A2.bracken.S -w 05_abundance/A2.kreport -l S -t 10
bracken -d database/K2/ -i annotation/A3.kreport -o 05_abundance/A3.bracken.S -w 05_abundance/A3.kreport -l S -t 10
bracken -d database/K2/ -i annotation/B1.kreport -o 05_abundance/B1.bracken.S -w 05_abundance/B1.kreport -l S -t 10
bracken -d database/K2/ -i annotation/B2.kreport -o 05_abundance/B2.bracken.S -w 05_abundance/B2.kreport -l S -t 10
bracken -d database/K2/ -i annotation/B3.kreport -o 05_abundance/B3.bracken.S -w 05_abundance/B3.kreport -l S -t 10
bracken -d database/K2/ -i annotation/C1.kreport -o 05_abundance/C1.bracken.S -w 05_abundance/C1.kreport -l S -t 10
bracken -d database/K2/ -i annotation/C2.kreport -o 05_abundance/C2.bracken.S -w 05_abundance/C2.kreport -l S -t 10
bracken -d database/K2/ -i annotation/C3.kreport -o 05_abundance/C3.bracken.S -w 05_abundance/C3.kreport -l S -t 10
#汇总各样本丰度结果
conda install bioconda::kraken-biom
kraken-biom 05_abundance/*.kreport --max D  -o  ./05_abundance/S.biom #转换格式
biom  convert -i  ./05_abundance/S.biom -o  ./05_abundance/S.count.tsv.tmp  --to-tsv --header-key taxonomy #汇总
sed 's/; g__\([^;]\+\); s__/; g__\1; s__\1 /'   ./05_abundance/S.count.tsv.tmp >  ./05_abundance/S.taxID.count.tsv #修改格式
sed  '/^#/! s/^[0-9]\+\t\(.*[A-Za-z]\+__\([^;]\+\)\)$/\2\t\1/'  ./05_abundance/S.taxID.count.tsv  >  ./05_abundance/S.taxName.count.tsv
sed '1d; 2s/^#//' ./05_abundance/S.taxName.count.tsv |awk  -F "\t" -v 'OFS=\t'  '{$NF = ""; { print $0 }}' | sed 's/\t$//' > ./05_abundance/S.count.tsv
-----------------------------------------------------------------------------------------------------------------------
R包安装
conda install -c r r  #conda安装R、与浏览器打开的RStudio Serve不是同一个
conda install r-stringi 
conda install -c conda-forge r-tidyr
conda install -c bioconda r-pheatmap
conda install -c conda-forge r-ggalluvial 
BiocManager::install("phyloseq"),
BiocManager::install("PCAtools")
在Linux系统中,使用conda安装的R语言包和使用R环境中的install.packages()安装的R包可能会位于不同的路径下,会导致在命令行中无法直接访问到R环境中安装的包
解决方法：
确定R环境中包的路径 .libPaths()
将所有路径添加至环境变量(冒号分隔不同路径)、这样浏览器中访问R中安装的R包即可在命令行中使用
vim ~/.bashrc
export R_LIBS="/icdc/home/yangjing/R/x86_64-pc-linux-gnu-library/4.3:/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library"
source ~/.bashrc 
-------------------------------------------------  
Rscript script/draw_taxonBarplot.R   05_abundance/S.count.tsv  20   05_abundance/S.count.out #使用conda安装R及R包 (丰度前20的物种)
Rscript script/draw_taxonBarplot1.R   05_abundance/S.count.tsv  20   05_abundance/S.count.out #两个脚本的区别再与柱状图之间是否连线
mkdir species_abundance &&  mv S* species_abundance
#绘制krona环形图
conda install bioconda::krona
kreport2krona.py -r A1.kreport -o A1.kreport.krona #生成数据
ktImportText -o A1.kreport.krona.html A1.kreport.krona #绘图
###################################################################################################
#物种多样性分析
mkdir 06_diversity
ln -s ../05_abundance/species_abundance/S.biom
Rscript  ../script/diversity_alpha.R  S.biom  S.alpha  
Rscript  ../script/diversity_beta.R  S.biom  sample.txt  S.beta  #PCA和PCoA分析及绘图
###################################################################################################
#差异物种分析
mkdir 07_differences
#使用 lefse进行差异物种分析
conda install bioconda::lefse
ln -s  ../05_abundance/species_abundance/S.taxName.count.tsv
awk -F "\t" -v 'OFS=\t' '{$1=$NF; $NF=""; {print $0}}' S.taxName.count.tsv |sed 1d|sed 's/; /|/g'   > S.lefse.tmp
head -n 1 S.lefse.tmp |sed 's/\t/\n/g' |sed '1d'  |while read sp ;do awk '$2=="'$sp'" {printf "\t"$1 }' sample.txt ;done | awk '{print "class"$0 }' > S.lefse.header
cat S.lefse.header S.lefse.tmp   >   S.lefse
lefse_format_input.py   S.lefse  S.lefse.in -c 1 -s -1 -u 2 -o 1000000
lefse_run.py  S.lefse.in  S.lefse.out  
lefse_plot_res.py S.lefse.out S.lefse.LDA.pdf --format pdf
lefse_plot_cladogram.py S.lefse.out  S.lefse.cladogram.pdf --format pdf --labeled_start_lev 1
-------------------------------------------------
#使用DESeq2/edgeR进行差异物种分析
同时拷贝PerlLib和R文件夹及run_DE_analysis.pl脚本至script目录下
ln -s ../05_abundance/species_abundance/S.count.tsv
sed 's/ /_/g' S.count.tsv > S.count.matrix
perl ../script/run_DE_analysis.pl --matrix S.count.matrix --method DESeq2 --samples_file sample.txt --output DESeq2_out  
###################################################################################################
#组装(08_assemble)
megahit -1 02_clean_data/A1_1.fq.gz -2 02_clean_data/A1_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/A1_megahit --out-prefix A1
megahit -1 02_clean_data/A2_1.fq.gz -2 02_clean_data/A2_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/A2_megahit --out-prefix A2
megahit -1 02_clean_data/A3_1.fq.gz -2 02_clean_data/A3_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/A3_megahit --out-prefix A3
megahit -1 02_clean_data/B1_1.fq.gz -2 02_clean_data/B1_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/B1_megahit --out-prefix B1
megahit -1 02_clean_data/B2_1.fq.gz -2 02_clean_data/B2_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/B2_megahit --out-prefix B2
megahit -1 02_clean_data/B3_1.fq.gz -2 02_clean_data/B3_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/B3_megahit --out-prefix B3
megahit -1 02_clean_data/C1_1.fq.gz -2 02_clean_data/C1_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/C1_megahit --out-prefix C1
megahit -1 02_clean_data/C2_1.fq.gz -2 02_clean_data/C2_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/C2_megahit --out-prefix C2
megahit -1 02_clean_data/C3_1.fq.gz -2 02_clean_data/C3_2.fq.gz --min-contig-len 500  --out-dir 08_assemble/C3_megahit --out-prefix C3
#组装结果统计评估（省略）
conda install bioconda::quast
quast.py ./*.fa
###################################################################################################
#contig物种注释（基于unigene做物种注释同理）省略
diamond blastp -d nr.dmnd -q contig.fa -f 6 --id 90 -k 1 -e 1e-5 -c 1  -o result.tab #蛋白比对
wget https://software-ab.cs.uni-tuebingen.de/download/megan6/MEGAN_Ultimate_unix_6_25_9.sh #下载软件
sh MEGAN_Ultimate_unix_6_25_9.sh
wget https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip #下载GI号对应的taxanomy的mapping文件
unzip megan-map-Feb2022.db.zip
../MEGAN/megan/tools/blast2lca -i all_nr_match.txt -f BlastTab -ms 50 -me 0.000001 -a2t megan-map-Feb2022.db -o nr_blast_otu_tax.tsv #使用MEGAN LCA算法进行物种信息注释
###################################################################################################
#基因预测
mkdir 09_gene_prediction
prokka  --outdir A1_prokka --prefix A1 --addgenes --addmrna --locustag  A1 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/A1_megahit/A1.contigs.fa
prokka  --outdir A2_prokka --prefix A2 --addgenes --addmrna --locustag  A2 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/A2_megahit/A2.contigs.fa
prokka  --outdir A3_prokka --prefix A3 --addgenes --addmrna --locustag  A3 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/A3_megahit/A3.contigs.fa
prokka  --outdir B1_prokka --prefix B1 --addgenes --addmrna --locustag  B1 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/B1_megahit/B1.contigs.fa
prokka  --outdir B2_prokka --prefix B2 --addgenes --addmrna --locustag  B2 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/B2_megahit/B2.contigs.fa
prokka  --outdir B3_prokka --prefix B3 --addgenes --addmrna --locustag  B3 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/B3_megahit/B3.contigs.fa
prokka  --outdir C1_prokka --prefix C1 --addgenes --addmrna --locustag  C1 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/C1_megahit/C1.contigs.fa
prokka  --outdir C2_prokka --prefix C2 --addgenes --addmrna --locustag  C2 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/C2_megahit/C2.contigs.fa
prokka  --outdir C3_prokka --prefix C3 --addgenes --addmrna --locustag  C3 --kingdom  Bacteria --metagenome --cpus 64  ../08_assemble/C3_megahit/C3.contigs.fa
###################################################################################################
#基因预测结果去冗余(10_unigene)
mkdir 10_unigene
cat ../10_unigene/*_prokka/*.ffn > all.trans.fa
cat ../10_unigene/*_prokka/*.faa > all.pep.fa
conda install bioconda::seqtk
seqtk subseq  all.trans.fa all.cds.id > all.cds.fa
conda install bioconda::cd-hit
cd-hit-est -i all.cds.fa -o  all.cds.cdhit -c 0.95 -aS 0.9 -M 1000000 -T 64
cp  all.cds.cdhit  unigene_cds.fasta
awk '$1 ~/^>/{print $1}'  all.cds.cdhit | sed 's/^>//' > unigene.id
seqtk  subseq  all.pep.fa unigene.id > unigene_pep.fasta
awk '{if($1~/^>/){printf $1 $2} else if($NF ~ /*$/){print "\t"$3}}' all.cds.cdhit.clstr |sed 's/>//g; s/...$//'|awk '{print $2"\t"$1}'  > map_id.txt
chmod +x ../script/map_fasta_ids
perl ../script/map_fasta_ids map_id.txt  unigene_cds.fasta
perl ../script/map_fasta_ids map_id.txt  unigene_pep.fasta 
###################################################################################################
#eggnog注释(11_eggnog)
conda install bioconda::eggnog-mapper #下载软件
http://eggnog6.embl.de/download/emapperdb-5.0.2/ #数据库网址
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz #下载数据库
nohup gzip -d eggnog.db.gz &
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz #下载数据库
tar -zxvf eggnog.taxa.tar.gz
mkdir /icdc/home/yangjing/miniconda3/lib/python3.11/site-packages/data
nohup download_eggnog_data.py -y -f &   #自动下载数据库速度太慢、将前两步的数据库和下一步数据库移动到data目录下
wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz 
gunzip  eggnog_proteins.dmnd.gz
emapper.py -i 10_unigene/unigene_pep.fasta --output unigene_pep.eggnogmap -m diamond --cpu 64 --output_dir 11_eggnog/ --evalue 1e-5  #时间较长建议挂载后台
需要在script文件目录下存放emcp目录及其中kegg_info.RData、cog_funclass.tab文件、相关R包安装方法见install.R文件
Rscript ../script/emcp/emapperx.R  unigene_pep.eggnogmap.emapper.annotations ../10_unigene/unigene_pep.fasta #对结果进行统计绘图并构建orgdb
###################################################################################################
#uniprot注释(12_uniprot)
nohup wget  https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz  &  
#下载uniref90序列文件43G（2024-03-27）
nohup wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz &
nohup gzip -d idmapping_selected.tab.gz &
#下载idmapping文件11G（2024-03-27）
diamond  makedb --in  ./uniref90.fasta.gz  --db  uniref90
#比对
nohup diamond blastp  --db  ../database/uniprot/uniref90 --query ../10_unigene/unigene_pep.fasta --out unigene_pep.uniref90.m6  --outfmt 6 --max-target-seqs 1 --evalue 1e-5  &
perl ../script/uniref90_idmapping.pl unigene_pep.uniref90.m6  ../database/uniprot/idmapping_selected.tab > unigene_pep.uniref90.GOanno #基于idmapping 提取GO注释信息
Rscript  ../script/emcp/GOmapperx.R  unigene_pep.uniref90.GOanno # 构建orgdb，并统计GO 注释
###################################################################################################
#CARDA数据库(13_card_rgi)
conda create --name rgi rgi
conda activate rgi
rgi load  -i card.json --local 
#在数据库中运行此命令
rgi main --input_sequence ../10_unigene/unigene_pep.fasta --output_file rgi_out --input_type protein --alignment_tool DIAMOND --num_threads 64 --local --clean --include_loose
cut -f 1,6,9,11,15,16,17 rgi_out.txt > rgi_out.txt.info
conda deactivate
###################################################################################################
#碳水化合物酶数据库注释(14_cazydb)
#软件
conda create -n run_dbcan python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
conda activate run_dbcan
conda install run-dbcan==2.0.11  #只有这一个版本
#数据库
wget https://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07262023.fa #下载数据库
diamond makedb --in CAZyDB.07262023.fa -d CAZy
wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt  && hmmpress dbCAN.txt  #需下载完整不然会报错
wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa  && diamond makedb --in tcdb.fa -d tcdb 
wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm 
wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm 
wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm
#运行
run_dbcan.py  --out_dir cazy_diamond --db_dir ../database/CAZy_dbCAN2/ --tools diamond  --dia_cpu 64 ../10_unigene/unigene_pep.fasta  protein
#hmm比对用：run_dbcan --out_dir cazy_hmmer --db_dir /app/db  --tools hmmer  --hmm_cpu 64 unigene_pep.fasta  protein 
#结果整理
sed '1d' cazy_diamond/overview.txt | cut -f 1,4 | sed 's/+/\t/g' | awk '{for(i=2;i<=NF; i++){print $1"\t"$i}}' | awk '$2 ~ /^[A-Z]/'|sed 's/\t\(\([A-Z]\+\).*\)$/\t\1\t\2/' > cazy_diamond.family.txt 
le cazy_diamond.family.txt |cut -f 1,3|sort -u |cut -f 2 |sort |uniq -c
###################################################################################################
#基因丰度估计(15_gene_abundance)
mkdir 15_gene_abundance
ln -s ../10_unigene/unigene_cds.fasta
salmon index -t  unigene_cds.fasta  -i  unigene_index #建索引
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/A1_1.fq.gz -2 02_clean_data/A1_2.fq.gz -o  quants/A1.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/A2_1.fq.gz -2 02_clean_data/A2_2.fq.gz -o  quants/A2.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/A3_1.fq.gz -2 02_clean_data/A3_2.fq.gz -o  quants/A3.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/B1_1.fq.gz -2 02_clean_data/B1_2.fq.gz -o  quants/B1.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/B2_1.fq.gz -2 02_clean_data/B2_2.fq.gz -o  quants/B2.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/B3_1.fq.gz -2 02_clean_data/B3_2.fq.gz -o  quants/B3.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/C1_1.fq.gz -2 02_clean_data/C1_2.fq.gz -o  quants/C1.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/C2_1.fq.gz -2 02_clean_data/C2_2.fq.gz -o  quants/C2.quant
salmon  quant --validateMappings  --meta -p 64   -i  unigene_index  -l IU   -1 02_clean_data/C3_1.fq.gz -2 02_clean_data/C3_2.fq.gz -o  quants/C3.quant
sh combine.sh #合并表格
-------------------------------------------------
#combine.sh
quants=` ls quants/ | awk '{print "quants/"$1 }' | tr '\n'  ' ' `
names=` ls quants/ | sed 's/.quant$//'| tr '\n'   ' '  `
salmon quantmerge  --quants $quants --names  $names  --column tpm -o unigenens.tpm
salmon quantmerge --quants $quants --names  $names  --column numreads -o unigenens.count
-------------------------------------------------
Rscript  ../script/draw_abundance.R unigenens.tpm  sample.txt  unigenens.tpm_fig #绘图（箱线图、相关性热图、密度图、 pca图）
###################################################################################################
####差异基因分析
ln -s ../15_gene_abundance/unigenens.count
sed '1s/^Name\t//' unigenens.count > unigenens.count.matrix
perl ../script/run_DE_analysis.pl --matrix unigenens.count.matrix --method DESeq2 --samples_file  sample.txt  --output DE_out
##差异基因富集分析
ln -s ../11_eggnog/unigene_pep.eggnogmap.emapper.annotations
ln -s ../11_eggnog/R_Library/ ##链接eggnog注释时创建的R_Librar
Rscript ../script/enrich_analysis.R ./DE_out/unigenens.count.matrix.A_vs_C.DESeq2.DE_results ./unigene_pep.eggnogmap.emapper.annotations ./DE_out/unigenens.count.matrix.A_vs_C.DESeq2.DE_results.enrich
###################################################################################################

###分箱(17_bin)
find ../08_assemble -maxdepth 2 -name '*.contigs.fa' -exec ln -s {} . \; #链接组装结果
find ../02_clean_data -name '*.fq.gz' -exec sh -c 'gzip -dc "$0" > "$(basename "$0" .fq.gz).fastq"' {} \; #解压所有clean_data
-------------------------------------------------------------------------------软件安装
conda install -y mamba 
mamba create -y --name metawrap-env --channel ursky metawrap-mg=1.3.2
conda activate metawrap-env
conda install -y blas=2.5=mkl
#配置数据库
cd MY_CHECKM_FOLDER
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz
cd ../
checkm data setRoot /MY_CHECKM_FOLDER
-------------------------------------------------------------------------------
#可基于三种方法Bin（此处以2种方法演示）
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  A1.contigs.fa -o A1_metawrap_bin A1_1.fastq A1_2.fastq  
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  A2.contigs.fa -o A2_metawrap_bin A2_1.fastq A2_2.fastq
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  A3.contigs.fa -o A3_metawrap_bin A3_1.fastq A3_2.fastq
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  B1.contigs.fa -o B1_metawrap_bin B1_1.fastq B1_2.fastq  
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  B2.contigs.fa -o B2_metawrap_bin B2_1.fastq B2_2.fastq
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  B3.contigs.fa -o B3_metawrap_bin B3_1.fastq B3_2.fastq
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  C1.contigs.fa -o C1_metawrap_bin C1_1.fastq C1_2.fastq  
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  C2.contigs.fa -o C2_metawrap_bin C2_1.fastq C2_2.fastq
metaWRAP binning  -t 64 --metabat2 --maxbin2 --run-checkm  -a  C3.contigs.fa -o C3_metawrap_bin C3_1.fastq C3_2.fastq
#Bin提纯、多个软件结果整合（区别多样本结果整合）
metaWRAP  bin_refinement -t  64 -o A1_metawrap_refine -A A1_metawrap_bin/maxbin2_bins/ -B A1_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o A2_metawrap_refine -A A2_metawrap_bin/maxbin2_bins/ -B A2_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o A3_metawrap_refine -A A3_metawrap_bin/maxbin2_bins/ -B A3_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o B1_metawrap_refine -A B1_metawrap_bin/maxbin2_bins/ -B B1_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o B2_metawrap_refine -A B2_metawrap_bin/maxbin2_bins/ -B B2_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o B3_metawrap_refine -A B3_metawrap_bin/maxbin2_bins/ -B B3_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o C1_metawrap_refine -A C1_metawrap_bin/maxbin2_bins/ -B C1_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o C2_metawrap_refine -A C2_metawrap_bin/maxbin2_bins/ -B C2_metawrap_bin/metabat2_bins/
metaWRAP  bin_refinement -t  64 -o C3_metawrap_refine -A C3_metawrap_bin/maxbin2_bins/ -B C3_metawrap_bin/metabat2_bins/

#将所有分箱fasta文件存放在bin目录下、运行dRep去冗余(dRep无法下载)
conda create -n gtdbtk
conda activate drep
conda install bioconda::drep #安装失败
dRep dereplicate  out_dRep --length 50000 -comp 75 -con 25 -g ./bins/*.fasta 

#物种注释（18_bin_annotation）
conda create -n gtdbtk
conda activate gtdbtk
conda install bioconda::gtdbtk
nohup wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz & #数据库下载（2024-04-18）101G
#conda env config vars set GTDBTK_DATA_PATH="/icdc/Database/release220" #数据库挂载
#gtdbtk check_install #检查是否挂载成功
gtdbtk classify_wf --genome_dir  metawrap_70_10_bins  --out_dir  out_gtdbtk --cpus 60 --extension fa --tmpdir ./ --mash_db /icdc/Database/release220
#out_gtdbtk/gtdbtk.bac120.summary.tsv ：bin 物种分类信息
#out_gtdbtk/gtdbtk.bac120.classify.tree ：合并其他物种构建的tree    
