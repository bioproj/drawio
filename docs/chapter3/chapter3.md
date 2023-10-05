学习资料
+ [github.com/biovcnet/topic-metagenomics](https://github.com/biovcnet/biovcnet.github.io/wiki/TOPIC%3A-Metagenomics)

## Lesson 2 -- Taxonomic Classification using k-mers

## Lesson 3 -- Taxonomic Classification using MinHash sketches

## Lesson 4 -- Read mapping and read taxonomic classification

## Lesson 5 -- Genome and Metagenome Assembly

需要软件
```bash
conda install -c bioconda bbmap
conda install -c bioconda spades
```
测试数据
+ [topic-metagenomics](https://github.com/biovcnet/topic-metagenomics)


解压文件
```bash
cd data
# unzip files, concatenate R1's and R2's, and rezip
for file in *.fastq.zip; do unzip ${file}; done
# 合并reads1和reads2
cat SRR5780888_1.fastq SRR5780889_1.fastq | gzip > dvh_1.fastq.gz
cat SRR5780888_2.fastq SRR5780889_2.fastq | gzip > dvh_2.fastq.gz
# rezip individual files
for file in SRR*.fastq; do gzip ${file}; done
```
此时得到文件
```
├── dvh_1.fastq.gz
└── dvh_2.fastq.gz
```
使用了 BBMap 工具中的 `bbduk.sh` 脚本，对名为 dvh_1.fastq.gz 和 dvh_2.fastq.gz 的输入文件进行了质量修剪和适配器修剪，并将修剪后的序列写入名为 trimmed.fq.gz 的输出文件中
```bash
bbduk.sh \
  in=dvh_1.fastq.gz  \
  in2=dvh_2.fastq.gz  \
  out=trimmed.fq.gz \
  ktrim=r \
  k=23 \
  mink=11 \
  hdist=1 \
  tbo \
  tpe \
  minlen=70 \
  ref=adapters \
  ordered \
  ow=t
```
+ ktrim=r：使用右端修剪模式，即从右端开始修剪低质量碱基。
+ k=23：设置 k-mer 大小为 23，用于指定修剪过程中的匹配碱基数。
+ mink=11：设置最小 k-mer 大小为 11，即修剪过程中的最小匹配碱基数。
+ hdist=1：设置允许的最大汉明距离为 1，即在修剪过程中允许的最大不匹配碱基数。
+ tbo：开启两端输出模式，将修剪后的序列写入两个输出文件。
+ tpe：开启成对模式，确保两个输出文件中的序列依然是成对的。
+ minlen=70：设置最小序列长度为 70，小于该长度的序列将被丢弃。
+ ref=adapters：指定参考序列文件为 adapters，用于识别和修剪测序时可能存在的适配器序列。
+ ordered：按顺序处理输入文件，确保输入文件中的序列顺序与输出文件中的序列顺序相同。
+ ow=t：设置覆盖写入模式为真，即如果输出文件已经存在，则覆盖写入。

此时得到文件
```
trimmed.fq.gz
```

使用 BBMap 工具中的 bbduk.sh 脚本，对名为 trimmed.fq.gz 的输入文件进行了过滤，并将过滤后的序列写入名为 filtered.fq.gz 的输出文件中。
```bash
bbduk.sh in=trimmed.fq.gz \
  out=filtered.fq.gz \
  k=31 \
  ref=artifacts,phix \
  ordered \
  cardinality \
  ow=t
```
+ in=trimmed.fq.gz：指定输入文件为 trimmed.fq.gz。
+ out=filtered.fq.gz：指定输出文件为 filtered.fq.gz，过滤后的序列将写入该文件。
+ k=31：设置 k-mer 大小为 31，用于过滤过程中的匹配碱基数。
+ ref=artifacts,phix：指定过滤时使用的参考序列文件为 artifacts 和 phix，用于识别和过滤可能存在的序列污染或引+ 物序列。
+ ordered：按顺序处理输入文件，确保输入文件中的序列顺序与输出文件中的序列顺序相同。
+ cardinality：计算过滤后的序列的基数（distinct sequences）并输出。
+ ow=t：设置覆盖写入模式为真，即如果输出文件已经存在，则覆盖写入。
此时得到文件
```
filtered.fq.gz
```

使用 BBMap 工具中的 bbduk.sh 脚本，对名为 filtered.fq.gz 的输入文件进行了质量修剪，并将修剪后的序列写入名为 dvh_qtrimmed.fq.gz 的输出文件中
```bash
bbduk.sh in=filtered.fq.gz \
  out=dvh_qtrimmed.fq.gz \
  qtrim=r \
  trimq=10 \
  minlen=70 \
  ordered \
  maxns=0 \
  maq=8 \
  entropy=.95 \
  ow=t
```
+ in=filtered.fq.gz：指定输入文件为 filtered.fq.gz。
+ out=dvh_qtrimmed.fq.gz：指定输出文件为 dvh_qtrimmed.fq.gz，质量修剪后的序列将写入该文件。
+ qtrim=r：使用右端修剪模式，即从右端开始修剪低质量碱基。
+ trimq=10：指定质量阈值为 10，即低于该阈值的碱基将被修剪。
+ minlen=70：设置最小序列长度为 70，小于该长度的序列将被丢弃。
+ ordered：按顺序处理输入文件，确保输入文件中的序列顺序与输出文件中的序列顺序相同。
+ maxns=0：设置最大 N 碱基数为 0，即不允许序列中存在 N 碱基。
+ maq=8：设置最小平均质量阈值为 8，小于该阈值的序列将被丢弃。
+ entropy=.95：设置最小熵阈值为 0.95，小于该阈值的序列将被丢弃。
+ ow=t：设置覆盖写入模式为真，即如果输出文件已经存在，则覆盖写入。
此时得到文件
```
dvh_qtrimmed.fq.gz
```

使用了BBMap工具中的tadpole.sh脚本，对名为dvh_qtrimmed.fq.gz的输入文件进行了序列组装，并将组装后的contigs序列保存为名为tadpole_contigs.fasta的输出文件
```
tadpole.sh \
  in=dvh_qtrimmed.fq.gz \
  out=tadpole_contigs.fasta \
  k=124 \
  ow=t \
  prefilter=2 \
  prepasses=auto
```
+ in=dvh_qtrimmed.fq.gz：指定输入文件为dvh_qtrimmed.fq.gz。
+ out=tadpole_contigs.fasta：指定输出文件为tadpole_contigs.fasta，组装后的contigs序列将保存在该文件中。
+ k=124：指定k-mer大小为124，用于组装过程中的碱基匹配。
+ ow=t：设置覆盖写入模式为真，即如果输出文件已经存在，则覆盖写入。
+ prefilter=2：设置预过滤器阈值为2，用于在组装前过滤低频错误k-mers。
+ prepasses=auto：设置预处理的自动预处理次数，即根据输入数据自动确定预处理次数。

此时得到文件
```
tadpole_contigs.fasta
```
使用了SPAdes工具中的spades.py脚本，对名为dvh_qtrimmed.fq.gz的输入文件进行了序列组装，并将组装结果保存在名为dvh_spades的输出目录中
```
spades.py \
  -o dvh_spades \
  --12 dvh_qtrimmed.fq.gz \
  --only-assembler
```
+ -o dvh_spades：指定输出目录为dvh_spades，组装结果将保存在该目录中。
+ --12 dvh_qtrimmed.fq.gz：指定配对的序列文件为dvh_qtrimmed.fq.gz。
+ --only-assembler：设置仅进行组装步骤，不进行错误校正和拼接步骤

使用了statswrapper.sh脚本，对dvh_spades/*.fasta目录下的FASTA格式文件进行统计分析
```
statswrapper.sh \
  dvh_spades/*.fasta \
  tadpole_contigs.fasta > dvh_stats.txt
```

| n_scaffolds | n_contigs | scaf_bp | contig_bp | gap_pct | scaf_N50 | scaf_L50 | ctg_N50 | ctg_L50 | scaf_N90 | scaf_L90 | ctg_N90 | ctg_L90 | scaf_max | ctg_max | scaf_n_gt50K | scaf_pct_gt50K | gc_avg  | gc_std  | filename                                                        |
| ----------- | --------- | ------- | --------- | ------- | -------- | -------- | ------- | ------- | -------- | -------- | ------- | ------- | -------- | ------- | ------------ | -------------- | ------- | ------- | --------------------------------------------------------------- |
| 215         | 215       | 3742024 | 3742024   | 0.000   | 14       | 88810    | 14      | 88810   | 45       | 19130    | 45      | 19130   | 284425   | 284425  | 27           | 74.346         | 0.63322 | 0.05737 | /workspaces/topic-metagenomics/data2/dvh_spades/before_rr.fasta |
| 105         | 105       | 3742672 | 3742672   | 0.000   | 9        | 165421   | 9       | 165421  | 22       | 52191    | 22      | 52191   | 285836   | 285836  | 22           | 90.559         | 0.63321 | 0.05773 | /workspaces/topic-metagenomics/data2/dvh_spades/contigs.fasta   |
| 99          | 105       | 3743107 | 3742672   | 0.012   | 8        | 168161   | 9       | 165421  | 20       | 52191    | 22      | 52191   | 408455   | 285836  | 20           | 90.552         | 0.63321 | 0.05884 | /workspaces/topic-metagenomics/data2/dvh_spades/scaffolds.fasta |
| 172         | 172       | 52664   | 52664     | 0.000   | 71       | 271      | 71      | 271     | 153      | 252      | 153     | 252     | 3188     | 3188    | 0            | 0.000          | 0.61376 | 0.05661 | /workspaces/topic-metagenomics/data2/tadpole_contigs.fasta      |


完整代码
```bash
for prefix in `ls *_1.fastq.gz | cut -f1 -d'_' | sort -u`; do 

  echo ${prefix}

  R1=( ${prefix}*_1.fastq.gz )
  R2=( ${prefix}*_2.fastq.gz )
    
  #Trim adapters
  # 'ordered' means to maintain the input order as produced by clumpify.sh
  bbduk.sh in=${R1} in2=${R2} out=trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ordered ow=t

  #Remove synthetic artifacts and spike-ins by kmer-matching
  # 'cardinality' will generate an accurate estimation of the number of unique kmers in the dataset using the LogLog algorithm
  bbduk.sh in=trimmed.fq.gz out=filtered.fq.gz k=31 ref=artifacts,phix ordered cardinality ow=t
  
  #Quality-trim and entropy filter the remaining reads.
  # 'entropy' means to filter out reads with low complexity
  # 'maq' is 'mininum average quality' to filter out overall poor reads
  bbduk.sh in=filtered.fq.gz out=${prefix}_qtrimmed.fq.gz qtrim=r trimq=10 minlen=70 ordered maxns=0 maq=8 entropy=.95 ow=t

  # Assembly using tadpole
  tadpole.sh in=${prefix}_qtrimmed.fq.gz out=tadpole_contigs.fasta k=124 ow=t prefilter=2 prepasses=auto
  
  # Assembly quality-trimmed reads using SPAdes
  spades.py -o ${prefix}_spades --12 ${prefix}_qtrimmed.fq.gz --only-assembler
  
  # calculate assembly statistics
  statswrapper.sh ${prefix}_spades/*.fasta tadpole_contigs.fasta > ${prefix}_stats.txt
  
  cat ${prefix}_stats.txt
  
  # remove extra files
  # rm trimmed.fq.gz filtered.fq.gz
  
done
```

## Lesson 6 -- Binning Metagenome-assembled Genomes

## Lesson 7 -- Bin Evaluation
