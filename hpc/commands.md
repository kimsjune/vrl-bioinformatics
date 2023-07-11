## Week 1

### relative vs. absolute paths

Go to your home directory, and make a subdirectory called `hello`.
```bash
cd /home/$USER
mkdir hello
```
Now, go to your scratch directory and make a subdirectory called `bye`.
```bash
cd /scratch/$USER
mkdir bye
```
What happens if you try to `cd` into `hello` folder from the current scratch directory?  
To go to the `hello` directory, you have to put in the full/absolute path.
```bash
cd /home/$USER/hello
pwd
```
If I am already in my home directory, I don't _have to_ write the full directory.
```
cd /home/$USER
cd hello
pwd
```

### / in front of directory
`/` refers to the root directory. We write `/scratch/$USER`, not `scratch/$USER` because the latter means that we are looking for a `scratch/$USER` directory _within_ the current directory (wherever that is). Such directory does not exist anywhere else but in the _root_ (`/`) directory. You can `cd /` and `ls` to see what's in the root directory. You won't have permission to do anything. If you try to do anything funny in the root directory, it might get reported to Compute Canada.


## Week 2  

### to open and edit the bashrc file

```bash
nano ~/.bashrc
```

### what to write in the file

```bash
module load StdEnv/2020 star
export WS="/scratch/$USER/workshop"
export ACC="def-YOUR_PI_NAME"
```  

`control+x`, `y` and `enter` to save and exit.  

### download genome fasta and annotation

```bash
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
```

### gunzip the two files

```bash
gunzip -k *
```

### index job script

```bash
#!/bin/bash
#SBATCH --account=def-fdick
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=STAR_index
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /scratch/$USER/workshop/stargenome --genomeFastaFiles /scratch/$USER/workshop/mm10.fa  --sjdbGTFfile /scratch/$USER/workshop/gencode.vM1.annotation.gtf --sjdbOverhang 74
```

### submit the index job

```bash
sbatch --account=$ACC index.sh
```


## Week 3

### pool reads together (NOT a batch job)

```bash
#!/bin/bash
cat ./NTC_D1_L001/* ./NTC_D1_L002/* ./NTC_D1_L003/* ./NTC_D1_L004/* > ./NTC_D1.fastq.gz
cat ./NTC_D2_L001/* ./NTC_D2_L002/* ./NTC_D2_L003/* ./NTC_D2_L004/* > ./NTC_D2.fastq.gz
cat ./NTC_D3_L001/* ./NTC_D3_L002/* ./NTC_D3_L003/* ./NTC_D3_L004/* > ./NTC_D3.fastq.gz
cat ./NTC_G1_L001/* ./NTC_G1_L002/* ./NTC_G1_L003/* ./NTC_G1_L004/* > ./NTC_G1.fastq.gz
cat ./NTC_G2_L001/* ./NTC_G2_L002/* ./NTC_G2_L003/* ./NTC_G2_L004/* > ./NTC_G2.fastq.gz
cat ./NTC_G3_L001/* ./NTC_G3_L002/* ./NTC_G3_L003/* ./NTC_G3_L004/* > ./NTC_G3.fastq.gz
cat ./DKO_D1_L001/* ./DKO_D1_L002/* ./DKO_D1_L003/* ./DKO_D1_L004/* > ./DKO_D1.fastq.gz
cat ./DKO_D2_L001/* ./DKO_D2_L002/* ./DKO_D2_L003/* ./DKO_D2_L004/* > ./DKO_D2.fastq.gz
cat ./DKO_D3_L001/* ./DKO_D3_L002/* ./DKO_D3_L003/* ./DKO_D3_L004/* > ./DKO_D3.fastq.gz
cat ./DKO_G1_L001/* ./DKO_G1_L002/* ./DKO_G1_L003/* ./DKO_G1_L004/* > ./DKO_G1.fastq.gz
cat ./DKO_G2_L001/* ./DKO_G2_L002/* ./DKO_G2_L003/* ./DKO_G2_L004/* > ./DKO_G2.fastq.gz
cat ./DKO_G3_L001/* ./DKO_G3_L002/* ./DKO_G3_L003/* ./DKO_G3_L004/* > ./DKO_G3.fastq.gz
```

### align job script

```bash
#!/bin/bash
#SBATCH --account=def-fdick
#SBATCH --time=30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=STAR_align
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_D1.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_D1
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_D2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_D2
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_D3.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_D3
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_G1.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_G1
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_G2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_G2
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/NTC_G3.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/NTC_G3
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_D1.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_D1
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_D2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_D2
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_D3.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_D3
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_G1.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_G1
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_G2.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_G2
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/DKO_G3.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/DKO_G3

```

### a for loop to concatenate all lanes together

```bash
for i in NTC DKO
do
        for j in D1 D2 D3 G1 G2 G3
        do
cat << EOF >> catAll.sh
cat ./${i}_${j}_L001/* ./${i}_${j}_L002/* ./${i}_${j}_L003/* ./${i}_${j}_L004/* > ./../reads/${i}_${j}.fastq.gz
EOF
done
done
```

save as `loop_cat.sh`, which means you should `nano loop_cat.sh` in the first place.

### change permissions to execute the _loop_ script, then run<sup>*</sup>
You are running `loop_cat.sh` to _generate the `cat` command_ for all your samples. The `cat` commands are written to `catAll.sh`. 

```bash
chmod +x ./loop_cat.sh
./loop_cat.sh
```

### confirm, change permissions, then run the actual `cat` commands

```bash
head catAll.sh
chmod +x ./catAll.sh
./catAll.sh
```

With real, full-size libraries, you would have to submit this as a slurm job.

### a for loop to align all samples

```bash
for i in NTC DKO
do
        for j in D1 D2 D3 G1 G2 G3
        do
cat << EOF >> alignAll.sh
STAR --runThreadN 8 --genomeDir /scratch/$USER/workshop/stargenome --readFilesIn /scratch/$USER/workshop/reads/${i}_${j}.fastq.gz --readFilesCommand zcat --outFileNamePrefix ./sam/${i}_${j}
EOF
done
done
```

save as `loop_align.sh`.
### change permission to execute the _loop_ script, then run

```bash
chmod +x loop_align.sh
./loop_align.sh
```

### unlike `cat` commands above, the aligning job should be submitted to slurm
This is because it will take too long.  
We need to add the following `SBATCH` lines to `alignAll.sh`. Remember that this is the actual `STAR` command file.  
Then submit by `sbatch --account=$ACC alignAll.sh`. 

```bash
#!/bin/bash
#SBATCH --account=def-fdick
#SBATCH --time=30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=STAR_align
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
```


## Week 4


### checking how well our samples have been aligned 


```bash
cd /scratch/$USER/workshop/sam
less NTC_D1Log.final.out
```
We want to see high overall mapping rates _and_ unique mapping rates.  

### you can also just look at .sam files because they are just text

```bash
head -n 100 ./sam/NTC_D1Aligned.out.sam
```

### for `htseq-count` and single-end libraries, we don't need to bother converting .sam to sorted .bam file, so we can go straight to _counting_.  
We (`htseq-count` program) are _counting_ how many reads align to exons of a given gene in our annotation file (.gtf)  
To use this program we need to install HTSeq, which is a python program. This is slightly different to `module load`.

```bash
module load python/3
virtualenv --no-download ENV
source ENV/bin/activate
pip install HTSeq --no-index
htseq-count --version
```
### htseq-count command and batch job submission script

I think it's better to use `-i gene_name` instead of the default `gene_id`.
```bash
#!/bin/bash
#SBATCH --account=def-fdick
#SBATCH --mem=32G
#SBATCH --job-name=count_sam
#SBATCH --output=%x.out
#SBATCH --time=30:00
htseq-count -i gene_name -f sam -s no -q ./sam/NTC_D1Aligned.out.sam ./sam/NTC_D2Aligned.out.sam ./sam/NTC_D3Aligned.out.sam ./sam/NTC_G1Aligned.out.sam ./sam/NTC_G2Aligned.out.sam ./sam/NTC_G3Aligned.out.sam ./sam/DKO_D1Aligned.out.sam ./sam/DKO_D2Aligned.out.sam ./sam/DKO_D3Aligned.out.sam ./sam/DKO_G1Aligned.out.sam ./sam/DKO_G2Aligned.out.sam ./sam/DKO_G3Aligned.out.sam ./gencode.vM1.annotation.gtf > count_sam.csv
```
