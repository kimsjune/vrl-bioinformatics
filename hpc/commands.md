## Week 2  

### to open and edit the bashrc file

```{bash}
nano ~/.bashrc
```

### what to write in the file

```{bash}
module load StdEnv/2020 star
export WS="/scratch/$USER/workshop"
export ACC="def-YOUR_PI_NAME"
```  

`control+x`, `y` and `enter` to save and exit.  

### download genome fasta and annotation

```{bash}
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.fa.gz  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz
```

### gunzip the two files

```{bash}
gunzip -k *
```

### index job script

```{bash}
#!/bin/bash
#SBATCH --account=def-fdick
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=16G
#SBATCH --job-name=STAR_index
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /scratch/skim823/workshop/stargenome --genomeFastaFiles /scratch/skim823/workshop/mm10.fa  --sjdbGTFfile /scratch/skim823/workshop/gencode.vM1.annotation.gtf --sjdbOverhang 74
```

### submit the index job

```{bash}
sbatch --account=$ACC index.sh
```