# SRA-TOOLS (obtener información)
```r
PASO 1 : obtener el link de descarga
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
PASO 2: descargar, crear el stk.tar.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz
PASO 3 : darle permiso al archivo
Chmod 777 stk.tar.gz
PASO 4:Añadimos la ruta al directorio y buscamos el archivo bin
export PATH=$PATH:$PWD/sratoolkit.3.1.1-ubuntu64/bin
ls;
paso 5: ubicar los archivos en directorio actual
mv ERR*/* .
paso 6: todos los archivos en SRA a FASTQ, separando las lecturas paired-end en dos archivos
fasterq-dump --split-files *.sra
paso 7: comprime los archvios fastq
gzip *fastq 
paso 8: 
fastqc *
```
#Código 2: Descargar archivo .sra y dividirlo en dos archivos .fastq

prefetch -h 
prefetch --max-size 50G --option-file sra_accessions_1.txt  #;
rm -r ERR12389866/ ERR12543675/  ##Borrar las carpetas pero con cuidado
fasterq-dump --split-files *.sra #separar en forward y reverse;
gzip *fastq ;
fastqc *


# Ensamblaje por mapeo

```r
##Tratar de que todos los archivos descargados Fastq, referencia, sam y los indexados esten dentro de una misma carpeta
#0# descargar de NCBI
prefetch --max-size 50G --option-file accessions_mpox.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra ;
gzip *fastq ;
fastqc * ;

#1# indexar el genoma de referencia#
bwa index reference.fasta ;

##Ejecutar todo el paso 2 y 3 en una sola corrida porqeu es un bucle########################
#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 4 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 4 ${prefix}.bam ;
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls ;
########################################

#4# extraer genoams consenso #
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 
ls ;

#5# annotation #
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

conda deactivate ;
cp */*.ffn ffn/ ; 
ls ; 
```

# Instalar Prokka
```r
#creamos el ambiente
conda create -n prokka-env
#activamos
conda activate
#instalr prokka dentro del ambiente
sudo apt install prokka
#ver si esta prokka
prokka -h
```
