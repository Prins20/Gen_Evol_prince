# GENOMICA EVOLUTIVA
Coleccion de codigos
# Descargar Genomas - (sratoolkit ) - CODIGO 1
```r
PASO 1 : obtener el link de descarga
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

PASO 2: descargar, crear el archivo de salida stk.tar.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz

PASO 3 : darle permiso al archivo
chmod 777 stk.tar.gz

PASO 4:Añadimos la ruta al directorio y el archivo lo llevamos al bin
export PATH=$PATH:$PWD/sratoolkit.3.1.1-ubuntu64/bin

PASO 5 : Verificar que estea el programa
fastq-dump --version
prefetch --version


 # OTRO :Forma de descarga sra-toolkit
 paso 1: Instalar, programa
 sudo apt install sra-toolkit

 paso 2: Verificar que estea el programa
 fastq-dump --version
 prefetch --version
```


# Descargar archivo .sra y dividirlo en dos archivos .fastq - CODIGO 2
```r
paso 1: Llamar al programa
prefetch -h
paso 2: Descargar datos SRA con un tamaño maximo de 50G
prefetch --max-size 50G --option-file sra_accessions_1.txt
paso 3: Eliminar las carpetas
rm -r ERR12389866/ ERR12543675/
paso 4: Convertir lo sarchivos SRA a FASTQ, separar el forward y reverse
fasterq-dump --split-files *.sra
paso 5: Comprimir los FASTQ generados
gzip *fastq
paso 6: realizar control de calidad de  FASTQ 
fastqc *fastq.gz
```

# Ensamblaje por mapeo - CODIGO 3

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
ls 
```


#  Intalar y correr PROKKA-CODIGO 4
```r
# instalacion de Prokka #
conda create -n prokka_env ;
conda activate prokka_env ;
conda install -c conda-forge -c biocondaconda install conda-forge::r-base prokka ;

# analisis en Prokka
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

# instalacion de artemis #
conda create -n art
conda activate art
conda install bioconda::artemis

conda install conda-forge::r-base

# desactivar conda 
conda deactivate ;
cp */*.ffn ffn/ ; 
ls ; 

# Instalar Prokka

#creamos el ambiente
conda create -n prokka-env
#activamos
conda activate
#instalr prokka dentro del ambiente
sudo apt install prokka
#ver si esta prokka
prokka -h
```
# ANOTACIÓN- CODIGO 5

```r
# si desea correr el comando desde un archico ".sh" debe añardir la "shebang" al inicio #
# 5.1: crear un archivo con extension *.sh, ejemplo "comando_1.sh"#
# 5.2: pegar el siguiente contenido o a través de cat o a través de "nano" #

#!/usr/bin
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;conda install bioconda::blast

# 5.3: dar permiso al archivo generado
chmod 777 comando_1.sh

# 5.4: correr el programa
./comando_1.sh
```
# VERIFIC AR LA CALIDA DDE ANOTACIÓN- CODIGO 6
```r
# 6.1: para obervar los headers de cada contig en todos los archivos *.fa
grep ">" *.fa

# 6.2: para observar todo el contenido de todos los archivos *.fa
cat *.fa

# 6.3: para obervar las 10 primeras lineas de cada archivo *.fa
cat OQ603638.fa | head -n 10
cat OQ603651.fa | head -n 10
cat SRR30716253.fa | head -n 10
cat SRR30716253.fa | head -nconda install bioconda::blast 10

# 6.4: para observar los resultados de la anotación en ARTEMIS, debe contar con el archivo *fa original y el  archivo *.gff. Cargue primero el genoma y luego el archivo de anotación
conda activate art
art
file >> open file manager >> cargar el genoma en extension *fa
file >> read and entry >> cargar el archivo *gff
explorar

# 6.5: identifique las regiones inferidas por el programa (CDS), identifique si se identificó la identidad de esas regiones o si algunas aparecen como "hypothetical"
```


# ensamblaje con nanopore (PROGRAMAS A INTALAR) - CODIGO 7

```r
# 7.1 : instalar los programas
# 7.1.1 : NanoPlot : calidad de secuencias Nanopore
conda activate base
pip install NanoPlot
pip install NanoPlot --upgrade
# 7.1.2 : Nanofilt : Filtrado por calidad de lecturas Nanopore
pip install nanofilt
pip install nanofilt --upgrade
# 7.1.3 : Flye: de-novo assembly
git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install


 7.1.4 : Minimap2 : polishing (parte 1)
git clone https://github.com/lh3/minimap2
cd minimap2 && make

# 7.1.5 : Racon : polishing (parte 2)# correr linea por linea
git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
 sudo apt  install cmake #solo si te lo piden#
make
cd build/bin/ 
export PATH=$PATH:$HOME/bin
cp racon $HOME/bin
chmod +x $HOME/bin/racon
# 7.1.6 : Requerimientos de MEDAKA (Pyabpoa, bcftools, samtools (v1.11), minimap2)

pip install pyabpoa
sudo apt install bcftools
conda install -c bioconda samtools==1.11

# 7.1.7 : MEDAKA, secuencias consenso (si MEDAKA no funciona correctamente, instala los programas requeridos)

conda install -c conda-forge –c bioconda medaka

or

pip install medaka

```
#ENSAMBLAJE NANOPORE (PIPELINES)-CODIGO 8
material de apoyo > https://denbi-nanopore-training-course.readthedocs.io/en/stable/index.html
```r
# 8.1: descargar la informacion (códigos SRR17110067 y SRR17110070)
mkdir sra_files ;
prefetch --max-size 50G --option-file accessions.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra 
gzip *.fastq ;
mkdir sra_files ;
mv *.sra sra_files/ ;

# 8.2: inspeccionar las longitudes de los reads ##
zcat SRR17110067.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20
zcat SRR17110070.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20

# 8.3: NanoPlot
NanoPlot -t 2 -o SRR17110067_QC --fastq SRR17110067.fastq.gz
NanoPlot -t 2 -o SRR17110070_QC --fastq SRR17110070.fastq.gz

# 8.4: NanoFilt
gunzip -c SRR17110067.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110067.trim.fastq.gz ;
gunzip -c SRR17110070.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110070.trim.fastq.gz ;
ls -lh ;

# 8.5: Flye
flye -o SRR17110067.genoma --nano-raw SRR17110067.trim.fastq.gz --threads 4 ;
flye -o SRR17110070.genoma --nano-raw SRR17110070.trim.fastq.gz --threads 4 ;
ls -lh ;

# 8.6 : Minimap2 + Racon (Polishing)
minimap2 -x ava-ont -t 4 SRR17110067.genoma/assembly.fasta SRR17110067.trim.fastq.gz > overlaps1.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps1.paf SRR17110067.genoma/assembly.fasta > SRR17110067.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.genoma/assembly.fasta SRR17110070.trim.fastq.gz > overlaps2.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps2.paf SRR17110070.genoma/assembly.fasta > SRR17110070.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110067.racon1.fasta SRR17110067.trim.fastq.gz > overlaps3.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps3.paf SRR17110067.racon1.fasta > SRR17110067.racon2.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.racon1.fasta SRR17110070.trim.fastq.gz > overlaps4.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps4.paf SRR17110070.racon1.fasta > SRR17110070.racon2.fasta ;

# 8.7 : Medaka (consensus)
medaka_consensus -i SRR17110070.trim.fastq.gz -d SRR17110070.racon2.fasta -o medaka_SRR17110070 -t 4 ;
medaka_consensus -i SRR17110067.trim.fastq.gz -d SRR17110067.racon2.fasta -o medaka_SRR17110067 -t 4 ;

# 8.8 : QUAST
quast.py -o quast_results -m 0 consensus.fasta

# 8.9 : Bandage

```
# BLAST - CODIGO 9
```r
# 9.1 : instalacion a traves de CONDA
conda install bioconda::blast
or
conda install -c conda-forge -c bioconda -c defaults blast

# 9.2 : http://www.mgc.ac.cn/VFs/
# 9.3 : Default webpage accessible to all users worldwide
# 9.4 : Download
# 9.5 : DNA sequences of full dataset
# 9.6 : Protein sequences of full dataset
# 9.6 :
gzip -d VFDB_setB_nt.fas.gz 
gzip -d VFDB_setB_pro.fas.gz

# 9.7 : run BLAST+
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl ;
blastn -db VFDB_setB_nt.fas -query GCA_001183825.1.fasta -perc_identity 90 -outfmt 6 -num_threads 4 > blast.csv ;
head blast.csv ;
cat blast.csv ;

# 9.8 : headers
sed '1i query.acc.ver subject.acc.ver perc.identity alignment.length mismatches gap.opens q.start q.end s.start s.end evalue bit.score' blast.csv | tr " " "\t" > blast.2.csv

# 9.9 : revisar resultados
head blast.2.csv
cat blast.2.csv

# 9.10 : instalar R
conda install -c conda-forge -c bioconda -c defaults r-base

# 9.11 : analizar los datos en R
#leer la data resultante de blast#
data <- read.csv("blast.2.csv", sep="\t", header=TRUE)
#conocer el número de filas y columnas de la tabla resultante#
dim(data)
#conocer las filas asignadas a una columna determinada#
length(data$subject.acc.ver)
#conocer el número de elementos únicos de esa columna#
length(unique(data$subject.acc.ver))
length(unique(data$query.acc.ver))
#conocer estadísticos básicos en un solo paso#
summary(data$query.acc.ver)
summary(data$alignment.length)
#obtener un boxplot de los porcentajes de identidad#
boxplot(data$perc.identity)
boxplot(data$perc.identity, xlab="genoma", ylab="% identidad")
summary(data$perc.identity)
data.frame(names(data))
#obtenet un plot longitud de alineamiento vs %identidad#
plot(data$alignment.length, data$perc.identity, xlab="length", ylab="% identity", main="BLASTn VFDB vs Chlamydia", pch=16, col="blue", cex=2)

# 9.12 : instalar bedtools desde conda para extraer las regiones "blasteadas"

conda install conda install -c conda-forge -c bioconda -c defaults bedtools


bedtools getfasta -fi  GCA_001183825.1.fasta -bed extract.txt -fo virulence.fasta

# 9.13 : Generar un archivo "bed" con R

#leer nuevamente en R#
data <- read.csv("blast.2.csv", sep="\t", header=TRUE)
#crear un objeto con las columnas requeridas para un archivo "bed"#
head(data)
seq <- data.frame(genome=data$query.acc.ver, start=data$q.start, end=data$q.end)
head(seq)
#generar el archivo "bed"#
write.table(seq, "extract.txt", sep="\t", row.names = F, col.names =F, quot=F)

# 9.14 : emplear "bedtools getfasta", identificar los argumentos # 

bedtools getfasta -fi  GCA_001183825.1.fasta -bed extract.txt -fo virulence.fasta

# 9.15 : extraer información de los headers de VFDB 
grep ">" VFDB_setB_nt.fas | sed -e 's/]\ \[/*/g' | sed -e 's/]//g' | sed -e 's/\ \[/*/g'| sed -e 's/)\ /*/g' | sed -e 's/*(/*/g' | head -n 10 > headers.txt

# 9.16 : traducir las secuencias con virtual ribosome #
https://services.healthtech.dtu.dk/services/VirtualRibosome-2.0/

```
# ORTHO. ANI -CODIGO 10
```r
# 10.1: instalar NCBI-DATASETS
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

# 10.2: emplear el archivo "accessions.txt"

# 10.3: descargar los genomas con la lista sugerida. Emplear el codigo "command_ncbidatasets.sh" disponible en https://github.com/Vjimenez-vasquez/NCBI-DATASETS

./command_ncbidatasets.sh accessions.txt 

# 10.4: guardar los genomas en una nueva carpeta y comprimirla

# 10.5: ingresar a la pagina de OAT para correr el algoritmo ORTHO-ANI
https://www.ezbiocloud.net/tools/orthoani

# 10.6: descargar BLAST+
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
descargar "ncbi-blast-2.16.0+-win64.exe" para instarlar en el sistema
o
descargar "ncbi-blast-2.16.0+-x64-win64.tar.gz" para descargar los ejecutables directamente

# 10.7: entender el algoritmo
https://help.ezbiocloud.net/orthoani-genomic-similarity/
https://pypi.org/project/orthoani/
```
# PANGENOME ANALYSIS -CODIGO 11
```r
conda install -c conda-forge -c defaults -c bioconda roary
conda install -c conda-forge -c defaults -c bioconda snp-sites
conda install -c conda-forge -c defaults -c bioconda raxml
conda install -c conda-forge -c defaults -c bioconda figtree

# 11.1: annotation (PROKKA)
conda activate prokka_env
mkdir -p annotation ;
mkdir -p ffn ;
for r1 in *fasta
do
prefix=$(basename $r1 .fasta)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Bacteria ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;
conda deactivate ;
cp */*.ffn ffn/ ; 
ls ;

# 11.2: inferring clusters, core genes and accesory genes (ROARY)
# https://github.com/sanger-pathogens/Roary #
# roary -p 4 -f roary_output -g 200000 -z -r -e -n -v -cd 80 -i 90 annotation/*.gff ; #
roary -p 4 -f roary_output -g 200000 -r -e -n -v -cd 80 -i 90 annotation/*.gff ;
cp roary_output/core_gene_alignment.aln . ;
ls -lh ; 

# 11.3: SNPs alignment (SNP-SITES)
snp-sites -m -o snp1.phy core_gene_alignment.aln ; 
snp-sites -m -c -o snp2.phy core_gene_alignment.aln ; 
ls -lh ;

# 11.4: phylogeny (RAXML)
raxmlHPC-PTHREADS -p 1254512 -m GTRCAT -s snp2.phy -n nwk -# 20 -T 4 ;
mv RAxML_bestTree.nwk raw_tree.nwk ;
rm RAxML_* ;
mkdir phylogeny ;
mv snp1.phy snp2.phy snp2.phy.reduced raw_tree.nwk core_gene_alignment.aln phylogeny/ ;

# 11.5: cargar el programa "pangenome_command_2.R" en R o R-Studio
#ingresar la ruta correcta en cada caso (donde se encuentran los archivos "gene_presence_absence.csv" y "metadata_1.tsv")
setwd("")
dir()

#ingresar la ruta correcta hasta donde se encuentra el archivo pangenome_command_2.R
source("../../pangenome_command_2.R")

pangenome <- pres_abs(metadata = "metadata_1.tsv", roary_output = "gene_presence_absence.csv", last_column = "3", output = "out_5.tsv")
head(pangenome)
class(pangenome)
pangenome[,1:10]

# 11.6 : visualizacion en microreact
https://microreact.org/
cargar el arbol enraizado (formato .nwk) y la metadata final (out_5.tsv)
```


































