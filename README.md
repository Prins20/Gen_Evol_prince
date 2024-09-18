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
