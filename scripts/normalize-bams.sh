#!/bin/bash

# Caminho da pasta que contém os arquivos BAM
INPUT_DIR="data/mapped_reads"
OUTPUT_DIR="data/sorted_reads"

# Cria o diretório de saída, se não existir
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Itera sobre todos os arquivos .bam na pasta
for bam_file in "$INPUT_DIR"/*.bam; do
  # Verifica se o arquivo existe
  if [[ -f "$bam_file" ]]; then
    # Define o nome do arquivo de saída para o sort com o sufixo _sorted
    sorted_bam="$(basename "${bam_file%.bam}_sorted.bam")"
    
    # Executa o comando samtools sort
    echo "Executando samtools sort em $bam_file"
    tools/samtools/samtools sort -T "${bam_file%.bam}" -o "$OUTPUT_DIR/$sorted_bam" "$bam_file"
    
    # Executa o comando samtools index
    echo "Executando samtools index em $OUTPUT_DIR/$sorted_bam"
    tools/samtools/samtools index "$OUTPUT_DIR/$sorted_bam"
  fi
done
