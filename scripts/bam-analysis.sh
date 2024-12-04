#!/bin/bash

# Diretório contendo os arquivos BAM
INPUT_DIR="data/sorted_reads"
# Diretório para armazenar os resultados
OUTPUT_DIR="analysis/flag_depth_stats"

# Cria o diretório de saída, se não existir
mkdir -p "$OUTPUT_DIR"

# Loop para percorrer todos os arquivos BAM no diretório de entrada
for bam_file in "$INPUT_DIR"/*.bam; do
    # Extrai o nome base do arquivo BAM sem a extensão
    base_name=$(basename "$bam_file" .bam)
    
    # Executa o comando flagstat e redireciona para o arquivo de saída
    tools/samtools/samtools flagstat "$bam_file" >> "$OUTPUT_DIR/${base_name}.flag"
    
    # Executa o comando depth e redireciona para o arquivo de log de profundidade
    tools/samtools/samtools depth "$bam_file" >> "$OUTPUT_DIR/log_depth_${base_name}.txt"
done
