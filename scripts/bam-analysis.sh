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

    OUTPUT_LOG_DEPTH="$OUTPUT_DIR/log_depth_${base_name}.txt"    
    
    # Executa o comando depth e redireciona para o arquivo de log de profundidade
    tools/samtools/samtools depth "$bam_file" >> $OUTPUT_LOG_DEPTH

    # Conta as linhas do arquivo de log e faz o cálculo
    line_count=$(wc -l < "$OUTPUT_LOG_DEPTH")

    # Realiza o cálculo (Valor do wc -l / 12157105) * 100
    result=$(echo "scale=4; ($line_count / 12157105) * 100" | bc)

    # Cria o arquivo de resultado e exporta o cálculo de cobertura por extensão
    OUTPUT_RESULT="$OUTPUT_DIR/${base_name}_extension_coverage.txt"
    echo "$result" > "$OUTPUT_RESULT"

    # Cria o arquivo de resultado e exporta o cálculo de cobertura em profundidade
    awk '{total += $3; count++} END { print total/count}' "$OUTPUT_LOG_DEPTH" > "$OUTPUT_DIR/${base_name}_depth_coverage.txt"

done
