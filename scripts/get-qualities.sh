#!/bin/bash

# Arquivo contendo os links
READS_DIR="data/reads/"

# Diretório de destino
DEST_DIR="analysis/qualities"

# Cria o diretório de destino, se ainda não existir
mkdir -p "$DEST_DIR"

# Loop sobre todos os arquivos .fastq na pasta de reads
for file in "$READS_DIR"*.fastq; do
    if [[ -f "$file" ]]; then
        # Extrai o nome do arquivo sem extensão
        base_name=$(basename "$file" .fastq)
        
        # Define o diretório de saída específico para cada arquivo
        output_dir="$DEST_DIR/$base_name"
        
        # Cria o diretório de saída
        mkdir -p "$output_dir"
        
        # Executa o FastQC com as opções fornecidas
        tools/fastqc "$file" --extract --delete -o "$output_dir"
    else
        echo "Nenhum arquivo .fastq encontrado no diretório $READS_DIR."
    fi
done
