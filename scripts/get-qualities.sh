#!/bin/bash

# Arquivo contendo os links
READS_DIR="data/reads/"

# Diretório de destino
DEST_DIR="analysis/qualities"

# Verifica se foi passado um arquivo como parâmetro
if [[ -n "$1" && -f "$1" ]]; then
    # Se um arquivo for passado como parâmetro, usa esse arquivo
    FILES=("$1")
else
    # Caso contrário, percorre todos os arquivos .fastq no diretório
    FILES=("$READS_DIR"*.fastq)
fi

# Cria o diretório de destino, se ainda não existir
mkdir -p "$DEST_DIR"

# Loop sobre todos os arquivos especificados
for file in "${FILES[@]}"; do
    if [[ -f "$file" ]]; then
        # Extrai o nome do arquivo sem extensão
        base_name=$(basename "$file" .fastq)
        
        # Define o diretório de saída específico para cada arquivo
        output_dir="$DEST_DIR/$base_name"
        
        # Cria o diretório de saída
        mkdir -p "$output_dir"
        
        # Executa o FastQC com as opções fornecidas
        tools/fastqc/fastqc "$file" --extract --delete -o "$output_dir"
    else
        echo "Nenhum arquivo .fastq encontrado no diretório $READS_DIR ou o arquivo especificado não é válido."
    fi
done
