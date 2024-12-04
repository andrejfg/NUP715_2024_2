#!/bin/bash

# Arquivo contendo os links
LINKS_FILE="read-links.txt"

# Diretório de destino
DEST_DIR="data/reads/"

# Garantir que o diretório de destino exista
mkdir -p "$DEST_DIR"

# Verificar se o arquivo de links existe
if [[ ! -f "$LINKS_FILE" ]]; then
    echo "Arquivo $LINKS_FILE não encontrado!"
    exit 1
fi

# Loop para processar cada link no arquivo
while IFS= read -r link; do
    # Ignorar linhas vazias
    [[ -z "$link" ]] && continue

    # Extrair o número de identificação do link
    id=$(echo "$link" | grep -oP '(?<=acc=)[A-Za-z0-9]+')

    # Nome do arquivo de destino
    output_file="${DEST_DIR}${id}.fastq.gz"

    # Fazer o download
    echo "Baixando $link..."
    curl -L -o "$output_file" "$link"

    # Verificar se o download foi bem-sucedido
    if [[ -f "$output_file" ]]; then
        echo "Download concluído: $output_file"

        # Descompactar o arquivo
        echo "Descompactando $output_file..."
        gunzip -f "$output_file"
    else
        echo "Erro ao baixar $link"
    fi
done < "$LINKS_FILE"

echo "Processo concluído!"
