#!/bin/bash

# Verifica se o input foi passado como argumento
if [ "$#" -ne 1 ]; then
    echo "Nenhum arquivo especificado, processando todos os arquivos em data/reads/"
    # Define o diretório das reads
    READS_DIR="data/reads/"
else
    # Define o arquivo ou diretório passado como argumento
    INPUT_PATH="$1"
    
    if [ -f "$INPUT_PATH" ]; then
        # Se for um arquivo, processa ele diretamente
        FILES=("$INPUT_PATH")
    elif [ -d "$INPUT_PATH" ]; then
        # Se for um diretório, usa esse diretório para buscar os arquivos .fastq
        FILES=("$INPUT_PATH"/*.fastq)
    else
        # Se o caminho não for válido (nem arquivo nem diretório)
        echo "Erro: O caminho fornecido não é válido."
        exit 1
    fi
fi

# Cria o diretório de saída se não existir
filtered_dir_name="data/filtered_reads"
mkdir -p "$filtered_dir_name"

# Inicializa variáveis com valores padrão
trim_qual_right=20
trim_qual_window=5
trim_qual_step=1
trim_qual_type=min
min_len=30
out_format=3
config_file="prinseq.config"

# Carrega as configurações do arquivo de configuração, se existir
if [ -f "$config_file" ]; then
    while IFS='=' read -r key value; do
        case "$key" in
            trim_qual_right) trim_qual_right="$value" ;;
            trim_qual_window) trim_qual_window="$value" ;;
            trim_qual_step) trim_qual_step="$value" ;;
            trim_qual_type) trim_qual_type="$value" ;;
            min_len) min_len="$value" ;;
            out_format) out_format="$value" ;;
        esac
    done < "$config_file"
fi

# Processa os arquivos fornecidos
for file in "${FILES[@]}"; do
    if [ -f "$file" ]; then
        # Extrai o nome do arquivo sem a extensão .fastq
        base_name=$(basename "$file" .fastq)
        
        # Define os diretórios para os arquivos 'good' e 'bad'
        output_dir="$filtered_dir_name/$base_name"
        out_good_path="$output_dir/good"
        out_good="$out_good_path/${base_name}_filtered_Q${trim_qual_right}_${min_len}"
        out_bad_path="$output_dir/bad"
        out_bad="$out_bad_path/${base_name}_bad"
        
        # Cria os diretórios de saída
        mkdir -p "$out_good_path" "$out_bad_path"
        
        # Executa o comando prinseq-lite
        tools/prinseq-lite -verbose -fastq "$file" \
            -out_good "$out_good" \
            -out_bad "$out_bad" \
            -trim_qual_right "$trim_qual_right" \
            -trim_qual_window "$trim_qual_window" \
            -trim_qual_step "$trim_qual_step" \
            -trim_qual_type "$trim_qual_type" \
            -min_len "$min_len" \
            -out_format "$out_format"
    else
        echo "Arquivo inválido ou inexistente: $file"
    fi
done
