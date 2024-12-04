#!/bin/bash

# Diretório de referência e de saída
REFERENCE_DIR="data/reference"
READS_DIR="data/reads"
OUTPUT_DIR="data/mapped_reads"

# Verifica se o diretório de saída existe, se não, cria
if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Se parâmetros forem passados para o script, usa-os como arquivos de leitura
if [ $# -gt 0 ]; then
  READ_FILES=("$@")
else
  # Caso contrário, procura arquivos .fastq na pasta data/reads
  READ_FILES=("$READS_DIR"/*.fastq)
fi

# Índice de referência para o Bowtie2
REFERENCE_INDEX="$REFERENCE_DIR/S288C_reference_index"

# Loop para realizar o alinhamento para cada arquivo de leitura
for READ_FILE in "${READ_FILES[@]}"; do
  # Extrai o nome do arquivo sem o caminho e a extensão
  BASENAME=$(basename "$READ_FILE" .fastq)
  
  # Define o nome do arquivo de saída SAM
  OUTPUT_FILE_SAM="$OUTPUT_DIR/${BASENAME}_bowtie2.sam"
  OUTPUT_FILE_BAM="$OUTPUT_DIR/${BASENAME}_bowtie2.bam"
  
  # Executa o Bowtie2
  echo "Alinhando $READ_FILE para $OUTPUT_FILE_SAM"
  bowtie2 --very-sensitive -x "$REFERENCE_INDEX" -U "$READ_FILE" -S "$OUTPUT_FILE_SAM" -p 4
  samtools view -b -S $OUTPUT_FILE_SAM > $OUTPUT_FILE_BAM
  
  # Verifica se o comando foi bem-sucedido
  if [ $? -eq 0 ]; then
    echo "Alinhamento de $READ_FILE concluído com sucesso!"
  else
    echo "Erro no alinhamento de $READ_FILE"
  fi
done
