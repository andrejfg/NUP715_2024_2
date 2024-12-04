#!/bin/bash

# Caminho da pasta de entrada com os arquivos .bam
INPUT_DIR="data/sorted_reads"

# Caminho da pasta de saída onde os arquivos .bcf serão salvos
OUTPUT_DIR="analysis/variants"

# Caminho do arquivo de referência
reference_file="data/reference/S288C_reference.fa"

# Criação da pasta de saída, caso não exista
mkdir -p "$OUTPUT_DIR"

# Percorrendo todos os arquivos .bam no diretório de entrada
for bam_file in "$INPUT_DIR"/*.bam; do
  # Extraindo o nome do arquivo base sem o caminho e sem a extensão .bam
  base_name=$(basename "$bam_file" .bam)

  # Caminho completo do arquivo de saída .bcf
  output_file="$OUTPUT_DIR/${base_name}_calls"

  # Executando o comando bcftools mpileup e bcftools call
  tools/bcftools/bin/bcftools mpileup -Ou -f "$reference_file" "$bam_file" | tools/bcftools/bin/bcftools call -mv -Ob -o "$output_file.bcf"
  tools/bcftools/bin/bcftools view "$output_file.bcf" | tools/bcftools/bin/vcfutils.pl varFilter > "$output_file.vcf"
  bgzip -c "$output_file.vcf" > "$output_file.vcf.gz"
  tools/bcftools/bin/bcftools index -f "$output_file.vcf.gz"
  
  echo "Processado: $bam_file -> $output_file"
done
