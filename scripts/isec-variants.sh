#!/bin/bash

# Caminho da pasta de entrada com os arquivos .vcf.gz
INPUT_DIR="analysis/variants"

# Caminho da pasta de saída onde os arquivos .bcf serão salvos
OUTPUT_DIR="analysis/variants_isec/"

# Criação da pasta de saída, caso não exista
mkdir -p "$OUTPUT_DIR"

# Monta a lista de arquivos .vcf.gz no diretório de entrada
VCF_FILES=$(ls "$INPUT_DIR"/*.vcf.gz)

# Chama o comando bcftools isec para processar todos os arquivos em uma única execução
tools/bcftools/bin/bcftools isec $VCF_FILES -p "$OUTPUT_DIR" -n =1
