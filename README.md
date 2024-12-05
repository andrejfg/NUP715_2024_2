### NUP715_2024_2

Este repositório contém todos os scripts necessários para realizar um estudo de variantes genéticas relacionadas à tolerância ao etanol em *Saccharomyces cerevisiae*, com foco na produção de biocombustíveis. A estrutura do projeto está organizada em diretórios e scripts específicos para facilitar a execução do pipeline.

---

### Estrutura do Repositório

- **`analysis/`**  
  Contém arquivos resultantes das análises geradas pelos pipelines.

- **`data/`**  
  Diretório para armazenar arquivos de entrada como **FASTQ**, **SAM**, e **BAM**.

- **`scripts/`**  
  Scripts individuais que compõem as etapas do pipeline.

- **`tools/`**  
  Ferramentas e dependências necessárias para executar os pipelines.

- **`execute-all.sh`**  
  Script principal que executa todas as etapas do pipeline em sequência.

- **`prinseq.config`**  
  Arquivo de configuração com parâmetros para filtros do PRINSEQ.

- **`read-links.txt`**  
  Lista de URLs para download das reads que serão analisadas.

---

### Como Executar

1. Certifique-se de ter permissão de superusuário (root).  
   O script `prepare-env.sh` instala dependências, exigindo privilégios administrativos.
   ```bash
   sudo bash prepare-env.sh
   ```


2. Execute o pipeline:  
   ```bash
   bash execute-all.sh
   ```

3. Para ajustar os parâmetros dos filtros, edite o arquivo `prinseq.config`. Caso esteja ausente, serão usados os parâmetros padrões descritos abaixo:

   ```
   trim_qual_right=20
   trim_qual_window=5
   trim_qual_step=1
   trim_qual_type=min
   min_len=30
   out_format=3
   ```

---

### Descrição dos Scripts

- **`scripts/install-dependencies.sh`**  
  Instala todas as dependências necessárias para executar o pipeline (ex: PRINSEQ, SAMtools, bowtie2).

- **`scripts/get-reads.sh`**  
  Faz o download das reads especificadas no arquivo `read-links.txt`.

- **`scripts/get-qualities.sh`**  
  Avalia a qualidade das reads usando PRINSEQ ou ferramentas similares.

- **`scripts/filter-reads.sh`**  
  Aplica filtros nas reads conforme os parâmetros definidos em `prinseq.config`.

- **`scripts/index-reference.sh`**  
  Indexa a sequência de referência para alinhamento utilizando ferramentas como bowtie2.

- **`scripts/align_reads.sh`**  
  Realiza o alinhamento das reads contra a sequência de referência, gerando arquivos BAM.

- **`scripts/normalize-bams.sh`**  
  Normaliza os arquivos BAM para uniformizar coberturas.

- **`scripts/bam-analysis.sh`**  
  Executa análises nos arquivos BAM, incluindo métricas de cobertura e mapeamento.

- **`scripts/identify-variants.sh`**  
  Identifica variantes genéticas com base nos alinhamentos processados.

- **`scripts/isec-variants.sh`**  
  Filtra e intersecta variantes para análises específicas.

---

### Exemplo de Arquivos

- **`read-links.txt`**  
  URLs de reads para download:
  ```
  https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR14457781
  https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1570792
  ```

- **`prinseq.config`**  
  Configurações padrão:
  ```
  trim_qual_right=30
  trim_qual_window=5
  trim_qual_step=1
  trim_qual_type=min
  min_len=30
  out_format=3
  ```

---

### Contribuições

Contribuições são bem-vindas! Abra uma *issue* ou envie um *pull request* para sugestões de melhorias.
