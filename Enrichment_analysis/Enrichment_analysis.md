cd "/home/zhaoling/Documents/Enrichment_analysis"


4. KEGG analysis

Reference species: mmu

1.1.int_ref_Genes508_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu_KEGG_1.tsv


2. 2.int_ref_Genes158_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu_KEGG_2.tsv

3. 3.int_ref_Genes133_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu_KEGG_3.tsv


5. GO analysis

Reference species: mmu

1.1.int_ref_Genes508_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_1.tsv


2. 2.int_ref_Genes158_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_2.tsv


3. 3.int_ref_Genes133_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_3.tsv


6. KEGG analysis

Reference species: ocu

1.1.int_ref_Genes508_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu_KEGG_1.tsv


2. 2.int_ref_Genes158_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu_KEGG_2.tsv

3. 3.int_ref_Genes133_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu_KEGG_3.tsv


7. GO analysis

Reference species: ocu

1.1.int_ref_Genes508_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_1.tsv

2. 2.int_ref_Genes158_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_2.tsv

3. 3.int_ref_Genes133_adj.fa

  sudo docker run --rm \
  -v "/home/zhaoling/Documents/Enrichment_analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_3.tsv



8. Check that the analyses produced output

After each run, check the output folder as follows:

ls -lh KOBAS_output_mmu_KEGG_1
ls -lh KOBAS_output_mmu_KEGG_1
ls -lh KOBAS_output_mmu_GO_1
ls -lh KOBAS_output_ocu_KEGG_1
ls -lh KOBAS_output_ocu_GO_1

To list all generated files:

find . -type f -name "*KOBAS*" -print
find "/home/zhaoling/Documents/Enrichment_analysis" -maxdepth 2 -type f
