# KEGG analysis
## Reference species: mmu
KEGG Organisms ==== https://www.kegg.jp/kegg/tables/br08606.html
### 1.int_ref_Genes508_adj
```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"

sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu
```

### 2.int_ref_Genes158_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa\
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu
```

### 3.int_ref_Genes133_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa\
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu
```

# GO analysis
## Reference species: mmu
### 1.int_ref_Genes508_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/GO_output_mmu
```

### 2.int_ref_Genes158_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/GO_output_mmu
```

---

# KEGG analysis
## Reference species: ocu
KEGG Organisms ==== https://www.kegg.jp/kegg/tables/br08606.html
### 1.int_ref_Genes508_adj
```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"

sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu
```
### 2.int_ref_Genes158_adj
```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"

sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu
```

### 3.int_ref_Genes133_adj
```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"

sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu
```


# GO analysis
## Reference species: mmu
### 1.int_ref_Genes508_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/GO_output_ocu
```

### 2.int_ref_Genes158_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/GO_output_ocu
```

### 3.int_ref_Genes133_adj

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/GO_output_ocu
```
