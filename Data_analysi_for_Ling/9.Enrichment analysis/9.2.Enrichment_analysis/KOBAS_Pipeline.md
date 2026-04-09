# Enrichment analysis using `**KOBAS**`(https://pmc.ncbi.nlm.nih.gov/articles/PMC3125809/)

![KOBAS workflow](KOBAS.png)

## 0. Prepare the input files
Make sure the FASTA protein files have been stored inside a specific folder, like so:

```text
/home/danilosantoro/Documents/Enrichment analysis
```

Expected files according to the Venn diagram:

```text
1.int_ref_Genes508_adj.fa
2.int_ref_Genes158_adj.fa
3.int_ref_Genes133_adj.fa
```

They must be protein FASTA files, for example:

```text
>gene1
MVEWCLPQDIDLEGVEFKS...
>gene2
MKAMPWNWTCLLSHLLMV...
```

### 1. Install Docker on Ubuntu
1.1. Open the terminal and run:

```bash
sudo apt update
sudo apt install -y docker.io
```

1.2. Start Docker session as follows:

```bash
sudo systemctl start docker
sudo systemctl enable docker
```

1.3. Optional, tough strongly recommended: allow Docker without always typing `sudo` over and over again

```bash
sudo usermod -aG docker $USER
newgrp docker
```

<mark>***Note: If that does not refresh permissions immediately, log out and log back in***</mark>.

### 2. Download the KOBAS Docker image

```bash
sudo docker pull agbase/kobas:3.0.3_3
```

If Docker works without `sudo`, then `sudo` can be removed from all later commands.

### 3. Move on to the working folder

```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"
```

Check the input files as follows:

```bash
ls -lh
```

### 4. Create separate output folders.
This is important to keep each analysis separated.

```bash
mkdir -p KOBAS_output_mmu_KEGG_1
mkdir -p KOBAS_output_mmu_KEGG_2
mkdir -p KOBAS_output_mmu_KEGG_3

mkdir -p KOBAS_output_mmu_GO_1
mkdir -p KOBAS_output_mmu_GO_2
mkdir -p KOBAS_output_mmu_GO_3

mkdir -p KOBAS_output_ocu_KEGG_1
mkdir -p KOBAS_output_ocu_KEGG_2
mkdir -p KOBAS_output_ocu_KEGG_3

mkdir -p KOBAS_output_ocu_GO_1
mkdir -p KOBAS_output_ocu_GO_2
mkdir -p KOBAS_output_ocu_GO_3
```

***Reference species codes***? **Take a look at** https://www.kegg.jp/kegg/tables/br08606.html
Looking for closely related species, we have the following KEGG organism codes:

- Mus musculus (house mouse) === > `mmu`

- Oryctolagus cuniculus (rabbit) === > `ocu`

### 4. KEGG analysis
Reference species: `mmu`

***1.1.int_ref_Genes508_adj.fa***

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
  -o /work-dir/KOBAS_output_mmu_KEGG_1
```

***2. 2.int_ref_Genes158_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu_KEGG_2
```

***3. 3.int_ref_Genes133_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d K \
  -o /work-dir/KOBAS_output_mmu_KEGG_3
```

### 5. GO analysis
Reference species: `mmu`

***1.1.int_ref_Genes508_adj.fa***

```bash
cd "/home/danilosantoro/Documents/Enrichment analysis"

sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_1
```

***2. 2.int_ref_Genes158_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_2
```
***3. 3.int_ref_Genes133_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s mmu \
  -d G \
  -o /work-dir/KOBAS_output_mmu_GO_3
```
### 6. KEGG analysis
Reference species: `ocu`

***1.1.int_ref_Genes508_adj.fa***
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
  -o /work-dir/KOBAS_output_ocu_KEGG_1
```

***2. 2.int_ref_Genes158_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu_KEGG_2
```

***3. 3.int_ref_Genes133_adj.fa***
```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d K \
  -o /work-dir/KOBAS_output_ocu_KEGG_3
```

### 7. GO analysis
Reference species: `ocu`

***1.1.int_ref_Genes508_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/1.int_ref_Genes508_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_1
```

***2. 2.int_ref_Genes158_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/2.int_ref_Genes158_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_2
```

***3. 3.int_ref_Genes133_adj.fa***

```bash
sudo docker run --rm \
  -v "/home/danilosantoro/Documents/Enrichment analysis:/work-dir" \
  agbase/kobas:3.0.3_3 \
  -j \
  -i /work-dir/3.int_ref_Genes133_adj.fa \
  -t fasta:pro \
  -s ocu \
  -d G \
  -o /work-dir/KOBAS_output_ocu_GO_3
```

### 8. Check that the analyses produced output
After each run, check the output folder as follows:

```bash
ls -lh KOBAS_output_mmu_KEGG_1
ls -lh KOBAS_output_mmu_KEGG_1
ls -lh KOBAS_output_mmu_GO_1
ls -lh KOBAS_output_ocu_KEGG_1
ls -lh KOBAS_output_ocu_GO_1
```

To list all generated files:

```bash
find . -type f -name "*KOBAS*" -print
find "/home/danilosantoro/Documents/Enrichment analysis" -maxdepth 2 -type f
```

