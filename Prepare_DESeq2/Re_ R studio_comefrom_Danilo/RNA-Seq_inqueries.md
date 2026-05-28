Hi, Ling!
I took a glance at your folder and found the following list:

(/scratch/project_2002674/conda_envs/lepus_jupiter) [dansanto@puhti-login15 RNA-Seq_PRJNA826339]$ ls -lht
total 86M
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 16:30 variants
-rw-rw---- 1 lingzhao project_2002674  987 Mar  5 16:30 10_variant_calling.sh
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 15:48 counts
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 15:45 ref
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 13:09 logs
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 12:12 qc
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 11:49 aln
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 10:05 star_index
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  5 09:39 trimmed
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  4 14:58 ref_src
-rw-rw---- 1 lingzhao project_2002674  236 Mar  3 15:20 config.sh
-rw-rw---- 1 lingzhao project_2002674 1.4K Mar  3 15:20 00_setup.sh
drwxrws--- 3 lingzhao project_2002674 4.0K Mar  1 19:54 qc_lepus_europaeus
drwxrws--- 2 lingzhao project_2002674 4.0K Mar  1 19:54 logs_lepus_europaeus
-rw-rw---- 1 lingzhao project_2002674  199 Mar  1 19:53 11_multiqc.sh
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 27 15:34 variants_lepus_europaeus
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 27 10:15 counts_lepus_europaeus
-rw-rw---- 1 lingzhao project_2002674  565 Feb 27 09:19 09_featurecounts.sh
-rw-rw---- 1 lingzhao project_2002674 1.6K Feb 26 22:17 temp.bed
-rw-rw---- 1 lingzhao project_2002674  39M Feb 26 21:35 genes.bed
-rw-rw---- 1 lingzhao project_2002674  751 Feb 26 14:45 08_post_alignment_qc.sh
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 26 13:52 aln_lepus_europaeus
-rw-rw---- 1 lingzhao project_2002674 1010 Feb 26 11:53 07_star_align.sh
drwx--S--- 2 lingzhao project_2002674 4.0K Feb 26 11:23 star_index_lepus_europaeus
-rwxrwx--- 1 lingzhao project_2002674  117 Feb 26 10:48 submit_interactive.sh
-rw-rw---- 1 lingzhao project_2002674  829 Feb 26 10:29 06_star_index.sh
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 24 16:06 trimmed_lepus_europaeus
-rw-rw---- 1 lingzhao project_2002674  970 Feb 24 15:08 05_fastq.sh
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 24 11:42 ref_lepus_europaeus
-rw-rw---- 1 lingzhao project_2002674  662 Feb 24 11:41 04_prepare_reference.sh
drwxrws--- 2 dansanto project_2002674 4.0K Feb 24 11:30 fastq
drwxrws--- 2 lingzhao project_2002674 4.0K Feb 23 11:12 tmp
drwxrws--- 2 dansanto project_2002674 4.0K Feb 13 10:56 sra
-rw-rw---- 1 lingzhao project_2002674   96 Feb 13 09:43 list.txt
-rw-rw---- 1 lingzhao project_2002674  48M May 10  2023 subread-2.0.6-Linux-x86_64.tar.gz
drwxr-s--- 6 lingzhao project_2002674 4.0K May  6  2023 subread-2.0.6-Linux-x86_64


Among these ones, I dove deeper into the folders called qc_lepus_europaeus and qc (I think that this ones corresponds for lepus timidus, right?). I strongly recommed you to create different repository, one for lepus europaeus and one for lepus timidus. This will avoid to mix up the data because the folders are not complitely clear to me. Within these ones, I found the multiqc_report.html only for lepus europaes and not for lepus timidus (see below). Can you try to produce it? Moreover, the same file should be produced also for the trimmed reads if it is possibile (not mandatory this last one).

# Lepus europaues

(/scratch/project_2002674/conda_envs/lepus_jupiter) [dansanto@puhti-login15 RNA-Seq_PRJNA826339]$ ls -lht qc_lepus_europaeus/
total 9.6M
drwxrwx--- 2 lingzhao project_2002674 4.0K Mar  1 19:54 multiqc_data
-rw-rw---- 1 lingzhao project_2002674 5.0M Mar  1 19:54 multiqc_report.html

# Lepus timidus

(/scratch/project_2002674/conda_envs/lepus_jupiter) [dansanto@puhti-login15 RNA-Seq_PRJNA826339]$ ls -lht qc/
total 4.6M
-rw-rw---- 1 lingzhao project_2002674  533 Mar  5 12:12 SRR18740842.flagstat.txt
-rw-rw---- 1 lingzhao project_2002674  179 Mar  5 12:12 star_mapping_summary.tsv


I also looked for the read counts files since we will use this one on DESeq2. I downloaded the content inside both counts_lepus_europaeus file and counts (I assumed that this one is related with Lepus timidus, right?) using the following commands:

### Download files Lepus europeaus

rsync -avP \
dansanto@puhti.csc.fi:/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/qc_lepus_europaeus/ \
"/home/danilosantoro/Documents/Hare_RNA-Seq pipeline/Ling_results_02-03-26/1.QC/"

rsync -avP \
dansanto@puhti.csc.fi:/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/counts_lepus_europaeus/ \
"/home/danilosantoro/Documents/Hare_RNA-Seq pipeline/Ling_results_02-03-26/9.featurecounts/"


### Download files Lepus timidus

rsync -avP \
dansanto@puhti.csc.fi:/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/qc/ \
"/home/danilosantoro/Documents/Hare_RNA-Seq pipeline/Ling_results_02-03-26/2.Lepus_timidus/1.QC/"

rsync -avP \
dansanto@puhti.csc.fi:/scratch/project_2002674/RNAseq_hares/scripts/RNA-Seq_trial/RNA-Seq_PRJNA826339/counts/ \
"/home/danilosantoro/Documents/Hare_RNA-Seq pipeline/Ling_results_02-03-26/2.Lepus_timidus/9.featurecounts/"

Inside these ones, I took the gene_counts_featureCounts.txt since this one is the correct file to use in DESeq2. After polishing the file with a specifi code in the Rstudio (I will give it to you), I observed that lepus europaes is ok, but the lepus timidus contains a list of genes starting from gene_1 till gene_23778. In this case, we lost the information about the gene_ID contained in the Lepus_timidus_annotation.gff file I forwarded to you by Jaakko and also in this case it is possibile (I'm not sure) that the same gene mapp back to different regions, which in turn will have more read counts for the same gene. I'm sending you the headings of both Lepus europaeus and Lepus timidus.

### Lepus europaeus

Geneid	SRR18740835	SRR18740836	SRR18740837	SRR18740838	SRR18740839	SRR18740840	SRR18740841	SRR18740842
LOC133763884	122	72	54	122	88	30	48	64
TRIM24	8430	4624	3753	12849	7607	1898	4458	6117
SVOPL	2	1	1	0	0	1	1	1
ATP6V0A4	0	0	4	0	13	3	24	0
TMEM213	0	0	0	0	0	0	0	0
KIAA1549	918	635	853	847	833	519	566	566
LOC133767786	8	10	0	14	21	18	29	30
ZC3HAV1L	1255	1326	411	1673	1001	1030	1395	936
ZC3HAV1	8372	1963	504	29305	1192	177	3239	1577
IFT56	2004	1904	1847	6569	1331	335	832	1193
UBN2	7220	5557	5339	11225	5504	2742	3855	4795
TRNAR-CCU	16	13	10	19	7	0	4	2
FMC1	1557	1353	823	2213	1384	500	1043	1927
LUC7L2	25671	15700	14987	50222	26020	6302	13488	24029
LOC133760428	4	0	2	13	0	2	0	3
KLRG2	125	66	72	272	48	62	49	89
CLEC2L	10	0	2	13	11	8	0	6
HIPK2	7005	4758	7471	6114	4893	3008	4575	2936
TBXAS1	44	1	24	86	26	0	26	3
PARP12	9485	4362	1291	22378	4117	990	6796	2439
LOC133749760	820	792	375	2611	316	298	1039	361
LOC133763922	10	11	7	13	3	3	4	0
KDM7A	4827	4257	4184	4480	2664	3036	3005	2306
LOC133763932	31	25	24	101	10	18	59	27
SLC37A3	4947	4817	3541	5756	3157	2705	3327	3429
LOC133768655	9	6	2	2	10	4	0	0
RAB19	19	5	3	75	14	0	14	35
MKRN1	11878	12224	8052	19366	7083	4241	5602	7748
DENND2A	2496	3293	295	3844	5546	313	2199	5030

### Lepus timidus

Geneid	SRR18740835	SRR18740836	SRR18740837	SRR18740838	SRR18740839	SRR18740840	SRR18740841	SRR18740842
gene_1	6491	3655	2870	9750	5662	1437	3412	4646
gene_2	0	0	0	0	0	0	0	0
gene_3	1000	844	852	2773	639	143	358	607
gene_4	2122	1422	1388	2996	1217	701	973	1182
gene_5	1227	1047	642	1715	1065	422	847	1515
gene_6	15165	8511	7899	25766	16379	3424	6689	12970
gene_7	10	0	2	10	8	8	0	6
gene_8	0	0	0	0	0	0	0	0
gene_9	40	1	13	77	22	0	26	3
gene_10	4	0	2	8	0	0	0	0
gene_11	18	5	3	71	14	0	13	35
gene_12	300	934	321	522	399	363	366	446
gene_13	11020	11029	8038	15354	8867	3761	3638	5076
gene_14	0	0	4	0	0	2	0	0
gene_15	2	1	0	0	0	0	0	0
gene_16	8	0	1	2	2	3	1	0
gene_17	5898	3230	1997	8378	4291	1605	2288	3993
gene_18	11	27	15	14	14	20	4	10
gene_19	6813	6441	3989	10614	5869	1292	2432	3790
gene_20	0	4	2	0	3	0	2	2
gene_21	2	0	1	0	0	0	0	3
gene_22	1	0	0	0	0	0	2	7
gene_23	26	1	2	2	0	0	6	20
gene_24	0	0	0	2	0	0	22	2
gene_25	0	0	0	0	0	0	0	0
gene_26	0	0	0	0	0	0	0	0
gene_27	0	0	0	0	0	0	0	0
gene_28	25	983	323	10	4	36	47	10
gene_29	0	0	0	0	0	0	0	0


As you can noticed, the Geneid for lepus europaes is good because you have the gene name, while for Lepus timidus we have number starting from gene_1 to gene_23778. I will send you the list so that you can look at them if you want.

Therefore, I want to ask you if you can take a look at the code you used to be sure what what we are observing is not an artefacts, especially for lepus timidus.

Anyways, thank you so much and you are in the right path!

Thank you for your work.

Dani