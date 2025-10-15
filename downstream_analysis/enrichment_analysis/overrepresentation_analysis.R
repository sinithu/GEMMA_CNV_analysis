#Data (gene lists) preparation and overrepresentation analysis

#Generating clean gene lists of protein-coding genes (excluding HLA genes, as these may lie in regions of the genome that are difficult to detect)
#Case DELs
case_dels <- read.table("/path/genes_biotype_DEL_cases.txt")
case_dels_genes <- unlist(strsplit(case_dels$V1, ","))
case_dels_biotypes <- unlist(strsplit(case_dels$V2, ","))
df_case_dels <- data.frame(
  Gene = case_dels_genes,
  Biotype = case_dels_biotypes)
df_case_dels_protco <- df_case_dels[df_case_dels$Biotype=="protein_coding",]
df_case_dels_protco_hla <- df_case_dels_protco[!grepl("^HLA", df_case_dels_protco$Gene), ]
genes_case_dels_protco_hla <- df_case_dels_protco_hla$Gene
case_genelist_DEL <- unique(genes_case_dels_protco_hla)
#Case DUPs
case_dups <- read.table("/path/genes_biotype_DUP_cases.txt")
case_dups_genes <- unlist(strsplit(case_dups$V1, ","))
case_dups_biotypes <- unlist(strsplit(case_dups$V2, ","))
df_case_dups <- data.frame(
  Gene = case_dups_genes,
  Biotype = case_dups_biotypes)
df_case_dups_protco <- df_case_dups[df_case_dups$Biotype=="protein_coding",]
df_case_dups_protco_hla <- df_case_dups_protco[!grepl("^HLA", df_case_dups_protco$Gene), ]
genes_case_dups_protco_hla <- df_case_dups_protco_hla$Gene
case_genelist_DUP <- unique(genes_case_dups_protco_hla)
#Control DELs
control_dels <- read.table("/path/genes_biotype_DEL_controls.txt")
control_dels_genes <- unlist(strsplit(control_dels$V1, ","))
control_dels_biotypes <- unlist(strsplit(control_dels$V2, ","))
df_control_dels <- data.frame(
  Gene = control_dels_genes,
  Biotype = control_dels_biotypes)
df_control_dels_protco <- df_control_dels[df_control_dels$Biotype=="protein_coding",]
df_control_dels_protco_hla <- df_control_dels_protco[!grepl("^HLA", df_control_dels_protco$Gene), ]
genes_control_dels_protco_hla <- df_control_dels_protco_hla$Gene
control_genelist_DEL <- unique(genes_control_dels_protco_hla)
#Control DUPs
control_dups <- read.table("/path/Desktop/genes_biotype_DUP_controls.txt")
control_dups_genes <- unlist(strsplit(control_dups$V1, ","))
control_dups_biotypes <- unlist(strsplit(control_dups$V2, ","))
df_control_dups <- data.frame(
  Gene = control_dups_genes,
  Biotype = control_dups_biotypes)
df_control_dups_protco <- df_control_dups[df_control_dups$Biotype=="protein_coding",]
df_control_dups_protco_hla <- df_control_dups_protco[!grepl("^HLA", df_control_dups_protco$Gene), ]
genes_control_dups_protco_hla <- df_control_dups_protco_hla$Gene
control_genelist_DUP <- unique(genes_control_dups_protco_hla)

#Generating ENTREZID gene lists for overrepresentation analysis
library("clusterProfiler")
library("org.Hs.eg.db")
#DELs
case_genes_DEL <- bitr(case_genelist_DEL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
control_genes_DEL <- bitr(control_genelist_DEL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
case_entrez_DEL <- case_genes_DEL$ENTREZID
control_entrez_DEL <- control_genes_DEL$ENTREZID
#DUPs
case_genes_DUP <- bitr(case_genelist_DUP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
control_genes_DUP <- bitr(control_genelist_DUP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
case_entrez_DUP <- case_genes_DUP$ENTREZID
control_entrez_DUP <- control_genes_DUP$ENTREZID

#Generating subsets of gene lists to include variants exclusive to cases and controls (i.e. excluding shared genes between cases and controls)
#DUPs
shared_genes_DUP <- intersect(case_genelist_DUP, control_genelist_DUP)
case_genelist_DUP_exclusive <- setdiff(case_genelist_DUP, shared_genes_DUP)
control_genelist_DUP_exclusive <- setdiff(control_genelist_DUP, shared_genes_DUP)
case_genes_DUP_exclusive <- bitr(case_genelist_DUP_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
control_genes_DUP_exclusive <- bitr(control_genelist_DUP_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
case_entrez_DUP_exclusive <- case_genes_DUP_exclusive$ENTREZID
control_entrez_DUP_exclusive <- control_genes_DUP_exclusive$ENTREZID
#DELs
shared_genes_DEL <- intersect(case_genelist_DEL, control_genelist_DEL)
case_genelist_DEL_exclusive <- setdiff(case_genelist_DEL, shared_genes_DEL)
control_genelist_DEL_exclusive <- setdiff(control_genelist_DEL, shared_genes_DEL)
case_genes_DEL_exclusive <- bitr(case_genelist_DEL_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
control_genes_DEL_exclusive <- bitr(control_genelist_DEL_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
case_entrez_DEL_exclusive <- case_genes_DEL_exclusive$ENTREZID
control_entrez_DEL_exclusive <- control_genes_DEL_exclusive$ENTREZID


#Reactome overrepresentation analysis
library("ReactomePA")
case_reactome_DEL <- enrichPathway(gene=case_entrez_DEL, readable=TRUE)
df_case_reactome_DEL <- head(case_reactome_DEL, n=10)

control_reactome_DEL <- enrichPathway(gene=control_entrez_DEL, readable=TRUE)
df_control_reactome_DEL <- head(control_reactome_DEL, n=10)

case_reactome_DUP <- enrichPathway(gene=case_entrez_DUP, readable=TRUE)
df_case_reactome_DUP <- head(case_reactome_DUP, n=10)

control_reactome_DUP <- enrichPathway(gene=control_entrez_DUP)
df_control_reactome_DUP <- head(control_reactome_DUP, n=10)

#Reactome ORA for exclusive genes
case_reactome_DEL_exclusive <- enrichPathway(gene=case_entrez_DEL_exclusive, readable=TRUE)
df_case_reactome_DEL_exclusive <- head(case_reactome_DEL_exclusive, n=10)

control_reactome_DEL_exclusive <- enrichPathway(gene=control_entrez_DEL_exclusive, readable=TRUE)
df_control_reactome_DEL_exclusive <- head(control_reactome_DEL_exclusive, n=10)

case_reactome_DUP_exclusive <- enrichPathway(gene=case_entrez_DUP_exclusive, readable=TRUE)
df_case_reactome_DUP_exclusive <- head(case_reactome_DUP_exclusive, n=10)

control_reactome_DUP_exclusive <- enrichPathway(gene=control_entrez_DUP_exclusive)
df_control_reactome_DUP_exclusive <- head(control_reactome_DUP_exclusive, n=10)


#KEGG overrepresentation analysis
case_kegg_DEL <- enrichKEGG(case_entrez_DEL)
df_case_kegg_DEL <- head(case_kegg_DEL, n=10)

control_kegg_DEL <- enrichKEGG(control_entrez_DEL)
df_control_kegg_DEL <- head(control_kegg_DEL, n=10)

case_kegg_DUP <- enrichKEGG(case_entrez_DUP)
df_case_kegg_DUP <- head(case_kegg_DUP, n=10)

control_kegg_DUP <- enrichKEGG(control_entrez_DUP)
df_control_kegg_DUP <- head(control_kegg_DUP, n=10)

#KEGG ORA for exclusive genes
case_kegg_DEL_exclusive <- enrichKEGG(case_entrez_DEL_exclusive)
df_case_kegg_DEL_exclusive <- head(case_kegg_DEL_exclusive, n=10)

control_kegg_DEL_exclusive <- enrichKEGG(control_entrez_DEL_exclusive)
df_control_kegg_DEL_exclusive <- head(control_kegg_DEL_exclusive, n=10)

case_kegg_DUP_exclusive <- enrichKEGG(case_entrez_DUP_exclusive)
df_case_kegg_DUP_exclusive <- head(case_kegg_DUP_exclusive, n=10)

control_kegg_DUP_exclusive <- enrichKEGG(control_entrez_DUP_exclusive)
df_control_kegg_DUP_exclusive <- head(control_kegg_DUP_exclusive, n=10)


#GO (Cellular Component) overrepresentation
case_go_cc_DEL <- enrichGO(gene = case_entrez_DEL, OrgDb= org.Hs.eg.db,
                           ont = "CC",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_cc_DEL <- head(simplify(case_go_cc_DEL), n=10)

control_go_cc_DEL <- enrichGO(gene = control_entrez_DEL, OrgDb= org.Hs.eg.db,
                              ont = "CC",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_cc_DEL <- head(simplify(control_go_cc_DEL), n=10)

case_go_cc_DUP <- enrichGO(gene = case_entrez_DUP, OrgDb= org.Hs.eg.db,
                           ont = "CC",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_cc_DUP <- head(simplify(case_go_cc_DUP), n=10)

control_go_cc_DUP <- enrichGO(gene = control_entrez_DUP, OrgDb= org.Hs.eg.db,
                              ont = "CC",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_cc_DUP <- head(simplify(control_go_cc_DUP), n=10)

#GO CC ORA for exclusive genes
case_go_cc_DEL_exclusive <- enrichGO(gene = case_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "CC",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_cc_DEL_exclusive <- head(simplify(case_go_cc_DEL_exclusive), n=10)

control_go_cc_DEL_exclusive <- enrichGO(gene = control_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "CC",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_cc_DEL_exclusive <- head(simplify(control_go_cc_DEL_exclusive), n=10)

case_go_cc_DUP_exclusive <- enrichGO(gene = case_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "CC",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_cc_DUP_exclusive <- head(simplify(case_go_cc_DUP_exclusive), n=10)

control_go_cc_DUP_exclusive <- enrichGO(gene = control_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "CC",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_cc_DUP_exclusive <- head(simplify(control_go_cc_DUP_exclusive), n=10)


#GO (Molecular Function) overrepresentation analysis
case_go_mf_DEL <- enrichGO(gene = case_entrez_DEL, OrgDb= org.Hs.eg.db,
                           ont = "MF",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_mf_DEL <- head(simplify(case_go_mf_DEL), n=10)

control_go_mf_DEL <- enrichGO(gene = control_entrez_DEL, OrgDb= org.Hs.eg.db,
                              ont = "MF",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_mf_DEL <- head(simplify(control_go_mf_DEL), n=10)

case_go_mf_DUP <- enrichGO(gene = case_entrez_DUP, OrgDb= org.Hs.eg.db,
                           ont = "MF",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_mf_DUP <- head(simplify(case_go_mf_DUP), n=10)

control_go_mf_DUP <- enrichGO(gene = control_entrez_DUP, OrgDb= org.Hs.eg.db,
                              ont = "MF",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_mf_DUP <- head(simplify(control_go_mf_DUP), n=10)

#GO MF ORA for exclusive genes
case_go_mf_DEL_exclusive <- enrichGO(gene = case_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "MF",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_mf_DEL_exclusive <- head(simplify(case_go_mf_DEL_exclusive), n=10)

control_go_mf_DEL_exclusive <- enrichGO(gene = control_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "MF",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_mf_DEL_exclusive <- head(simplify(control_go_mf_DEL_exclusive), n=10)

case_go_mf_DUP_exclusive <- enrichGO(gene = case_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "MF",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_mf_DUP_exclusive <- head(simplify(case_go_mf_DUP_exclusive), n=10)

control_go_mf_DUP_exclusive <- enrichGO(gene = control_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "MF",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_mf_DUP_exclusive <- head(simplify(control_go_mf_DUP_exclusive), n=10)
df_control_go_mf_DUP_exclusive


#GO (Biological Process) overrepresentation analysis
case_go_bp_DEL <- enrichGO(gene = case_entrez_DEL, OrgDb= org.Hs.eg.db,
                           ont = "BP",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_bp_DEL <- head(simplify(case_go_bp_DEL), n=10)

control_go_bp_DEL <- enrichGO(gene = control_entrez_DEL, OrgDb= org.Hs.eg.db,
                              ont = "BP",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_bp_DEL <- head(simplify(control_go_bp_DEL), n=10)

case_go_bp_DUP <- enrichGO(gene = case_entrez_DUP, OrgDb= org.Hs.eg.db,
                           ont = "BP",  #define category: BP, MF, CC, ALL
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           readable = TRUE)
df_case_go_bp_DUP <- head(simplify(case_go_bp_DUP), n=10)

control_go_bp_DUP <- enrichGO(gene = control_entrez_DUP, OrgDb= org.Hs.eg.db,
                              ont = "BP",  #define category: BP, MF, CC, ALL
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.05,
                              readable = TRUE)
df_control_go_bp_DUP <- head(simplify(control_go_bp_DUP), n=10)

#GO BP ORA for exclusive genes
case_go_bp_DEL_exclusive <- enrichGO(gene = case_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "BP",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_bp_DEL_exclusive <- head(simplify(case_go_bp_DEL_exclusive), n=10)

control_go_bp_DEL_exclusive <- enrichGO(gene = control_entrez_DEL_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "BP",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_bp_DEL_exclusive <- head(simplify(control_go_bp_DEL_exclusive), n=10)

case_go_bp_DUP_exclusive <- enrichGO(gene = case_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                     ont = "BP",  #define category: BP, MF, CC, ALL
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.05,
                                     readable = TRUE)
df_case_go_bp_DUP_exclusive <- head(simplify(case_go_bp_DUP_exclusive), n=10)

control_go_bp_DUP_exclusive <- enrichGO(gene = control_entrez_DUP_exclusive, OrgDb= org.Hs.eg.db,
                                        ont = "BP",  #define category: BP, MF, CC, ALL
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.05,
                                        readable = TRUE)
df_control_go_bp_DUP_exclusive <- head(simplify(control_go_bp_DUP_exclusive), n=10)


#Disease Ontology (DO) overrepresentation analysis
case_do_DEL <- enrichDO(case_entrez_DEL)
df_case_do_DEL <- head(case_do_DEL, n=10)

control_do_DEL <- enrichDO(control_entrez_DEL)
df_control_do_DEL <- head(control_do_DEL, n=10)

case_do_DUP <- enrichDO(case_entrez_DUP)
df_case_do_DUP <- head(case_do_DUP, n=10)

control_do_DUP <- enrichDO(control_entrez_DUP)
df_control_do_DUP <- head(control_do_DUP, n=10)

#DO ORA for exclusive genes
case_do_DEL_exclusive <- enrichDO(case_entrez_DEL_exclusive)
df_case_do_DEL_exclusive <- head(case_do_DEL_exclusive, n=10)

control_do_DEL_exclusive <- enrichDO(control_entrez_DEL_exclusive)
df_control_do_DEL_exclusive <- head(control_do_DEL_exclusive, n=10)

case_do_DUP_exclusive <- enrichDO(case_entrez_DUP_exclusive)
df_case_do_DUP_exclusive <- head(case_do_DUP_exclusive, n=10)

control_do_DUP_exclusive <- enrichDO(control_entrez_DUP_exclusive)
df_control_do_DUP_exclusive <- head(control_do_DUP_exclusive, n=10)


#Wikipathways overrepresentation analysis
case_wp_DEL <- enrichWP(case_entrez_DEL, organism = "Homo sapiens")
df_case_wp_DEL <- head(case_wp_DEL, n=10)

control_wp_DEL <- enrichWP(control_entrez_DEL, organism = "Homo sapiens")
df_control_wp_DEL <- head(control_wp_DEL, n=10)

case_wp_DUP <- enrichWP(case_entrez_DUP, organism = "Homo sapiens")
df_case_wp_DUP <- head(case_wp_DUP, n=10)

control_wp_DUP <- enrichWP(control_entrez_DUP, organism = "Homo sapiens")
df_control_wp_DUP <- head(control_wp_DUP, n=10)

#Wikipathways ORA for exclusive genes
case_wp_DEL_exclusive <- enrichWP(case_entrez_DEL_exclusive, organism = "Homo sapiens")
df_case_wp_DEL_exclusive <- head(case_wp_DEL_exclusive, n=10)

control_wp_DEL_exclusive <- enrichWP(control_entrez_DEL_exclusive, organism = "Homo sapiens")
df_control_wp_DEL_exclusive <- head(control_wp_DEL_exclusive, n=10)

case_wp_DUP_exclusive <- enrichWP(case_entrez_DUP_exclusive, organism = "Homo sapiens")
df_case_wp_DUP_exclusive <- head(case_wp_DUP_exclusive, n=10)

control_wp_DUP_exclusive <- enrichWP(control_entrez_DUP_exclusive, organism = "Homo sapiens")
df_control_wp_DUP_exclusive <- head(control_wp_DUP_exclusive, n=10)