# Definindo os vetores
snps_causais <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")
snps_selecionados <- c("SNP2", "SNP5", "SNP6", "SNP10", "SNP14", "SNP25", "SNP56")

# Verificando SNPs que estão em ambos os vetores
snps_presentes <- intersect(snps_causais, snps_selecionados)

# Verificando SNPs selecionados que não estão nos causais
snps_nao_presentes <- setdiff(snps_selecionados, snps_causais)

# Imprimindo os resultados
# Contando os SNPs presentes e não presentes
contagem_presentes <- length(snps_presentes)
contagem_nao_presentes <- length(snps_nao_presentes)

# Imprimindo os resultados
cat("SNPs presentes nos dois vetores:", snps_presentes, "\n")
cat("Quantidade de SNPs presentes nos dois vetores:", contagem_presentes, "\n")
cat("SNPs selecionados que não estão nos causais:", snps_nao_presentes, "\n")
cat("Quantidade de SNPs selecionados que não estão nos causais:", contagem_nao_presentes, "\n")