# Definindo a trilha de dados
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Calculando as correlações dos SNPs selecionados por cada kernel
gamma = 0.01
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "linear"
corr_linear <- validacao_cruzada(dados[[1]][c(snps_selec_ref[[1]],"fenotipo")],folds,gamma,cost,epsilon,kernel)[5]

gamma = 0.001
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"
corr_radial_0001 <-validacao_cruzada(dados[[1]][c(snps_selec_ref[[3]],"fenotipo")],folds,gamma,cost,epsilon,kernel)[5]

gamma = 0.01
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"
corr_radial_001 <-validacao_cruzada(dados[[1]][c(snps_selec_ref[[4]],"fenotipo")],folds,gamma,cost,epsilon,kernel)[5]

gamma = 0.1
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"
corr_radial_01 <-validacao_cruzada(dados[[1]][c(snps_selec_ref[[5]],"fenotipo")],folds,gamma,cost,epsilon,kernel)[5]

gamma = 1
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"
corr_radial_1 <-validacao_cruzada(dados[[1]][c(snps_selec_ref[[5]],"fenotipo")],folds,gamma,cost,epsilon,kernel)[5]


# Definindo os vetores
snps_causais <- c("SNP1", "SNP2", "SNP3", "SNP4", "SNP5", "SNP6", "SNP7", "SNP8")

# Função para processar e formatar a lista de SNPs selecionados
process_snps <- function(snps_causais, snps_selecionados) {
  # Retirando o termo "SNP" e deixando apenas o número
  numeros_causais <- as.numeric(sub("SNP", "", snps_causais))
  numeros_selecionados <- as.numeric(sub("SNP", "", snps_selecionados))
  
  # Ordenando a lista de SNPs selecionados
  numeros_selecionados_ordenados <- sort(numeros_selecionados)
  
  # Verificando quais SNPs selecionados são causais
  causais_presentes <- numeros_selecionados_ordenados[numeros_selecionados_ordenados %in% numeros_causais]
  
  # Criando a lista de SNPs com formatação em negrito para causais
  snp_list <- sapply(numeros_selecionados_ordenados, function(num) {
    if (num %in% causais_presentes) {
      return(paste0("\\textbf{", num, "}"))
    } else {
      return(as.character(num))
    }
  })
  
  # Agrupando os SNPs em linhas de cinco números cada
  snp_list_str <- paste(snp_list, collapse = ", ")
  snp_list_str <- gsub("(([^,]+, ){4}[^,]+), ", "\\1\\\\\\\\", snp_list_str)
  
  # Contando os SNPs presentes e não presentes
  contagem_presentes <- length(causais_presentes)
  
  return(list(snp_list_str = snp_list_str, contagem_presentes = contagem_presentes))
}

# Processando cada lista de SNPs selecionados
result_1  <- process_snps(snps_causais, snps_causais)
result_2  <- process_snps(snps_causais, snps_selec_ref[[1]])
result_3  <- process_snps(snps_causais, snps_selec_ref[[2]])
result_4  <- process_snps(snps_causais, snps_selec_ref[[3]])
result_5  <- process_snps(snps_causais, snps_selec_ref[[4]])
result_6  <- process_snps(snps_causais, snps_selec_ref[[5]])
result_7  <- process_snps(snps_causais, uniao_final)
result_8  <- process_snps(snps_causais, intersecao_final)
result_9  <- process_snps(snps_causais, valor_p_bruto_selecao)
result_10 <- process_snps(snps_causais, valor_p_corrigido_selecao)

# Formatando as correlações para duas casas decimais
corr_linear <- sprintf("%.2f", corr_linear)
corr_radial_0001 <- sprintf("%.2f", corr_radial_0001)
corr_radial_001 <- sprintf("%.2f", corr_radial_001)
corr_radial_01 <- sprintf("%.2f", corr_radial_01)
corr_radial_1 <- sprintf("%.2f", corr_radial_1)

# Criando o conteúdo do arquivo em formato LaTeX
latex_content <- "\\documentclass{article}\n\\usepackage{amsmath}\n\\usepackage{graphicx}\n\\usepackage{array}\n\\begin{document}\n"
latex_content <- paste0(latex_content, "\\begin{table}[h]\n\\caption{Rank CARSCORE para a Simulação 4.}\n\\label{tab 2}\n")
latex_content <- paste0(latex_content, "\\begin{tabular}{c|c|c|c|c}\n\\hline\n")
latex_content <- paste0(latex_content, "\\textbf{Método} & \\textbf{$\\gamma$} & \\textbf{\\begin{tabular}[c]{@{}c@{}}\\emph{SNPs} selecionados\\end{tabular}} & \\textbf{\\# \\emph{SNPs} (V)} & \\textbf{\\begin{tabular}[c]{@{}c@{}}Correlação\\\\média\\end{tabular}} \\\\\n\\hline\n")

# Função para adicionar uma linha na tabela LaTeX
add_to_table <- function(method, gamma, snps, snps_count, correlation) {
  latex_content <<- paste0(latex_content, method, " & ", gamma, " & \\begin{tabular}[c]{@{}c@{}}", snps, "\\end{tabular} & ", snps_count, " & ", correlation, " \\\\\n\\hline\n")
}

# Adicionando as listas de SNPs selecionados e informações adicionais na tabela
add_to_table("SNPs causais", "-", result_1$snp_list_str, "-", "-")
add_to_table("SMS Linear", "-", result_2$snp_list_str, paste(length(snps_selec_ref[[1]]), "(", result_2$contagem_presentes, ")", sep=""), corr_linear)
add_to_table("SMS Radial", "0,001", result_4$snp_list_str, paste(length(snps_selec_ref[[2]]), "(", result_4$contagem_presentes, ")", sep=""), corr_radial_0001)
add_to_table("SMS Radial", "0,01", result_5$snp_list_str, paste(length(snps_selec_ref[[3]]), "(", result_5$contagem_presentes, ")", sep=""), corr_radial_001)
add_to_table("SMS Radial", "0,1", result_5$snp_list_str, paste(length(snps_selec_ref[[4]]), "(", result_5$contagem_presentes, ")", sep=""), corr_radial_01)
add_to_table("SMS Radial", "1", result_6$snp_list_str, paste(length(snps_selec_ref[[5]]), "(", result_6$contagem_presentes, ")", sep=""), corr_radial_1)
add_to_table("União", "-", result_7$snp_list_str, paste(length(uniao_final), "(", result_7$contagem_presentes, ")", sep=""), "-")
add_to_table("Interseção", "-", result_8$snp_list_str, paste(length(intersecao_final), "(", result_8$contagem_presentes, ")", sep=""), "-")
add_to_table("Valor-p bruto", "-", result_9$snp_list_str, paste(length(valor_p_bruto_selecao), "(", result_9$contagem_presentes, ")", sep=""), "-")
add_to_table("Valor-p corrigido", "-", result_10$snp_list_str, paste(length(valor_p_corrigido_selecao), "(", result_10$contagem_presentes, ")", sep=""), "-")

latex_content <- paste0(latex_content, "\\end{tabular}\n\\end{table}\n\\end{document}")

# Escrevendo o conteúdo em um arquivo .tex
writeLines(latex_content, "snps_output_4.tex")

# Imprimindo mensagem de confirmação
cat("Arquivo 'snps_output_4.tex' foi gerado com sucesso.\n")
