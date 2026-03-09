#write this string to file
#define file name
sink("my_output_SMS_Corr_1_CAR.txt")

#Imprimindo as informacoes finais para o texto

cat('SMS Linear =',snps_selec_ref[[1]], "\n")

cat('SMS Radial Gamma 0.001 =',snps_selec_ref[[2]], "\n")

cat('SMS Radial Gamma 0.01 =',snps_selec_ref[[3]], "\n")

cat('SMS Radial Gamma 0.1 =',snps_selec_ref[[4]], "\n")

cat('SMS Radial Gamma 1 =',snps_selec_ref[[5]], "\n")

cat('Uniao =',uniao_final, "\n")

cat('Intersecao =',intersecao_final, "\n")

cat('Valor-p bruto =',valor_p_bruto_selecao, "\n")

cat('Valor-p corrigido =',valor_p_corrigido_selecao, "\n")

################# Medidas de Otimalidade dos SNPs causais###################

var_sel_causais<-c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8")
