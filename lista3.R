library(readxl)
setwd("C:/Users/Guto/Desktop/Planejamento_Experimentos")
dados = read_excel("dados_exer2.xlsx")
names(dados)
######################################################### Questão 2 ########################################################

layout(matrix(c(1,2,3,4,5,6), ncol=3))
boxplot(dados$Resposta~dados$`Pista(Linhas)`)
boxplot(dados$Resposta~dados$`Horario(Colunas)`)
boxplot(dados$Resposta~dados$Tratamentos)
boxplot(dados$Resposta~dados$Solado)
boxplot(dados$Resposta~dados$Firmeza)
boxplot(dados$Resposta~dados$Cadarço)

####################################################        (a)    ####################################################

dados$`Pista(Linhas)`= as.factor(dados$`Pista(Linhas)`)
dados$`Horario(Colunas)`= as.factor(dados$`Horario(Colunas)`)
dados$Tratamentos= as.factor(dados$Tratamentos)

modelo = aov(`Resposta` ~ `Solado` * `Firmeza` *`Cadarço` + `Horario(Colunas)` + `Pista(Linhas)` ,data = dados)
summary(modelo)

modelo$coefficients


####### Análise de Resíduos 
res = residuals(modelo)

mean(res)
summary(res)
plot(res)

qqnorm(res)
qqline(res)
hist(res)
ks.test(res,pnorm)


####################################################        (b)    ####################################################

trat=factor(dados$Tratamentos)
contr=c(1,rep(-1/11,11))
contrasts(trat)=contr
resposta=dados$Resposta
pista=dados$`Pista(Linhas)`
horario=dados$`Horario(Colunas)`

# Ajustando o modelo para comparação

modelo2=aov(resposta~factor(pista)+factor(horario) + trat)

summary(modelo2,split=list(trat=list("Padrão x Alternativos"=1)))

