library(tidyverse)
library(readxl)
base = read_excel("base_modelo_AVC.xlsx")


################################## Tratando a base a base ##################################


# Criando as variaveis Dummies para Idade e classificação de diabetes, considerando a categoria Idoso como referencia
  # Classe media = Income == 3,4,5 ; Classe alta = Income == 6,7,8 . "Olhar na base o que isso representa em intervalos de valores"
  # Pré diabetes : 'diabetes' == 1 ; diabetes: 'diabetes' == 2
  # Jovem até 39 anos ; Adulto até 64 anos ; idoso >= 65 anos. Conforme a tabela no link https://www.icpsr.umich.edu/web/NAHDAP/studies/34085/datasets/0001/variables/AGEG5YR?archive=NAHDAP

base = base %>% 
  mutate(`Classe_media` = ifelse(base$Income %in% c(3,4,5),1,0),
         `Classe_alta` = ifelse(base$Income %in% c(6,7,8),1,0),
         `diabetes` = ifelse(base$Diabetes_012 == 2,1,0),
         `pre_diabetes` = ifelse(base$Diabetes_012 == 1,1,0),
         `Jovem` = ifelse(base$Age %in% c(1,2,3,4),1,0),
         `Adulto` = ifelse(base$Age %in% c(5,6,7,8,9),1,0)) %>%
  select(-c(1,12,13))

############################  Verificando as dummies
#length(which(base$Age %in% c(1,2,3,4))) = 4025
#length(which(base$Age %in% c(5,6,7,8,9))) = 16020
#length(which(base$Age %in% c(10,11,12,13))) = 14261
nrow(base[which(base$Jovem == 1),])
nrow(base[which(base$Adulto == 1),])
nrow(base[intersect(which(base$Jovem == 0),which(base$Adulto == 0)),])

# Dessa maneira a categoria de referencia das variaveis explicativas idade, renda e presença de diabetes será "idoso",sem diabetes e de class baixa


# Análise exploratória da base

table(base$Stroke)

prop.table(base$Stroke)
par(mfrow = c(1,2))
boxplot(base$BMI,main = "IMC")
summary(base$BMI)
hist(base$BMI,main = "IMC")

for (i in c(1:3,5:ncol(base))) {
  print(colnames(base)[i])
  print(table(base[,i]))
}


####################################################    Definindo os primeiros Ajustes   ####################################################

# Ajuste com todas variáveis explicativas

ajuste.1 = glm(Stroke ~ .,data = base, family = binomial(link = "logit"))
summary(ajuste.1)

# Ajuste através do do método de seleção stepwise

ajuste.2 = step(ajuste.1, direction = "both", k =2)
summary(ajuste.2)

qchisq(0.95,34293)


# Vabase# Variavel pré-diabetes não foi selecionada para o modelo, logo ela se juntará a categoria de referencia

## Todas variáveis Significativas ao nível de 5%

anova(ajuste.2)

# Verificando Multicolinearidade

vif(ajuste.2)

# Baixa multicolinearidade entre as covariaveis.

# Verificando ajuste do modelo

D = qchisq(0.95,34293);D

# Deviance do ajuste é menor que o quantil 0.95 da distribuição da estatística de teste.
# Então o modelo passa no teste de bondade de ajuste.

# Definindo os valores previstos

pi = ajuste.2$fitted.values
avc_ajustado = ifelse(pi > 0.3,1,0)
table(avc_ajustado)
prop.table(table(avc_ajustado))


######################################################################################################################################################################
######################################################################################################################################################################
######################################################################################################################################################################

#Chances : Razão de Probabilidades

exp(ajuste.2$coefficients)

