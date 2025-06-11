#Exercício

getwd()
#setwd("C:/Users/guto/Documents/CodigosR")

# Descricao do codigo: este codigo foi desenvolvido e utilizado em minha graduacao em estatistica na disciplina de Estatistica Computacional 2.
# Nele, estão desenvolvidos diversas funcoes( algumas ja implementadas em R) utilizadas como ferramentas para fins estatísticos, como:
#   - Resolucao de sistemas triangulares
#   - Decomposicao de Matrizes(Chlesky e LU)
#   - Inversao de Matrizes
#   - Métodos Iterativos para resolucao de sistemas de equacoes
#   - Método de Newton para achar raizes de equacoes nao lineares
#   - Geracao de Variaveis aleatorias
#   - Integracao Numerica
#   - Processos de Reamostragem


#######  Pontos a melhorar no codigo: Comentar melhor cada método, inserir exemplos, mostrar codigos nativos que fazem a mesma coisa.#############


u = c(2,0,0,0,0,1,6,0,0,0,4,7,8,0,0,3,9,12,10,0,10,14,16,15,20)
n = sqrt(length(u))
U = matrix(u,n,n)
L = t(U)
b=c(-1,2,4,5,1)
A = U+L
B = A+L

############################### Resolução de Sistemas triângulares Superiores(Inversa da transposta = transposta da inversa.)

tri_sup = function(U,b,n){
  X = rep(0,n)
  X[n] = (U[n,n]^-1)*b[n]
  for (i in (n-1):1) {
    somatorio = 0
    for (j in 1:n) {
    somatorio = somatorio + U[i,j]*X[j]
    }
    X[i] = (U[i,i]^-1)*(b[i]-somatorio)
  }
  return(X)
}
#13/04: Código testado


############################### Resolução de Sistemas triângulares Inferiores
tri_inf = function(L,b,n){
  X= rep(0,n)
  X[1] = (L[1,1]^-1)*b[1]
  for (i in 2:n){
    somatorio = 0
    for (j in 1:(i-1)) {
      somatorio = somatorio + L[i,j]*X[j]  
    }
    X[i] = (L[i,i]^-1)*(b[i] - somatorio)
  }
  return(X)
}

#13/04: Código testado



############################### Cholesky

fator_cholesky = function(A,n){
  
  fatorcho=matrix(rep(0,n^2),n,n) #Preenchendo a matrix por colunas

  fatorcho[1,1]=(A[1,1])^1/2
  for (j in 2:n) {
    fatorcho[1,j]=A[1,j]/fatorcho[1,1]
  }
  if (n != 2) {
    for (i in 2:(n-1)) {
      somatorio1=0
      for (k in 1:(i-1)) {
        somatorio1= somatorio1+(fatorcho[k,i])^2
      }
      fatorcho[i,i]=((A[i,i]-somatorio1))^(1/2)
      for (j in (i+1):n) {
        somatorio2=0
        for (c in 1:(i-1)) {
          somatorio2= somatorio2+fatorcho[c,i]*fatorcho[c,j]
        }
      fatorcho[i,j]= (A[i,j]-somatorio2)/fatorcho[i,i]
      }  
    }
  } 
  somatorio3=0
  for (i in 1:(n-1)) {
    somatorio3=somatorio3+(fatorcho[i,n])^2
  }
  fatorcho[n,n] = (A[n,n]-somatorio3)^(1/2)
  
  return(fatorcho)
}

#13/04: Código testado



############################### Fatoração LDU(L*U) para a Matriz B

FatorLU = function(B){
  n = sqrt(length(B))
  #Definido a Matriz nula L
  fator_L = matrix(rep(0,n^2),n,n)
  #Definindo a Matriz Identidade
  fator_U = matrix(rep(0,n^2),n,n)
  diag(fator_U) = 1

  for (i in 1:n) {
    fator_L[i,1] = B[i,1]
  }
  for (j in 2:n) {
    fator_U[1,j] = B[1,j]/fator_L[1,1] 
  }
  if (n != 2) {
    for (i in 2:(n-1)) {
      for (j in 2:i) {
        somatorio = 0
        for (k in 1:(j-1)) {
          somatorio = somatorio + fator_L[i,k]*fator_U[k,j]
        }
        fator_L[i,j] = B[i,j]- somatorio
      }
      for (j in (i+1):n) {
        somatorio2=0
        for (k in 1:(j-1)) {
          somatorio2 = somatorio2 + fator_L[i,k]*fator_U[k,j]
        }
        fator_U[i,j] = (B[i,j]- somatorio2)/fator_L[i,i]
      }
    }
  } 
  for (j in 2:n) {
    somatorio3=0
    for (k in 1:(j-1)) {
      somatorio3 = somatorio3 + fator_L[n,k]*fator_U[k,j]
    }
    fator_L[n,j] = B[n,j] - somatorio3
  }
  return(list(fator_L,fator_U))
}

#28/03: Código testado


############################### Inversão de Matrizes Triângulares Superiores (Inverte Matriz)
## Inversa de t(A) = t(inversa de A)

Inv_Matriz = function(U){
  n = sqrt(length(U))
  X = matrix(rep(0,n^2),n,n)
#Utilizando um método de Cholesky para obter matrizes triângulares.
#A = t(U)*U , onde U é o fator de Cholesky.
  
  for (i in n:1) {
    X[i,i] = 1/U[i,i]
    if (i != n) {
      for (j in (i+1):n) {
        somatorio = 0
        for (k in (i+1):j) {
          somatorio = somatorio + X[k,j]*U[i,k]
        }
        X[i,j] = -somatorio/U[i,i]
      }
    }
  }
  return(X)
}

#28/03: Código testado


############################### Método de Jacobi


Jacobi = function(B,x,b,e){
  
  n = length(b)
  Matriz_Coeficietes = B
  J_U = matrix(rep(0,length(Matriz_Coeficietes)),n,n)
  J_L = matrix(rep(0,length(Matriz_Coeficietes)),n,n)
  J_D = matrix(rep(0,length(Matriz_Coeficietes)),n,n)
  J_U[upper.tri(J_U)] = Matriz_Coeficietes[upper.tri(Matriz_Coeficietes)]
  J_L[lower.tri(J_L)] = Matriz_Coeficietes[lower.tri(Matriz_Coeficietes)]
  diag(J_D) = diag(Matriz_Coeficietes)

  y = rep(0,n)
  erro = 1
  iteracoes = 0

  while (erro > e){
    
    y = rep(0,n)
    y = solve(J_D) %*% ((-1*(J_L+J_U)%*%x)+b)
    
    iteracoes = iteracoes + 1
    
    erro = sum((y-x)^2)
    
    x = y
  }
  return(list(y,erro,iteracoes))
}

#Código Testado em 15/04

############################### Método de Seidel

Seidel = function(B,x,b,e){
  
  n = length(b)
  
  S_U = matrix(rep(0,n^2),n,n)
  S_L = matrix(rep(0,n^2),n,n)
  S_D = matrix(rep(0,n^2),n,n)
  S_U[upper.tri(S_U)] = B[upper.tri(B)]
  S_L[lower.tri(S_L)] = B[lower.tri(B)]
  diag(S_D) = diag(B)
  
  y = rep(0,n)
  erro = 1
  iteracoes = 0
  
  while (erro > e) {
    y = -1 * solve(S_L + S_D) %*% S_U %*% x + solve(S_L + S_D) %*% b
    erro = sum((y-x)^2)
    x = y
    iteracoes = iteracoes + 1
  }
  
  return(list(y,erro,iteracoes))
  
}

#Código Testado em 15/04


############################### Método de Newton Simples

### Método usado para encontrar uma aproximação para raiz de funções.
#### "du" é a primeira derivada da função avaliada em X do passo anterior
#### "u" é a própria função avaliada em X do passo anterior.  


New_Raph = function(x,u,e){
  
  du = eval(D(u,"x"))
  erro = 1
  iteracoes = 0
  
  while (erro > e) {
    y = x-(du^-1)*eval(u)
    erro = abs((y-x))
    x = y
    iteracoes = iteracoes + 1
  }
  return(list(y,iteracoes))
}


############################### Método de Newton para um vetor de funções.

New_Raph_geral = function(resp_0,
                          vetor_variaveis,
                          vetor_expressoes,
                          e){
  
  n = length(resp_0)
  erro = 100
  iteracoes = 0
  
  
  while (erro > e) {
    
    #Tem que atribuir essas variaveis na mão caso queira mudar o problema
    
    x=resp_0[1]
    y=resp_0[2]
    #z=resp_0[3]
    
    J = matrix(rep(0,n^2),n,n)
    
    for (i in 1:n) {
      for (j in 1:n) {
        J[i,j] = eval(D(vetor_expressoes[i],vetor_variaveis[j]))
      }
    }
    
    U = matrix(rep(0,n),n,1)
    
    for (i in 1:n) {
      U[i,1] = eval(vetor_expressoes[i])
    }
    
    resp_1 = resp_0 - solve(J) %*% U
    
    erro = sum((resp_1-resp_0)^2)
    
    resp_0 = resp_1
    
    iteracoes = iteracoes + 1
  }
  return(list(resp_1,iteracoes))
}


############################### Geração de Váriaveis Aleatórias


############################### Método da Transformação inversa para VAs contínuas.
## Utilize a inversa em função de 'U'

trans_inversa = function(Inv_Fx){
  U = runif(1,0,1)   ##Gerando um valor entre 0 e 1
  X =eval(Inv_Fx)    ##Tomando a inversa e aplicando em U.
  return(X)
}


############################### Método da Aceitação - Rejeição

## Suponha g(.) uma fdp que sabemos gerar valores;
## Queremos gerar valores para um fdp f(.);
## Pelo método da aceitação-rejeição, temos:

### - Um valor será aceito para f(.) se f(.)/g(.) <= c,
###   onde c é o valor máximo que essa razão pode assumir.


AceRej_continua = function(Inv_Fx,fy,gy,c){
  
  X = -1000
  iteracoes = 0
  
  while (X == -1000) {
  
    Y = trans_inversa(Inv_Fx)     # 1)Gerando valor Y para uma dist.(g(.)) a partir do método da Inversa.
    
    U = runif(1,0,1)              # 2)Gerando um número aleatório U.
  
    r = eval(fy)/(c*eval(gy))     # 3)Validando se Y pode pertencer a dist. f(.).
  
    if (U <= r) {
      X=Y
    }
    iteracoes = iteracoes + 1
  }
  return(X)
}


############################### Método da Transformação Inversa para VAs Discretas
trans_inv_discreta = function(fp,n,contem_zero){
  
  if (contem_zero == F){
    
    U = runif(1,0,1)
    y = 1:n
    Y = eval(fp)
    Fy = rep(0,n)
    
    for (i in 1:n) {
      Fy[i] = sum(Y[1:i])
    }
    
    if (U< Fy[1]){
      X = 1
    }
    else{
      for (i in 1:(n-1)) {
        if (U >= Fy[i] && U < Fy[i+1]) {
          X = i+1         #Como Fy[1] é o valor acumulado se x=0, F[i+1]= P(X<=i).
        }
      }
    } 
  }
  
  
  if (contem_zero == T){
    
    U = runif(1,0,1)
    y = 0:n
    Y = eval(fp)
    Fy = rep(0,n)
    
    for (i in 1:n) {
      Fy[i] = sum(Y[1:i])
    }
    
    if (U< Fy[1]){
      X = 0
    }
    else{
      for (i in 1:(n-1)) {
        if (U >= Fy[i] && U < Fy[i+1]) {
          X = i         #Como Fy[1] é o valor acumulado se x=0, F[i+1]= P(X<=i).
        } 
      }
    }
  }
  return(X)
}


############################### Método da Aceitação-Rejeição para VAs Discretas

AceRej_Discreta = function(qy,nq,contem_zero,py,np,c){
  
  X = -1
  while ( X == -1 ) {
    y = trans_inv_discreta(qy,nq,contem_zero) 
    U = runif(1,0,1)
  
    if (y <= np) {
      if (U <= eval(py)/(eval(qy)*c)){
        X = y
      }
    }
  }
  return(X)
}  


############################### Integração Numérica

############################### Integral Retangular

Integral_retangular = function(a,b,n,fx){
  
  h = (b-a)/n
  
  integral_fx = 0
  
  for (i in 1:n) {
    
    x = a + (2*i-1)*h/2
    integral_fx = integral_fx + eval(fx)
    
  }
  return(h*integral_fx)  
}


############################### Integral Trapezio

Integral_trapezio = function(a, b, n, fx){
  
  h=(b-a)/n
  
  integral = 0
  
  for (i in 1:(n-1)) {
    
    x = (a+i*h)
    
    integral = integral + eval(fx)  
    
  }
  
  x = b
  
  integral = h*integral+eval(fx)/2
  return(integral)
}


############################### Um terço de Simpsom

Um_terco_Simp = function(a , b , n ,fx){
  
  h = (b-a)/n
  
  integral = 0
  
  for (i in 1:(n/2)) {
    x = a+(2*i-1)*h
    
    integral = integral + 4*eval(fx)
    
  }
  
  for (i in 1:((n/2)-1)) {
    
    x = a+2*i*h
    
    integral = integral + 2*eval(fx)
    
  }
  
  x = a
  integral = integral + eval(fx)
  x = b
  integral = integral + eval(fx)
  
  return(h/3*integral)
}

  
############################### Aproximação de La Place
  
  #Método consiste em tentar escrever a função original como o kernel de uma normal.
  
  #1 - Queremos calcular uma integral de g(x)
  #2 - La place calcula integrais no formato h(x)exp{Mg(x)}
  #3 - Precisamos achar c, tal que c é o ponto de máximo de g(x)
  #4 - Aproximamos a integral por:
  #
  #             raiz(2pi/M|g''(c)|) * [ Z(b,c,1/M|g''(c)|) - Z(a,c,1/M|g''(c)|) ]
  
  
  
  
  
  
aprox_Laplace = function(a,b,c,hx,M,gx){
    
    x = c
    gc = eval(gx)
    hc = eval(hx)
    dgx = D(gx,"x")
    ddgx = D(dgx,"x")
    ddgc = eval(ddgx)
    
    I = hc*exp(M*gc)*sqrt((2*pi)/(M*abs(ddgc)))*(pnorm(b,c,(1/M*abs(ddgc)))-pnorm(a,c,(1/M*abs(ddgc))))
    
    
    return(I)
    
  }

#-------------------------------------------------------------------------------------------------------------------------
#Monte Carlo Simples

## Calculando como espernça de uma uniforme(a,b)

MC_Unif = function(a,b,n,gx){
  
  x = runif(n,a,b)
  y = eval(gx)
  
  I = (b-a)*mean(y)
  
  return(I)
  
}
  
  
##  Se a integral já for um fdp bastar gerar valores no intervalo e calcular a média

x = rnorm(100000,0,1)

length(x[which(x<=1.96)])/length(x) - length(x[which(x<= -1.96)])/length(x) 


#Exemplo 3 X~Exp(1):
x = rexp(10000,1)

#a) E(X)
media_exp_1 = mean(x)
media_exp_1

#b)P(1<X<3)
P_X_3_1 = (length(x[which(x <=3)]) - length(x[which(x<=1)]) )/length(x)
P_X_3_1

  
#-------------------------------------------------------------------------------------------------------------------------
#Monte Carlo com Amostragem por importância

#Exemplo 4:

#E(X^2), X~ Beta(2,2)

#a)Unif(0,1)

gx = expression(x^2)
fx = expression(6*x*(1-x))
hx = expression(1)
x = runif(1000,0,1) 

esp_2 = mean(eval(gx)*eval(fx)/eval(hx))

var(x)+(mean(x)^2)

  
#-------------------------------------------------------------------------------------------------------------------------
# Reamostragem
  
  
# Jackknife , Bootstrap e algoritmo EM.

  install.packages("bootstrap")
  library(bootstrap)
  library(dplyr)

  x = c(2,3.7,5,1.8,4)
  
 
# Jackknife para amostra x  
  
  
   media_amostra = mean(x)
  
  #teta_chapeu_i = jackknife(x,mean)$jack.values
  
  teta_chapeu_i = NULL
  
  for (i in 1:length(x)) {
    
    teta_chapeu_i[i] = mean(x[-i]) 
    
  }
  
  
  tetas_til = length(x)*media_amostra - (length(x)-1)*teta_chapeu_i
  
  estimador_jack = mean(tetas_til)
  estimador_jack
  
  var_est_jack = 1/(length(x)*(length(x)-1)) *sum((tetas_til-media_amostra)^2)
  
  erro_pad = sqrt(var_est_jack)
  

# Realizando bootstrap para amostra x.
  
  # De forma automatica 
  
  medias_das_amostras = bootstrap(x,2000, mean)
  
  amostras_de_medias = medias_das_amostras$thetastar
  mean(amostras_de_medias)
  var(amostras_de_medias)
  sqrt(var(amostras_de_medias))
  
  # De forma manual
  
  amostras_de_medias = NULL
  B = 50
  for (i in 1:B) {
    
    amostras_de_medias[i] = mean(sample(x,length(x),T))
    
  }
  
  est_boot = mean(amostras_de_medias)
  est_boot
  
  var_est_boot = var(amostras_de_medias)
  var_est_boot
  
  sum((amostras_de_medias-est_boot)^2)/(length(amostras_de_medias)-1)
  

  
  #Exemplo 3
  
  amostra = c(-1.81,0.63,2.22,2.41,2.95,4.16,4.24,4.53,5.09)
  
  n = length(amostra)
  
  media_amostral = mean(amostra)
  var_amostral = var(amostra)
  
  
  #Para Jacknife
  
  
  jack_amostra = jackknife(amostra,var)
  amostra_de_estimadores = jack_amostra[[3]]
  pseudos_valores = n*var_amostral - (n-1)*amostra_de_estimadores
  estimador_da_media_do_estimador = mean(pseudos_valores)
  #Variância amostral de "estimadores ~ " dividido por N
  estimador_da_var_do_estimador = var(pseudos_valores)/n
  Intervalo = c(estimador_da_media_do_estimador+ qt(0.05,8)*sqrt(estimador_da_var_do_estimador),
                estimador_da_media_do_estimador- qt(0.05,8)*sqrt(estimador_da_var_do_estimador))
  Intervalo[1] = 0
  
  # Para bootstrap com B = 10000
  
  amostra_bootstrap_variancias = bootstrap(x,10000,var)[[1]]
  mean(amostra_bootstrap_variancias)
  mean(amostra)
  var(amostra_bootstrap_variancias)
  var(amostra)
  
  estimador_media = mean(amostra_bootstrap_variancias)
  estimador_vicio = estimador_media - media_amostral
  estimador_var = var(amostra_bootstrap_variancias)
  estimador_dp = sqrt(estimador_var)
  intervalo


#Algoritmo EM
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#Estimador Jackknife para media
  
##Amostra
x = c(2,3.7,5,1.8,4)
x
  
#Media amostra completa

teta_media = mean(x)

#Media das amostras jackknife

teta_medias_i = NULL

for (i in 1:5) {
  
  teta_medias_i[i] = mean(x[-i])
  
}

# Pseudo Valores - Tetas Til

teta_media_til = 5*teta_media - 4*teta_medias_i

# Estimativa para a media dos dados

est_media = mean(teta_media_til)

#Vicio estimado do estimado da media

vicio_media = teta_media - est_media
vicio_media

#Estimativa da variancia do estimador

est_var_media = 1/20 * sum((teta_media_til - teta_media)^2)
est_var_media
  
#Erro padrao para estimado da media

erro_pad_media = sqrt(est_var_media)

#Intervalos de Confiança

IC_inf_media = est_media - qt(0.975,5-1)*erro_pad_media
IC_sup_media = est_media + qt(0.975,5-1)*erro_pad_media





#Estimador Jackknife para variancia

##Amostra
x = c(2,3.7,5,1.8,4)
x

##Variancia da amostra completa
teta = var(x)

##variancias das amostras jacknife( 5 amostras com 4 elementos - o elemento i é retirado)
tetas_i = NULL

for (i in 1:5) {
  
  tetas_i[i] = var(x[-i])
  
}  

## Pseudo valores - Tetha Til

teta_til = 5*teta - 4*tetas_i

## Estimador da Variância da amostra - media dos teta til

var_chapeu = mean(teta_til) 

#vicio

vicio_var = teta - var_chapeu
vicio_var

# Estimativa para variância do estimador

est_var = 1/(5*4) * sum((teta_til-teta)^2)
est_var

# Erro padrão do estimador

erro_pad_var = sqrt(est_var)
erro_pad_var

#Intervalo de confiaça

IC_inf = var_chapeu - qt(0.975,4-1)*erro_pad_var
IC_sup = var_chapeu + qt(0.975,4-1)*erro_pad_var










  
  
  
  
  
  
  
  
