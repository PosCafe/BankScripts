library(lattice)

#-----------------------------------------------------#
# Agora faremos executaremos algumas linhas de codico #
# para entendermos o funcionamento basico da EOF      #
#-----------------------------------------------------#

# Criando dados de temperatura

T1 <- c(10,10.4,9.7,9.7,11.7,11,8.7,9.5,10.1,9.6,10.5,9.2,11.3,10.1,8.5)
T2 <- c(10.7,9.8,10,10.1,11.5,10.8,8.8,9.3,9.4,9.6,10.4,9,11.6,9.8,9.2)
x <- data.frame(T1, T2)

## Criando uma Matriz Covariancia

covar <- cov(x)

#  A matrix covariancia eh uma matriz simetrica, ou seja, 
#  quadrada, A = t(A) e diagonalizavel atraves de uma matriz ortogonal.

covar

t(covar)

## Calculando as EOFs, PCs, autovetores ou loadings
pca <- eigen(covar)

# Autovetores
eof <- pca$vectors
eof <- data.frame(e1 = eof[,1], e2 = eof[,2])
# Autovalores
aval <- pca$values

# Logo abaixo serao verificas algumas propriedades dos autovetores e autovalores

## EOF e uma matriz ortogonal
round(eof[,1]%*%eof[,2])

# A multiplicacao da matriz de covariancia pela transpota da matriz EOF e a propria matriz EOF
# resulta  nos autovalores
round(t(eof)%*%covar%*%eof, 2)
round(aval, 2)

## Variancia explicada

var.exp <- round(aval/sum(aval)*100)
p <- barplot(var.exp, ylim = c(0,100), col = c("RoyalBlue", "salmon"))
mtext(side = 1, line = 1, cex = 1.5,
      at = p,  text = c(expression(lambda[1]), expression(lambda[2])))
## Coeficientes de expansao

media <- apply(x, 2, mean)
anom <- t(t(x)-media)

coef <- t(t(eof)%*%t(anom))

## Analisando as EOFs utilizando sua serie temporal

plot(anom, pch = 16, cex = 2)
abline(v=0, col = "red")
abline(h=0, col = "red")
arrows(0, 0, eof[,1], eof[,2], col = "RoyalBlue", lwd = 3)
arrows(0, 0, -eof[,1], -eof[,2], col = "salmon", lwd = 3)

data.frame(x, coef)

#--------------------------------------------------#
# Agora aplicaremos a PCA em dados reais.          #
# Para este fim, utilizaremos os dados disponiveis #  
# no pacote "openair"                              #
#--------------------------------------------------#

library(openair)

# Escolha as variaveis que seram as componentes das EOFs.
# As colunas ja selecionadas abaixo sao apenas uma sugestao.

mydataNA <- na.omit(mydata)
raw <- mydataNA[,c(4,5,6)]

# O primeiro passo para calular uma EOF eh
# fazer a normalizacao do dado.

norm <- apply(raw, 2, function(x) { (x-mean(x))/sd(x) })

# Agora, devomos calcular a covariancia
covar <- cov(norm)

# A covariancia de uma matriz de dados normalizada, nada mais eh
# do que a propria matriz de correlacao.
# Duvida? Entao faca o seguinte teste:
corr <- cor(raw)

# Para calcular as EOFs utilizaremos a funcao eigen,
# como ja fizemos anteriormente

pca <- eigen(covar)

# EOFs
eof <- pca$vectors

# Autovalores
aval <- pca$values

# Variancia explicada
var.exp <- aval/sum(aval)*100
sum(var.exp)

plot(var.exp)

# Coeficiente de expansao ou serie temporal das EOFs
coef <- norm %*% eof

# Montando os data.frames das EOF e coef para facilitar a manipulacao

eof.df <- data.frame(e1 = eof[,1], e2 = eof[,2], e3 = eof[,3])
coef.df <- data.frame(time = mydataNA[,1], e1 = coef[,1], e2 = coef[,2], e3 = coef[,3])

# Plot da serie temporal das EOFs
xyplot(e1+e2+e3 ~ time, data = coef.df, auto.key = TRUE)

# Calculando os percentis de 20 e 80% dos coeficientes
qt <- apply(coef.df[,-1], 2, quantile, c(0.2, 0.8))

plot(coef.df[,1:2], type = 'l')
abline(h = qt[1,1], col = "blue", lwd = 2)
abline(h = qt[2,1], col = "red", lwd = 2)

# Buscando os valores do coef menores que o percentil de 20 %
# e maiores que o percentil de 80% para o primeiro autovetor apenas

flag20 <- coef.df$e1 < qt[1,1]
flag80 <- coef.df$e1 > qt[2,1]

range(coef.df[flag20,2])
range(coef.df[flag80,2])

# Verificando nas serie de dados original (sem NA)
# os padroes encontrados pelas EOF

mydata20 <- mydataNA[flag20,]
mydata80 <- mydataNA[flag80,]

raw20 <- raw[flag20,]
raw80 <- raw[flag80,]

mean.raw <- apply(raw, 2, mean)

mean.raw20 <- apply(raw20, 2, mean)
mean.raw80 <- apply(raw80, 2, mean)

anom.raw20 <- mean.raw20-mean.raw
anom.raw80 <- mean.raw80-mean.raw

