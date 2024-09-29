#algoritmo EHO aplicado no problema da mochila 0-1

import numpy as np


maxpeso=165
peso=[23, 31, 29, 44, 53, 38, 63, 85, 89, 82]
custo=[92, 57, 49, 68, 60, 43, 67, 84, 87, 72]
otimo=[1, 1, 1, 1, 0, 1, 0, 0, 0, 0]
#0.061153 seconds


#parametros:
nClan= 1 #Número de clãs (ou grupos de soluções).
nci=5 #Número de indivíduos em cada clã.
alpha=0.5 #Parâmetros de ajuste usados na atualização das posições das soluções.
beta=0.5  #Parâmetros de ajuste usados na atualização das posições das soluções.
MaxGen=0.5 #Número máximo de gerações.

Berat = np.zeros((nClan, nci)) #peso total das soluções
Beli = np.zeros((nClan, nci)) #peso total das soluções
Fitness = np.zeros((nClan, nci)) #peso total das soluções
Matriarch = np.zeros((nClan,len(peso)))  # Matriz para armazenar a melhor solução de cada população
FitMat = np.zeros(nClan)  # Vetor para armazenar o melhor fitness de cada população

X = np.random.rand(nClan * nci, len(peso))  # Gera uma matriz de números aleatórios contínuos entre 0 e 1
Y = np.round(X)  # Arredonda os valores para 0 ou 1
#print(Y)
for i in range(nClan): #itera em cada clan
    for j in range(nci):  #itera em cada elemento de cada clan
        Berat[i, j] = np.sum(Y[i * nci + j, :] * peso) #calcula o peso de cada solução
        Beli[i, j] = np.sum(Y[i * nci + j, :] * custo) #calcula o custo de cada solução
        while Berat[i, j] > maxpeso:  # Checa somente o peso total
            pss = np.where(Y[i * nci + j, :] == 1)[0]  # Encontra os índices dos itens selecionados
            r1 = np.random.randint(0, len(pss))  # Seleciona aleatoriamente um índice
            X[i * nci + j, pss[r1]] /= 2  # Reduz a quantidade do item selecionado pela metade
            Y[i * nci + j, :] = np.round(X[i * nci + j, :])  # Atualiza Y com valores arredondados
            Berat[i, j] = np.sum(Y[i * nci + j, :] * peso)  # Recalcula o peso total
            Beli[i, j] = np.sum(Y[i * nci + j, :] * custo)  # Recalcula o custo total
        Fitness[i,j] = Beli[i,j] #fitness igual o peso
    best = np.where(Fitness[i, :] == np.max(Fitness[i, :]))[0]  # Encontra os índices dos melhores fitness
    Matriarch[i, :] = X[(i - 1) * nci + best[0], :]  # Atualiza a Matriarca com a solução correspondente ao melhor fitness
    FitMat[i] = np.max(Fitness[i, :])  # Armazena o melhor fitness encontrado para a população i


#print(Y)
#print(Y)
#print("peso:")
#print(Berat)
#print("custo:")
#print(Beli)
#print("matriarca")
#print(np.round(Matriarch))
#print("Fitmat")
#print(FitMat)



bsf[0] = np.max(FitMat)  # Melhor solução inicial
inon = 0
t = 0
