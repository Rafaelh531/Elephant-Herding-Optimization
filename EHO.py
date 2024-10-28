import numpy as np

# Dados do problema
max_peso = 6404180
peso = np.array([382745, 799601, 909247, 729069, 467902, 44328, 34610, 698150, 823460, 903959, 
                 853665, 551830, 610856, 670702, 488960, 951111, 323046, 446298, 931161, 31385, 
                 496951, 264724, 224916, 169684])
custo = np.array([825594, 1677009, 1676628, 1523970, 943972, 97426, 69666, 1296457, 1679693, 1902996,
                  1844992, 1049289, 1252836, 1319836, 953277, 2067538, 675367, 853655, 1826027, 65731,
                  901489, 577243, 466257, 369261])

# Parâmetros
n_cla = 1  # Número de clãs (grupos de soluções)
n_individuos = 1000  # Número de indivíduos em cada clã
alpha = 0.5  # Parâmetro de ajuste para a atualização das posições
beta = 0.5   # Parâmetro de ajuste para a atualização das posições
max_geracoes = 300  # Número máximo de gerações

# Inicialização
dim = len(peso)  # Dimensão do problema (número de itens)

# Geração inicial aleatória de soluções (números entre 0 e 1)
X = np.random.rand(n_cla * n_individuos, dim)
# Arredondamento para 0 ou 1 (soluções binárias)
Y = np.round(X)

# Inicialização das soluções
peso_solucao = np.zeros((n_cla, n_individuos))
valor_solucao = np.zeros((n_cla, n_individuos))
aptidao = np.zeros((n_cla, n_individuos))
matriarca = np.zeros((n_cla, dim))
fitness_matriarca = np.zeros(n_cla)

for i in range(n_cla):  # Itera sobre cada clã
    for j in range(n_individuos):  # Itera sobre cada elefante (solução) do clã
        # Calcula o peso total da solução
        peso_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * peso)
        # Calcula o custo total da solução
        valor_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * custo)
        
        # Se o peso excede o limite, ajusta a solução
        while peso_solucao[i, j] > max_peso:
            itens_selecionados = np.where(Y[i * n_individuos + j, :] == 1)[0]  # Itens na solução
            item_aleatorio = np.random.choice(itens_selecionados)  # Seleciona item aleatório
            # Remove o item escolhido
            X[i * n_individuos + j, item_aleatorio] /= 2
            Y[i * n_individuos + j, :] = np.round(X[i * n_individuos + j, :])  # Arredonda
            # Recalcula o peso e o valor da solução ajustada
            peso_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * peso)
            valor_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * custo)
        
        # Fitness é igual ao valor total da solução
        aptidao[i, j] = valor_solucao[i, j]
    
    # A melhor solução do clã (matriarca) é a de maior fitness
    melhor = np.argmax(aptidao[i, :])  # Índice da solução de maior fitness
    matriarca[i, :] = X[i * n_individuos + melhor, :]  # Atualiza matriarca
    fitness_matriarca[i] = aptidao[i, melhor]  # Fitness da matriarca
    
    # Otimização
melhor_geral = [np.max(fitness_matriarca)]  # Armazena a melhor solução
t = 0  # Contador de gerações

while t < max_geracoes:
    # Atualiza as posições dos elefantes em cada clã
    for i in range(n_cla):
        # Calcula o centro do clã (média das posições)
        centro_cla = np.sum(X[i * n_individuos:(i + 1) * n_individuos, :], axis=0) / n_individuos
        
        # Atualiza as posições dos elefantes em direção à matriarca
        X[i * n_individuos:(i + 1) * n_individuos, :] += alpha * (np.tile(matriarca[i, :], (n_individuos, 1)) - X[i * n_individuos:(i + 1) * n_individuos, :]) * np.random.rand(n_individuos, dim)
        
        # Ajusta a posição dos elefantes
        for j in range(n_individuos):
            if np.array_equal(Y[i * n_individuos + j, :], np.round(matriarca[i, :])):
                X[i * n_individuos + j, :] += beta * centro_cla * np.random.rand(dim)
                
            # Normaliza a posição dos elefantes
            if np.min(X[i * n_individuos + j, :]) < 0:
                X[i * n_individuos + j, :] = (X[i * n_individuos + j, :] - np.min(X[i * n_individuos + j, :])) / (np.max(X[i * n_individuos + j, :]) - np.min(X[i * n_individuos + j, :]))
            else:
                X[i * n_individuos + j, :] = X[i * n_individuos + j, :] / np.max(X[i * n_individuos + j, :])
            Y[i * n_individuos + j, :] = np.round(X[i * n_individuos + j, :])
            
            # Recalcula o peso e o valor das soluções
            peso_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * peso)
            valor_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * custo)
            
            # Se o peso excede o limite, ajusta a solução
            while peso_solucao[i, j] > max_peso:
                itens_selecionados = np.where(Y[i * n_individuos + j, :] == 1)[0]
                item_aleatorio = np.random.choice(itens_selecionados)
                X[i * n_individuos + j, item_aleatorio] /= 2
                Y[i * n_individuos + j, :] = np.round(X[i * n_individuos + j, :])
                peso_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * peso)
                valor_solucao[i, j] = np.sum(Y[i * n_individuos + j, :] * custo)
            
            aptidao[i, j] = valor_solucao[i, j]

        # Substitui as piores soluções do clã
        pior = np.where(aptidao[i, :] == np.min(aptidao[i, :]))[0]
        if len(pior) == n_individuos:
            pior = np.arange(int(0.75 * n_individuos), n_individuos)
        X[i * n_individuos + pior, :] = np.tile(np.min(X[i * n_individuos:(i + 1) * n_individuos, :], axis=0), (len(pior), 1)) + (np.tile(np.max(X[i * n_individuos:(i + 1) * n_individuos, :], axis=0) - np.min(X[i * n_individuos:(i + 1) * n_individuos, :], axis=0) + 1, (len(pior), 1))) * np.random.rand(len(pior), dim)

        # Ajusta as novas soluções geradas
        for k in range(len(pior)):
            if np.min(X[i * n_individuos + pior[k], :]) < 0:
                X[i * n_individuos + pior[k], :] = (X[i * n_individuos + pior[k], :] - np.min(X[i * n_individuos + pior[k], :])) / (np.max(X[i * n_individuos + pior[k], :]) - np.min(X[i * n_individuos + pior[k], :]))
            else:
                X[i * n_individuos + pior[k], :] = X[i * n_individuos + pior[k], :] / np.max(X[i * n_individuos + pior[k], :])
            Y[i * n_individuos + pior[k], :] = np.round(X[i * n_individuos + pior[k], :])

            # Recalcula o peso e o valor das soluções
            peso_solucao[i, pior[k]] = np.sum(Y[i * n_individuos + pior[k], :] * peso)
            valor_solucao[i, pior[k]] = np.sum(Y[i * n_individuos + pior[k], :] * custo)

            # Se o peso excede o limite, ajusta a solução
            while peso_solucao[i, pior[k]] > max_peso:
                itens_selecionados = np.where(Y[i * n_individuos + pior[k], :] == 1)[0]
                item_aleatorio = np.random.choice(itens_selecionados)
                X[i * n_individuos + pior[k], item_aleatorio] /= 2
                Y[i * n_individuos + pior[k], :] = np.round(X[i * n_individuos + pior[k], :])
                peso_solucao[i, pior[k]] = np.sum(Y[i * n_individuos + pior[k], :] * peso)
                valor_solucao[i, pior[k]] = np.sum(Y[i * n_individuos + pior[k], :] * custo)
            
            aptidao[i, pior[k]] = valor_solucao[i, pior[k]]

        # Atualiza a matriarca se houver uma nova melhor solução no clã
        if np.max(aptidao[i, :]) > fitness_matriarca[i]:
            melhor = np.argmax(aptidao[i, :])
            matriarca[i, :] = X[i * n_individuos + melhor, :]
            fitness_matriarca[i] = np.max(aptidao[i, :])

    # Armazena a melhor solução da geração
    melhor_geral.append(np.max(fitness_matriarca))

    # Próxima geração
    t += 1

# Exibe a melhor solução final
melhor = np.argmax(fitness_matriarca)
solucao_final = np.round(matriarca[melhor, :])  # A solução final encontrada

print("Solução final:", solucao_final)

valor_total_solucao_final = np.sum(solucao_final * custo)

print("Valor total da solução final:", valor_total_solucao_final)
