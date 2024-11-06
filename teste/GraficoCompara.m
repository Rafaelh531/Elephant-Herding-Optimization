%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nome do Arquivo: GraficoCompara.m
% Descrição: Comparação entre os algoritmos binEHO e crEHO.
%
% Autor: Rafael Henrique da Rosa
% Data: 07/11/2024
% Disciplina: Otimização Combinatória
% Professor: Edgar Manuel Carreño Franco


% OBS: descomentar o problema que se deseja resolver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
close all;
tic;


%% Data

%P1
max_peso=165;
peso=[23 31 29 44 53 38 63 85 89 82];
custo=[92 57 49 68 60 43 67 84 87 72];
otimo = 309;

%P2
% max_peso=750;
% peso=[ 70 73 77 80 82 87 90 94 98 106 110 113 115 118 120];
% custo=[135 139 149 150 156 163 173 184 192 201 210 214 221 229 240];
% otimo = 1458;

%P3
%max_peso =6404180;
%peso=[ 382745 799601 909247 729069 467902 44328 34610 698150 823460 903959 853665 551830 610856 670702 488960 951111 323046 446298 931161 31385 496951 264724 224916 169684];
%custo=[825594 1677009 1676628 1523970 943972 97426 69666 1296457 1679693 1902996 1844992 1049289 1252836 1319836 953277 2067538  675367 853655 1826027 65731 901489 577243 466257 369261];
%otimo=[1 1 0 1 1 1 0 0 0 1 1 0 1 0 0 1 0 0 0 0 0 1 1 1]; 13.549.094
%otimo = 13549094;


%% Parâmetros Comuns
n_cla = 1; % Número de clãs (grupos de soluções)
n_individuos = 15; % Número de indivíduos em cada clã
max_geracoes = 50; % Número máximo de gerações


%% Parâmetros binEHO
DR = 0.2;

%% Parâmetros crEHO
alpha = 0.25; % Parâmetro de ajuste para a atualização das posições
beta = 0.4;  % Parâmetro de ajuste para a atualização das posições


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICIA O ALGORITMO binEHO

%% Inicialização
dim = size(peso, 2); % Dimensão do problema (número de itens)
DCount = ceil(DR*dim);
%disp(DR);
top = 0;

    % Geração inicial aleatória de soluções (números entre 0 e 1)
    X = rand(n_cla * n_individuos, dim);
    % Arredondamento para 0 ou 1 (soluções binárias)
    Y = round(X);
    XXXX = Y;
    %% Inicialização das soluções
    for i = 1:n_cla % Itera sobre cada clã
        for j = 1:n_individuos % Itera sobre cada elefante (solução) do clã
            % Calcula o peso total da solução
            peso_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
            % Calcula o custo total da solução
            valor_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);

            % Se o peso excede o limite, ajusta a solução
            while peso_solucao(i,j) > max_peso
                itens_selecionados = find(Y((i-1)*n_individuos+j, :) == 1); % Itens na solução
                item_aleatorio = ceil(rand * length(itens_selecionados)); % Seleciona item aleatório
                % Remove o item escolhido diretamente na matriz Y
                Y((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) = 0;
                % Recalcula o peso e o valor da solução ajustada
                peso_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
                valor_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
            end
            % Fitness é igual ao valor total da solução
            aptidao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
        end
        % A melhor solução do clã (matriarca) é a de maior fitness
        melhor = find(aptidao(i, :) == max(aptidao(i, :)));
        matriarca(i, :) = Y((i-1)*n_individuos + melhor(1), :); % Atualiza matriarca
        fitness_matriarca(i) = max(aptidao(i, :)); % Fitness da matriarca

    end

    %% Otimização
    melhor_geral(1) = max(fitness_matriarca); % Armazena a melhor solução
    t = 0; % Contador de gerações

    while t < max_geracoes

        % Atualiza as posições dos elefantes em cada clã
        for i = 1:n_cla

            e = matriarca(i, :);
            for j = 1:n_individuos
                if Y((i-1)*n_individuos+j, :) == round(matriarca(i, :)) %se for a matriarca
                    l = randi([1, dim]);
                    Y(j,l) = ~Y(j,l);
                else % se nao for a matriarca
                    Dimensionpool =   randi([1, dim], 1, DCount);
                    for k = 1:length(Dimensionpool)
                        Y(j,Dimensionpool(k)) = matriarca(i,Dimensionpool(k));
                    end
                end

                % Recalcula o peso e o valor das soluções
                peso_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
                valor_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* custo);

                % Se o peso excede o limite, ajusta a solução
                while peso_solucao(i,j) > max_peso
                    itens_selecionados = find(Y((i-1)*n_individuos+j, :) == 1); % Itens na solução
                    item_aleatorio = ceil(rand * length(itens_selecionados)); % Seleciona item aleatório
                    % Remove o item escolhido diretamente na matriz Y
                    Y((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) = 0;
                    % Recalcula o peso e o valor da solução ajustada
                    peso_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
                    valor_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
                end
                % Fitness é igual ao valor total da solução
                aptidao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
            end

            % Substitui as piores soluções do clã
            pior = find(aptidao(i, :) == min(aptidao(i, :)));
            if length(pior) == n_individuos
                pior = ceil(0.75 * n_individuos):n_individuos;
            end

            Y((i-1)*n_individuos+pior, :) = randi([0, 1], length(pior), dim);
            % Ajusta as novas soluções geradas
            for k = 1:length(pior)

                % Recalcula o peso e o valor das soluções
                peso_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* peso);
                valor_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* custo);

                % Se o peso excede o limite, ajusta a solução
                while peso_solucao(i, pior(k)) > max_peso
                    itens_selecionados = find(Y((i-1)*n_individuos+pior(k), :) == 1);
                    item_aleatorio = ceil(rand * length(itens_selecionados));
                    Y((i-1)*n_individuos+pior(k), itens_selecionados(item_aleatorio)) = 0;
                    peso_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* peso);
                    valor_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* custo);
                end

                % Atualiza a aptidão da solução ajustada
                aptidao(i, pior(k)) = valor_solucao(i, pior(k));
            end

            pior = find(aptidao(i, :) == min(aptidao(i, :)));
            pior = pior(1);
            Y((i-1)*n_individuos+pior, :) = e;

            % Atualiza a matriarca se houver uma nova melhor solução no clã
            if max(aptidao(i, :)) > fitness_matriarca(i)
                melhor = find(aptidao(i, :) == max(aptidao(i, :)));
                matriarca(i, :) = Y((i-1)*n_individuos + melhor(1), :); % Atualiza matriarca
                fitness_matriarca(i) = max(aptidao(i, :)); % Fitness da matriarca

            end
        end

        % Armazena a melhor solução da geração
        melhor_geral(t+2) = max(fitness_matriarca);

        % Próxima geração
        t = t + 1;
    end

    % Exibe a melhor solução final
    melhor = find(fitness_matriarca == max(fitness_matriarca));
    solucao_final = round(matriarca(melhor(1), :)); % A solução final encontrada
    custo_final = sum(solucao_final .* custo);
    top= top + custo_final;


% Exibe a melhor solução final
melhor = find(fitness_matriarca == max(fitness_matriarca));
solucao_final = round(matriarca(melhor(1), :)); % A solução final encontrada
custo_final = sum(solucao_final .* custo)
toc; % Tempo de execução

%limpa todas as variáveis que não sao o resultado do primeiro algoritmo e
%os dados do problema
graficobin = melhor_geral;
clearvars -except graficobin XXXX max_peso peso custo otimo n_cla n_individuos alpha beta max_geracoes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICIA O ALGORITMO crEHO
tic;

%% Inicialização
dim = size(peso, 2); % Dimensão do problema (número de itens)

% Geração inicial aleatória de soluções (números entre 0 e 1)
X = rand(n_cla * n_individuos, dim); 
% Arredondamento para 0 ou 1 (soluções binárias)
Y = round(X); 
Y = XXXX;
%% Inicialização das soluções
for i = 1:n_cla % Itera sobre cada clã
    for j = 1:n_individuos % Itera sobre cada elefante (solução) do clã
        % Calcula o peso total da solução
        peso_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* peso); 
        % Calcula o custo total da solução
        valor_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo); 
        
        % Se o peso excede o limite, ajusta a solução
        while peso_solucao(i,j) > max_peso
            itens_selecionados = find(Y((i-1)*n_individuos+j, :) == 1); % Itens na solução
            item_aleatorio = ceil(rand * length(itens_selecionados)); % Seleciona item aleatório
            % Remove o item escolhido
            X((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) = X((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) / 2;
            Y((i-1)*n_individuos+j, :) = round(X((i-1)*n_individuos+j, :)); % Arredonda
            % Recalcula o peso e o valor da solução ajustada
            peso_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
            valor_solucao(i,j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
        end
        
        % Fitness é igual ao valor total da solução
        aptidao(i,j) = valor_solucao(i,j); 
    end
    % A melhor solução do clã (matriarca) é a de maior fitness
    melhor = find(aptidao(i, :) == max(aptidao(i, :))); 
    matriarca(i, :) = X((i-1)*n_individuos + melhor(1), :); % Atualiza matriarca
    fitness_matriarca(i) = max(aptidao(i, :)); % Fitness da matriarca
end

%% Otimização
melhor_geral(1) = max(fitness_matriarca); % Armazena a melhor solução
t = 0; % Contador de gerações

while t < max_geracoes
    % Atualiza as posições dos elefantes em cada clã
    for i = 1:n_cla
        elite = matriarca(i, :);
        % Calcula o centro do clã (média das posições)
        centro_cla = sum(X((i-1)*n_individuos+(1:n_individuos), :)) / n_individuos;
        
        % Atualiza as posições dos elefantes em direção à matriarca
        X((i-1)*n_individuos+(1:n_individuos), :) = X((i-1)*n_individuos+(1:n_individuos), :) + alpha * (repmat(matriarca(i, :), n_individuos, 1) - X((i-1)*n_individuos+(1:n_individuos), :)) .* rand(n_individuos, dim);
        
        % Ajusta a posição dos elefantes
        for j = 1:n_individuos
            if Y((i-1)*n_individuos+j, :) == round(matriarca(i, :))
                % Atualiza baseado no centro do clã
                X((i-1)*n_individuos+j, :) = X((i-1)*n_individuos+j, :) + beta * centro_cla .* rand(1, dim);
            end
            % Normaliza a posição dos elefantes
            if min(X((i-1)*n_individuos+j, :)) < 0
                X((i-1)*n_individuos+j, :) = (X((i-1)*n_individuos+j, :) - min(X((i-1)*n_individuos+j, :))) / (max(X((i-1)*n_individuos+j, :)) - min(X((i-1)*n_individuos+j, :)));
            else
                X((i-1)*n_individuos+j, :) = X((i-1)*n_individuos+j, :) / max(X((i-1)*n_individuos+j, :));
            end
            Y((i-1)*n_individuos+j, :) = round(X((i-1)*n_individuos+j, :)); % Arredonda as soluções
            
            % Recalcula o peso e o valor das soluções
            peso_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
            valor_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
            
            % Se o peso excede o limite, ajusta a solução
            while peso_solucao(i, j) > max_peso
                itens_selecionados = find(Y((i-1)*n_individuos+j, :) == 1);
                item_aleatorio = ceil(rand * length(itens_selecionados));
                X((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) = X((i-1)*n_individuos+j, itens_selecionados(item_aleatorio)) / 2;
                Y((i-1)*n_individuos+j, :) = round(X((i-1)*n_individuos+j, :));
                peso_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* peso);
                valor_solucao(i, j) = sum(Y((i-1)*n_individuos+j, :) .* custo);
            end
            aptidao(i, j) = valor_solucao(i, j); % Atualiza a aptidão da solução
        end
        
        % Substitui as piores soluções do clã
        pior = find(aptidao(i, :) == min(aptidao(i, :)));
        if length(pior) == n_individuos
            pior = ceil(0.75 * n_individuos):n_individuos;
        end
        X((i-1)*n_individuos+pior, :) = repmat(min(X((i-1)*n_individuos+(1:n_individuos), :)),length(pior), 1) + repmat((max(X((i-1)*n_individuos+(1:n_individuos), :)) - min(X((i-1)*n_individuos+(1:n_individuos), :))) + 1,length(pior), 1) .* rand(length(pior), dim);
        
        % Ajusta as novas soluções geradas
        for k = 1:length(pior)
            if min(X((i-1)*n_individuos+pior(k), :)) < 0
                X((i-1)*n_individuos+pior(k), :) = (X((i-1)*n_individuos+pior(k), :) - min(X((i-1)*n_individuos+pior(k), :))) / (max(X((i-1)*n_individuos+pior(k), :)) - min(X((i-1)*n_individuos+pior(k), :)));
            else
                X((i-1)*n_individuos+pior(k), :) = X((i-1)*n_individuos+pior(k), :) / max(X((i-1)*n_individuos+pior(k), :));
            end
            Y((i-1)*n_individuos+pior(k), :) = round(X((i-1)*n_individuos+pior(k), :)); % Arredonda as novas soluções

            % Recalcula o peso e o valor das soluções
            peso_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* peso);
            valor_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* custo);
            
            % Se o peso excede o limite, ajusta a solução
            while peso_solucao(i, pior(k)) > max_peso
                itens_selecionados = find(Y((i-1)*n_individuos+pior(k), :) == 1);
                item_aleatorio = ceil(rand * length(itens_selecionados));
                X((i-1)*n_individuos+pior(k), itens_selecionados(item_aleatorio)) = X((i-1)*n_individuos+pior(k), itens_selecionados(item_aleatorio)) / 2;
                Y((i-1)*n_individuos+pior(k), :) = round(X((i-1)*n_individuos+pior(k), :));
                peso_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* peso);
                valor_solucao(i, pior(k)) = sum(Y((i-1)*n_individuos+pior(k), :) .* custo);
            end
            
            % Atualiza a aptidão da solução ajustada
            aptidao(i, pior(k)) = valor_solucao(i, pior(k));
        end
                    %elitismo

            pior = find(aptidao(i, :) == min(aptidao(i, :)));
            pior = pior(1);
            Y((i-1)*n_individuos+pior, :) = elite;


        % Atualiza a matriarca se houver uma nova melhor solução no clã
        if max(aptidao(i, :)) > fitness_matriarca(i)
            melhor = find(aptidao(i, :) == max(aptidao(i, :)));
            matriarca(i, :) = X((i-1)*n_individuos + melhor(1), :);
            fitness_matriarca(i) = max(aptidao(i, :));
        end
    end

    % Armazena a melhor solução da geração
    melhor_geral(t+2) = max(fitness_matriarca);

    % Próxima geração
    t = t + 1;
end

% Exibe a melhor solução final
melhor = find(fitness_matriarca == max(fitness_matriarca));
solucao_final = round(matriarca(melhor(1), :)); % A solução final encontrada
custo_final = sum(solucao_final .* custo)
toc; % Tempo de execução



figure; % Cria uma nova figura
plot(melhor_geral, 'LineWidth', 1.8, 'DisplayName', 'crEHO'); % Plota a primeira linha
hold on; % Mantém o gráfico para adicionar outra linha
plot(graficobin, 'LineWidth', 1.8, 'DisplayName', 'binEHO'); % Plota a segunda linha

% Configurações do gráfico
xlabel('Iterações'); % Rótulo do eixo X
ylabel('Valores'); % Rótulo do eixo Y
title('Comparação entre crEHO e binEHO - P2'); % Título do gráfico
legend; % Adiciona uma legenda para identificar as linhas
grid on; % Adiciona uma grade ao gráfico para facilitar a leitura
