clear all;
clc;

tic;

%% Data
%data001 = get(handles.uitable1,'UserData');
%data=data001{1};
%NamaBarang=data001{2};
%Kapasitas=str2num(get(handles.edit6,'string'));
%Modal=str2num(get(handles.edit7,'string'));


maxpeso=165;
peso=[ 23 31 29 44 53 38 63 85 89 82];
custo=[92 57 49 68 60 43 67 84 87 72];
%otimo=[1 1 1 1 0 1 0 0 0 0];
% 0.064184 seconds
dim=size(peso,2)

%% Parameter
nClan= 1; %Número de clãs (ou grupos de soluções).
nci=10; %Número de indivíduos em cada clã.
alpha=0.5; %Parâmetros de ajuste usados na atualização das posições das soluções.
beta=0.5;  %Parâmetros de ajuste usados na atualização das posições das soluções.
MaxGen=30; %Número máximo de gerações.

%% Inisialisasi
X=rand(nClan*nci,dim); % Gera uma matriz de números aleatórios contínuos entre 0 e 1
Y=round(X); %Arredonda os valores para 0 ou 1


%% Inisialisasi
for i=1:nClan %itera cada clan i
    for j=1:nci %itera cada elefante j do clan i
        Berat(i,j)=sum(Y((i-1)*nci+j,:).*peso); %calcua o peso de cada solução
        Beli(i,j)=sum(Y((i-1)*nci+j,:).*custo); %calcula o custo de cada solução
        while Berat(i,j)>maxpeso
            pss=find(Y((i-1)*nci+j,:)==1);
            r1=ceil(rand*length(pss));
            X((i-1)*nci+j,pss(r1))=X((i-1)*nci+j,pss(r1))/2;
            Y((i-1)*nci+j,:)=round(X((i-1)*nci+j,:));
            Berat(i,j)=sum(Y((i-1)*nci+j,:).*peso);
            Beli(i,j)=sum(Y((i-1)*nci+j,:).*custo);
        end
        Fitness(i,j)=Beli(i,j);
    end
    best=find(Fitness(i,:)==max(Fitness(i,:)));
    Matriarch(i,:)=X((i-1)*nci+best(1),:);
    FitMat(i)=max(Fitness(i,:));
end

%%
bsf(1)=max(FitMat);
inon=0;
t=0;

while t<MaxGen
    %Update posisi
    for i=1:nClan
        Xctr=sum(X((i-1)*nci+(1:nci),:))/nci;
        X((i-1)*nci+(1:nci),:)=X((i-1)*nci+(1:nci),:)+alpha*(repmat(Matriarch(i,:),nci,1)-X((i-1)*nci+(1:nci),:)).*rand(nci,dim);
        for j=1:nci
            if Y((i-1)*nci+j,:)==round(Matriarch(i,:))
                X((i-1)*nci+j,:)=X((i-1)*nci+j,:)+beta*Xctr.*rand(1,dim);
            end
            if min(X((i-1)*nci+j,:))<0
                X((i-1)*nci+j,:)=(X((i-1)*nci+j,:)-min(X((i-1)*nci+j,:)))/(max(X((i-1)*nci+j,:))-min(X((i-1)*nci+j,:)));
            else
                X((i-1)*nci+j,:)=X((i-1)*nci+j,:)/max(X((i-1)*nci+j,:));
            end
            Y((i-1)*nci+j,:)=round(X((i-1)*nci+j,:));
            Berat(i,j)=sum(Y((i-1)*nci+j,:).*peso);
            Beli(i,j)=sum(Y((i-1)*nci+j,:).*custo);
            while Berat(i,j)>maxpeso
                pss=find(Y((i-1)*nci+j,:)==1);
                r1=ceil(rand*length(pss));
                X((i-1)*nci+j,pss(r1))=X((i-1)*nci+j,pss(r1))/2;
                Y((i-1)*nci+j,:)=round(X((i-1)*nci+j,:));
                Berat(i,j)=sum(Y((i-1)*nci+j,:).*peso);
                Beli(i,j)=sum(Y((i-1)*nci+j,:).*custo);
            end
            Fitness(i,j)=Beli(i,j);
        end
        %Replace worst
        worst=find(Fitness(i,:)==min(Fitness(i,:)));
        if length(worst)==nci
            worst=ceil(0.75*nci):nci;
        end
        X((i-1)*nci+worst,:)=repmat(min(X((i-1)*nci+(1:nci),:)),length(worst),1)+...
            repmat((max(X((i-1)*nci+(1:nci),:))-min(X((i-1)*nci+(1:nci),:)))+1,length(worst),1).*rand(length(worst),dim);
        for k=1:length(worst)
            if min(X((i-1)*nci+worst(k),:))<0
                X((i-1)*nci+worst(k),:)=(X((i-1)*nci+worst(k),:)-min(X((i-1)*nci+worst(k),:)))/(max(X((i-1)*nci+worst(k),:))-min(X((i-1)*nci+worst(k),:)));
            else
                X((i-1)*nci+worst(k),:)=X((i-1)*nci+worst(k),:)/max(X((i-1)*nci+worst(k),:));
            end
            Y((i-1)*nci+worst(k),:)=round(X((i-1)*nci+worst(k),:));
            Berat(i,worst(k))=sum(Y((i-1)*nci+worst(k),:).*peso);
            Beli(i,worst(k))=sum(Y((i-1)*nci+worst(k),:).*custo);
            while Berat(i,worst(k))>maxpeso
                pss=find(Y((i-1)*nci+worst(k),:)==1);
                r1=ceil(rand*length(pss));
                X((i-1)*nci+worst(k),pss(r1))=X((i-1)*nci+worst(k),pss(r1))/2;
                Y((i-1)*nci+worst(k),:)=round(X((i-1)*nci+worst(k),:));
                Berat(i,worst(k))=sum(Y((i-1)*nci+worst(k),:).*peso);
                Beli(i,worst(k))=sum(Y((i-1)*nci+worst(k),:).*custo);
            end
            Fitness(i,worst(k))= Beli(i,worst(k));
        end
        %Update Matriarch
        if max(Fitness(i,:))>FitMat(i)
            best=find(Fitness(i,:)==max(Fitness(i,:)));
            Matriarch(i,:)=X((i-1)*nci+best(1),:);
            FitMat(i)=max(Fitness(i,:));
        end
    end
    %Plot
    bsf(t+2)=max(FitMat);
    if bsf(t+2)~=bsf(t+1)
        inon=t+1;
    end
    plot(0:length(bsf)-1,bsf,'r','LineWidth',2);
    line(inon,bsf(inon+1),'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',6);
    if inon<t*0.75
        text(inon,bsf(inon+1),sprintf(['\n\n Iter: ' num2str(inon) '\n Fit: ' sprintf('%12d',bsf(inon+1)) ',-']));
    else
        text(inon,bsf(inon+1),sprintf([' Iter: ' num2str(inon) '\nFit: ' sprintf('%12d',bsf(inon+1)) ',-''\n\n']),'HorizontalAlignment','Right');
    end
    ylim([min(bsf) max(bsf)*1.08]);
    %set(handles.axes1,'XColor',[1 1 1],'YColor',[1 11],'FontSize',8,'FontWeight','bold');
    xlabel('Iterasi'); ylabel('Total Profit (Fitness)');
    pause(0.00001);
    %Next iter
    t=t+1;
end
best=find(FitMat==max(FitMat))
FinalSol=round(Matriarch(best(1),:))
toc;
