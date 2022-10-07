

#Programa para resolucao numerica da Equacao de Poisson e 
#calculo do campo eletrico no espaco.


tempo1 = time(); #para calcular o tempo que o programa demora para rodar
#*********PARTE QUE PODE SER ALTERADA****************************************************

#Parametros:

#Distancias em m, carga em C, potencial em V.
x_carga = 0.5;    #posicao x da carga entre 0 e L
y_carga = 0.4;    #posicao y da carga entre 0 e L 
Q_eps0 = 4;       #carga sobre epsilon zero da carga pontual 
L = 1;            #limites do espaço
dx = dy = 0.05;   #definicao espacial
e = 10^-3;        #criterio de convergencia
nmax = 1000;      #limite de iteracoes
V1 = -1;          #potencial da placa 1, a esquerda
V2 = +1;          #potencial da placa 2, a direita 


#*************FIM************************************************************************






#*************AJUSTES INICIAIS PARA AMBOS OS CASOS **************************************

x = 0:dx:L;                            #discretizacao do espaco
y = 0:dy:L;                            #discretizacao do espaco
pho_eps0 = zeros(length(y),length(x)); #grid da densidade de carga
x_c = int32(x_carga/dx + 1) ;          #posicao x da carga no grid
y_c = int32(y_carga/dy + 1) ;          #posicao y da carga no grid
pho_eps0(y_c,x_c) = Q_eps0/dx^2;       #colocacao da densidade de carga sobre epsilon zero



#*************AJUSTE DAS CONDICOES DO PROBLEMA - SEM PLACA ******************************

#o potencial em todo os espaco:
V = zeros(length(y),length(x)); #LINHAS SAO O EIXO X E COLUNAS O EIXO Y

#Condicoes de contorno:
#primeira e ultima coluna:
V(:,1) = 0;
V(:,length(x)) = 0;
#primeira e ultima linha:
V(1,:) = 0;
V(length(y),:) = 0;

#****************************************************************************************
#*************CALCULO DO POTENCIAL V - SEM PLACA ****************************************

#calculo do potencial em todo o espaco pelo metodo de relaxacao de Jacobi:

n = 1;      #conta o numero de iteracoes
soma = e+1; #precisao maior do que o criterio de convergencia

while soma >= e || n <= nmax #enquanto a precisao nao tiver sido alcancada e 
                             #nao tiverem sido feitas nmax iteracoes
  
  Vold = V; #guarda o potencial anterior para comparar com o que sera calculado
  
  #Calculo do potencial:
  for i=2:length(y)-1   #linha = eixo y
    for j=2:length(x)-1   #coluna = eixo x
      V(i,j) = 0.25*( Vold(i+1,j)+Vold(i-1,j)+Vold(i,j+1)+Vold(i,j-1) + pho_eps0(i,j)*dx^2 ); #metodo de relaxacao
    endfor
  endfor
  
  #funcao de teste para precisao:
  soma = 0;
  for i=2:length(y)-1   #linha = eixo y
    for j=2:length(x)-1   #coluna = eixo x
      soma = soma + abs(V(i,j)-Vold(i,j));
    endfor
  endfor
  
  n = n + 1; #conta mais uma iteracao
  
endwhile

#*************GRAFICOS - SEM PLACA ******************************************************

#grafico 2D com curvas de nivel de V e campo eletrico
contour(y,x,V,50); hold on;                        # plota as curvas de nivel
plot(x_carga,y_carga,'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','r','MarkerSize',9);       # plota uma esfera na posicao da carga
cb = colorbar;                                     # adiciona uma barra de cor
set(cb,'FontSize',18)                              # tamanho dos numeros na barra de cor
rotate3d on;                                       # permite rodar o grafico
title('Solução da Eq. de Laplace','FontSize',22);  # coloca um título
xlabel('x(m)','FontSize',18);                      # coloca uma legenda no eixo x
ylabel('y(m)','FontSize',18);                      # coloca uma legenda no eixo y
set(gca,'FontSize',18)                             # tamanho dos numeros nos eixos
[Ex Ey] = gradient(-1*V,dy,dx);                    # calcula o campo eletrico
quiver(x,y,Ex,Ey);                                 # plota o E como um campo vetorial
hold off;

#grafico 3D de superficie
figure;                                            # cria uma nova figura
colormap(cool (64));                               # cores do grafico
surf(y,x,V);                                       # plota a superficie
h = legend("potencial");                           # adiciona uma legenda
set(h,'FontSize',16);                              # tamanho do legend
title("Solução da Eq. de Laplace","FontSize",22);  # coloca um título
xlabel('x(m)','FontSize',18);                      # coloca uma legenda no eixo x
ylabel('y(m)','FontSize',18);                      # coloca uma legenda no eixo y
zlabel('Potencial(V)','FontSize',18);              # coloca uma legenda no eixo z
set(gca,'FontSize',18)                             # tamanho dos numeros nos eixos


#grafico do campo eletrico pela distancia
x_E = 1:1:length(x)-x_c;                           #x para plot do campo eletrico (a direita da carga e no mesmo y)
x_E = x_E*dx;                                      #por problemas numericos tive que fazer esses dois passos
E_carga = Ex(y_c,x_c+1:length(x));                 #campo eletrico a direita da carga
E_esperado = Q_eps0/(4*pi) * x_E.^-2;              #campo eletrico proporcional ao inverso da distancia
E_comparacao = Q_eps0/(4*pi) * x_E.^-1;            #campo eletrico porporcional ao inverso do quadrado da distancia
#plot:
figure;                                            # cria uma nova figura
plot(log(x_E),log(E_carga),'r',"LineWidth",2,...   # plota o campo eletrico da carga
    log(x_E),log(E_esperado),'k',"LineWidth",2,... # plota o campo eletrico esperado
    log(x_E),log(E_comparacao),'g',"LineWidth",2); # plota o campo eletrico de comparacao
h = legend("E da carga pontual",...
          "E proporcional ao inverso da distancia",...
          "E porporcional ao inverso do quadrado da distancia");# adiciona uma legenda
ylim([-2 7]);
set(h,"location", "northeast",'FontSize',16);                              # tamanho do legend
title("Comportamento de E","FontSize",22);         # coloca um título
xlabel('log(x(m))','FontSize',18);                 # coloca uma legenda no eixo x
ylabel('log(|E|)','FontSize',18);                  # coloca uma legenda no eixo yset(gca,'FontSize',18)                             # tamanho dos numeros nos eixos
set(gca,'FontSize',18);                            # tamanho dos numeros nos eixos


#*************FIM - SEM PLACA ***********************************************************
#****************************************************************************************




#*************AJUSTE DAS CONDICOES DO PROBLEMA - COM PLACA ******************************

x1 = 0.3;       #posicao das placas
x2 = 0.7;       #posicao das placas
y1 = 0.2;       #posicao das placas
y2 = 0.8;       #posicao das placas

#condicoes de contorno para as placas e carga:
x_1 = int32(x1/dx + 1);  #calcula o indice do valor de x1
x_2 = int32(x2/dx + 1);  #calcula o indice do valor de x2
y_1 = int32(y1/dy + 1);  #calcula o indice do valor de y1
y_2 = int32(y2/dy + 1);  #calcula o indice do valor de y2

#o potencial em todo os espaco:
V = zeros(length(y),length(x)); #LINHAS SAO O EIXO X E COLUNAS O EIXO Y

#Condicoes de contorno:
#primeira e ultima coluna:
V(:,1) = 0;
V(:,length(x)) = 0;
#primeira e ultima linha:
V(1,:) = 0;
V(length(y),:) = 0;

#determina o potencial das placas:
for i=y_1:y_2
  V(i,x_1) = V1;
  V(i,x_2) = V2;
endfor

#****************************************************************************************
#*************CALCULO DO POTENCIAL V - COM PLACA ****************************************

#calculo do potencial em todo o espaco pelo metodo de relaxacao de Jacobi:

n = 1;      #conta o numero de iteracoes
soma = e+1; #precisao maior do que o criterio de convergencia

while soma >= e || n <= nmax #enquanto a precisao nao tiver sido alcancada e 
                             #nao tiverem sido feitas nmax iteracoes
  
  Vold = V; #guarda o potencial anterior para comparar com o que sera calculado
  
  #Calculo do potencial:
  for i=2:length(y)-1   #linha = eixo y
    for j=2:length(x)-1   #coluna = eixo x
      #condicao de (i,j) nao estar nas placas ou na carga:
      if ((j~=x_1 && j~=x_2) || (i<y_1 || i>y_2)) 
        V(i,j) = 0.25*( Vold(i+1,j)+Vold(i-1,j)+Vold(i,j+1)+Vold(i,j-1) + pho_eps0(i,j)*dx^2 ); #metodo de relaxacao
      endif
    endfor
  endfor
  
  #funcao de teste para precisao:
  soma = 0;
  for i=2:length(y)-1   #linha = eixo y
    for j=2:length(x)-1   #coluna = eixo x
      soma = soma + abs(V(i,j)-Vold(i,j));
    endfor
  endfor
  
  n = n + 1; #conta mais uma iteracao
  
endwhile

#*************GRAFICOS - COM PLACA ******************************************************

#grafico 2D com curvas de nivel de V e campo eletrico
figure;
contour(y,x,V,50); hold on;                        # plota as curvas de nivel
plot(x_carga,y_carga,'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','r','MarkerSize',9);       # plota uma esfera na posicao da carga
cb = colorbar;                                     # adiciona uma barra de cor
set(cb,'FontSize',18)                              # tamanho dos numeros na barra de cor
rotate3d on;                                       # permite rodar o grafico
title('Solução da Eq. de Laplace','FontSize',22);  # coloca um título
xlabel('x(m)','FontSize',18);                      # coloca uma legenda no eixo x
ylabel('y(m)','FontSize',18);                      # coloca uma legenda no eixo y
set(gca,'FontSize',18)                             # tamanho dos numeros nos eixos
[Ex Ey] = gradient(-1*V,dy,dx);                    # calcula o campo eletrico
quiver(x,y,Ex,Ey);                                 # plota o E como um campo vetorial
hold off;

#grafico 3D de superficie
figure;                                            # cria uma nova figura
colormap(cool (64));                               # cores do grafico
surf(y,x,V);                                       # plota a superficie
h = legend("potencial");                           # adiciona uma legenda
set(h,'FontSize',16);                              # tamanho do legend
title("Solução da Eq. de Laplace","FontSize",22);  # coloca um título
xlabel('x(m)','FontSize',18);                      # coloca uma legenda no eixo x
ylabel('y(m)','FontSize',18);                      # coloca uma legenda no eixo y
zlabel('Potencial(V)','FontSize',18);              # coloca uma legenda no eixo z
set(gca,'FontSize',18)                             # tamanho dos numeros nos eixos

#*************FIM - COM PLACA ***********************************************************
#****************************************************************************************



#*************MARCA O TEMPO TOTAL DO PROGRAMA: ******************************************

tempo2 = time(); #tempo que o programa roda
#tempo que o programa demorou:
disp('')
disp ('Tempo de execucao do programa:')
disp('(hh:mm:ss)')
disp(datestr((tempo2-tempo1)/(24*60*60), 'HH:MM:SS'))


##h = legend('');                                  # adiciona uma legenda
##set(h,'FontSize',20)                             # tamanho do legend
##set(gca,'FontSize',20)                           # tamanho dos numeros nos eixos
##legenda = strcat('',num2str(t));                 # gera uma string
##text(x,ylimites(2)*0.8,legenda,'FontSize',16);   # coloca uma texto para o tempo
