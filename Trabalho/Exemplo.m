 clear all
 close all
 clc

%% Dados topográficos

M = csvread('Trajeto13.csv',1,2);   % Leitura de arquivos no Excel (lembrar de colocar o código e o arquivo.csv na mesma pasta)


x_r  = M(:,2)*1000;             % Transformar x em kilometros para metros
y_r = M(:,1);                   % Posição em y


%% Tensão da linha: Exemplo: 138 kV -> Distância vertical mínima = 8,3m (Rodovia)

y_zs = y_r + 6.8;               % posição em y para definição de zona segura

plot(x_r,y_r,'b',x_r,y_zs,'r'); % plotar gráfico do relevo e zona segura

hold on                        

%% Posicionamento das torres

xt1 = 0;    
xt2 = 518;
xt3 = 1115;
xt4 = 1593;
xt5 = 2190;
xt6 = 2787;
xt7 = 3305;
xt8 = 3862;
xt9 = 4499;
xt10 = 4977;
xt11 = 5534;
xt12 = 6211;
xt13 = 6808;
xt14 = 7406;
xt15 = 8003;
xt16 = 8401;
xt17 = 8998;
xt18 = 8958;
xt19 = 9197;
altura_torres = 37;

i_t1 = find(x_r == xt1);
i_t2 = find(x_r == xt2);
i_t3 = find(x_r == xt3);
i_t4 = find(x_r == xt4);
i_t5 = find(x_r == xt5);
i_t6 = find(x_r == xt6);
i_t7 = find(x_r == xt7);
i_t8 = find(x_r == xt8);
i_t9 = find(x_r == xt9);
i_t10 = find(x_r == xt10);
i_t11 = find(x_r == xt11);
i_t12 = find(x_r == xt12);
i_t13 = find(x_r == xt13);
i_t14 = find(x_r == xt14);
i_t15 = find(x_r == xt15);
i_t16 = find(x_r == xt16);
i_t17 = find(x_r == xt17);
i_t18 = find(x_r == xt18);
i_t19 = find(x_r == xt19);

yt1 = y_r(i_t1)+ altura_torres;
yt2 = y_r(i_t2)+ altura_torres;
yt3 = y_r(i_t3)+ altura_torres;
yt4 = y_r(i_t4)+ altura_torres;
yt5 = y_r(i_t5)+ altura_torres;
yt6 = y_r(i_t6)+ altura_torres;
yt7 = y_r(i_t7)+ altura_torres;
yt8 = y_r(i_t8)+ altura_torres;
yt9 = y_r(i_t9)+ altura_torres;
yt10 = y_r(i_t10)+ altura_torres;
yt11 = y_r(i_t11)+ altura_torres;
yt12 = y_r(i_t12)+ altura_torres;
yt13 = y_r(i_t13)+ altura_torres;
yt14 = y_r(i_t14)+ altura_torres;
yt15 = y_r(i_t15)+ altura_torres;
yt16 = y_r(i_t16)+ altura_torres;
yt17 = y_r(i_t17)+ altura_torres;
yt18 = y_r(i_t18)+ altura_torres;
yt19 = y_r(i_t19)+ altura_torres;


T1 = [xt1 y_r(i_t1);xt1 yt1];
T2 = [xt2 y_r(i_t2);xt2 yt2];
T3 = [xt3 y_r(i_t3);xt3 yt3];
T4 = [xt4 y_r(i_t4);xt4 yt4];
T5 = [xt5 y_r(i_t5);xt5 yt5];
T6 = [xt6 y_r(i_t6);xt6 yt6];
T7 = [xt7 y_r(i_t7);xt7 yt7];
T8 = [xt8 y_r(i_t8);xt8 yt8];
T9 = [xt9 y_r(i_t9);xt9 yt9];
T10 = [xt10 y_r(i_t10);xt10 yt10];
T11 = [xt11 y_r(i_t11);xt11 yt11];
T12 = [xt12 y_r(i_t12);xt12 yt12];
T13 = [xt13 y_r(i_t13);xt13 yt13];
T14 = [xt14 y_r(i_t14);xt14 yt14];
T15 = [xt15 y_r(i_t15);xt15 yt15];
T16 = [xt16 y_r(i_t16);xt16 yt16];
T17 = [xt17 y_r(i_t17);xt17 yt17];
T18 = [xt18 y_r(i_t18);xt18 yt18];
T19 = [xt19 y_r(i_t19);xt19 yt19];

plot(T1(:,1), T1(:,2),'k');
plot(T2(:,1), T2(:,2),'k');
plot(T3(:,1), T3(:,2),'k');
plot(T4(:,1), T4(:,2),'k');
plot(T5(:,1), T5(:,2),'k');
plot(T6(:,1), T6(:,2),'k');
plot(T7(:,1), T7(:,2),'k');
plot(T8(:,1), T8(:,2),'k');
plot(T9(:,1), T9(:,2),'k');
plot(T10(:,1), T10(:,2),'k');
plot(T11(:,1), T11(:,2),'k');
plot(T12(:,1), T12(:,2),'k');
plot(T13(:,1), T13(:,2),'k');
plot(T14(:,1), T14(:,2),'k');
plot(T15(:,1), T15(:,2),'k');
plot(T16(:,1), T16 (:,2),'k');
plot(T17(:,1), T17(:,2),'k');
plot(T18(:,1), T18(:,2),'k');
%% Equações da catenária
%Cabo escolhido: PARROT rated strength:230.50kN ; Weight: 2894 Kg/Km

Tmax = 230500;       % rated strength
Tlimit = 0.25*Tmax;  % limite máximo de 25%
T0 = 0.24*Tmax;     % valor adotado no projeto (<25%)
ms = 2894/1000*9.81; % peso em N/m

fun = @(x0) yt1 - yt2 - T0/ms*(cosh(ms/T0*(xt1-x0))-cosh(ms/T0*(xt2-x0)));  % catenária escrita em função de x0
fun1 = @(x1) yt2 - yt3 - T0/ms*(cosh(ms/T0*(xt2-x1))-cosh(ms/T0*(xt3-x1)));
fun2 = @(x2) yt3 - yt4 - T0/ms*(cosh(ms/T0*(xt3-x2))-cosh(ms/T0*(xt4-x2)));
fun3 = @(x3) yt4 - yt5 - T0/ms*(cosh(ms/T0*(xt4-x3))-cosh(ms/T0*(xt5-x3)));
fun4 = @(x4) yt5 - yt6 - T0/ms*(cosh(ms/T0*(xt5-x4))-cosh(ms/T0*(xt6-x4)));
fun5 = @(x5) yt6 - yt7 - T0/ms*(cosh(ms/T0*(xt6-x5))-cosh(ms/T0*(xt7-x5)));
fun6 = @(x6) yt7 - yt8 - T0/ms*(cosh(ms/T0*(xt7-x6))-cosh(ms/T0*(xt8-x6)));
fun7 = @(x7) yt8 - yt9 - T0/ms*(cosh(ms/T0*(xt8-x7))-cosh(ms/T0*(xt9-x7)));
fun8 = @(x8) yt9 - yt10 - T0/ms*(cosh(ms/T0*(xt9-x8))-cosh(ms/T0*(xt10-x8)));
fun9 = @(x9) yt10 - yt11 - T0/ms*(cosh(ms/T0*(xt10-x9))-cosh(ms/T0*(xt11-x9)));
fun10 = @(x10) yt11 - yt12 - T0/ms*(cosh(ms/T0*(xt11-x10))-cosh(ms/T0*(xt12-x10)));
fun11 = @(x11) yt12 - yt13 - T0/ms*(cosh(ms/T0*(xt12-x11))-cosh(ms/T0*(xt13-x11)));
fun12 = @(x12) yt13 - yt14 - T0/ms*(cosh(ms/T0*(xt13-x12))-cosh(ms/T0*(xt14-x12)));
fun13 = @(x13) yt14 - yt15 - T0/ms*(cosh(ms/T0*(xt14-x13))-cosh(ms/T0*(xt15-x13)));
fun14 = @(x14) yt15 - yt16 - T0/ms*(cosh(ms/T0*(xt15-x16))-cosh(ms/T0*(xt16-x14)));
fun15 = @(x15) yt16 - yt17 - T0/ms*(cosh(ms/T0*(xt16-x17))-cosh(ms/T0*(xt17-x15)));
fun16 = @(x16) yt17 - yt18 - T0/ms*(cosh(ms/T0*(xt17-x18))-cosh(ms/T0*(xt18-x16)));
fun17 = @(x17) yt18 - yt19 - T0/ms*(cosh(ms/T0*(xt18-x19))-cosh(ms/T0*(xt19-x17)));

x0 = fzero(fun, (xt1+xt2)/2);    % achar a raiz x0 da função fun supondo que o valor inicial será próximo a xt1
x1 = fzero(fun1, (xt2+xt3)/2);
x2 = fzero(fun2, (xt4+xt3)/2);
x3 = fzero(fun3, (xt4+xt5)/2);
x4 = fzero(fun4, (xt6+xt5)/2);
x5 = fzero(fun5, (xt6+xt7)/2);
x6 = fzero(fun6, (xt7+xt8)/2);
x7 = fzero(fun7, (xt8+xt9)/2);
x8 = fzero(fun8, (xt9+xt10)/2);
x9 = fzero(fun9, (xt10+xt11)/2);
x10 = fzero(fun10, (xt11+xt12)/2);
x11 = fzero(fun11, (xt12+xt13)/2);
x12 = fzero(fun12, (xt13+xt14)/2);
x13 = fzero(fun6, (xt14+xt15)/2);
x14 = fzero(fun7, (xt15+xt16)/2);
x15 = fzero(fun8, (xt16+xt17)/2);
x16 = fzero(fun9, (xt17+xt18)/2);
x17 = fzero(fun10, (xt18+xt19)/2);

y0 = yt1 - T0/ms*(cosh(ms/T0*(xt1-x0))-1); % achar o valor de y0
y1 = yt2 - T0/ms*(cosh(ms/T0*(xt2-x1))-1);
y2 = yt3 - T0/ms*(cosh(ms/T0*(xt3-x2))-1);
y3 = yt4 - T0/ms*(cosh(ms/T0*(xt4-x3))-1);
y4 = yt5 - T0/ms*(cosh(ms/T0*(xt5-x4))-1);
y5 = yt6 - T0/ms*(cosh(ms/T0*(xt6-x5))-1);
y6 = yt7 - T0/ms*(cosh(ms/T0*(xt7-x6))-1);
y7 = yt8 - T0/ms*(cosh(ms/T0*(xt8-x7))-1);
y8 = yt9 - T0/ms*(cosh(ms/T0*(xt9-x8))-1);
y9 = yt10 - T0/ms*(cosh(ms/T0*(xt10-x9))-1);
y10 = yt11 - T0/ms*(cosh(ms/T0*(xt11-x10))-1);
y11 = yt12 - T0/ms*(cosh(ms/T0*(xt12-x11))-1);
y12 = yt13 - T0/ms*(cosh(ms/T0*(xt13-x12))-1);
y13 = yt14 - T0/ms*(cosh(ms/T0*(xt14-x13))-1);
y14 = yt15 - T0/ms*(cosh(ms/T0*(xt15-x14))-1);
y15 = yt16 - T0/ms*(cosh(ms/T0*(xt16-x15))-1);
y16 = yt17 - T0/ms*(cosh(ms/T0*(xt17-x16))-1);
y17 = yt18 - T0/ms*(cosh(ms/T0*(xt18-x17))-1);
%% Vetores com as coordenadas da catenária com a posição mais baixa na origem
catenaria_a = xt1-x0:xt2-x0;              % definição dos pontos em x da catenária 
catenaria_b = xt2-x1:xt3-x1;
catenaria_c = xt3-x2:xt4-x2;
catenaria_d = xt4-x3:xt5-x3;
catenaria_e = xt5-x4:xt6-x4;
catenaria_f = xt6-x5:xt7-x5;
catenaria_g = xt7-x6:xt8-x6;
catenaria_o = xt8-x7:xt9-x7;
catenaria_p = xt9-x8:xt10-x8;
catenaria_q = xt10-x9:xt11-x9;
catenaria_r = xt11-x10:xt12-x10;
catenaria_s = xt12-x11:xt13-x11;
catenaria_x = xt13-x12:xt14-x12;
catenaria_z = xt14-x13:xt15-x13;
catenaria_a1 = xt15-x14:xt16-x14;
catenaria_b1 = xt16-x15:xt17-x15;
catenaria_c1 = xt17-x16:xt18-x16;
catenaria_d1 = xt18-x17:xt19-x17;


catenaria_h = T0/ms*(cosh(ms/T0*(catenaria_a))-1);        % definição dos pontos m y da catenária
catenaria_i = T0/ms*(cosh(ms/T0*(catenaria_b))-1);
catenaria_j = T0/ms*(cosh(ms/T0*(catenaria_c))-1);
catenaria_k = T0/ms*(cosh(ms/T0*(catenaria_d))-1);
catenaria_l = T0/ms*(cosh(ms/T0*(catenaria_e))-1);
catenaria_m = T0/ms*(cosh(ms/T0*(catenaria_f))-1);
catenaria_n = T0/ms*(cosh(ms/T0*(catenaria_g))-1);
catenaria_t = T0/ms*(cosh(ms/T0*(catenaria_o))-1);
catenaria_u = T0/ms*(cosh(ms/T0*(catenaria_p))-1);
catenaria_v = T0/ms*(cosh(ms/T0*(catenaria_q))-1);
catenaria_w = T0/ms*(cosh(ms/T0*(catenaria_r))-1);
catenaria_y = T0/ms*(cosh(ms/T0*(catenaria_s))-1);
catenaria_f1 = T0/ms*(cosh(ms/T0*(catenaria_x))-1);
catenaria_g1 = T0/ms*(cosh(ms/T0*(catenaria_z))-1);
catenaria_h1 = T0/ms*(cosh(ms/T0*(catenaria_a1))-1);
catenaria_i1 = T0/ms*(cosh(ms/T0*(catenaria_b1))-1);
catenaria_j1 = T0/ms*(cosh(ms/T0*(catenaria_c1))-1);
catenaria_k1 = T0/ms*(cosh(ms/T0*(catenaria_d1))-1);


% Vetores com as coordenadas da catenária na posição real
catenaria_pos_real_a = catenaria_a + x0;
catenaria_pos_real_b = catenaria_b + x1;
catenaria_pos_real_c = catenaria_c + x2;
catenaria_pos_real_d = catenaria_d + x3;
catenaria_pos_real_e = catenaria_e + x4;
catenaria_pos_real_f = catenaria_f + x5;
catenaria_pos_real_g = catenaria_g + x6;
catenaria_pos_real_o = catenaria_o + x7;
catenaria_pos_real_p = catenaria_p + x8;
catenaria_pos_real_q = catenaria_q + x9;
catenaria_pos_real_r = catenaria_r + x10;
catenaria_pos_real_s = catenaria_s + x11;
catenaria_pos_real_x = catenaria_x + x12;
catenaria_pos_real_z = catenaria_z + x13;
catenaria_pos_real_a1 = catenaria_a1 + x14;
catenaria_pos_real_b1 = catenaria_b1 + x15;
catenaria_pos_real_c1 = catenaria_c1 + x16;
catenaria_pos_real_d1 = catenaria_d1 + x17;


catenaria_pos_real_h = catenaria_h + y0;
catenaria_pos_real_i = catenaria_i + y1;
catenaria_pos_real_j = catenaria_j + y2;
catenaria_pos_real_k = catenaria_k + y3;
catenaria_pos_real_l = catenaria_l + y4;
catenaria_pos_real_m = catenaria_m + y5;
catenaria_pos_real_n = catenaria_n + y6;
catenaria_pos_real_t = catenaria_t + y7;
catenaria_pos_real_u = catenaria_u + y8;
catenaria_pos_real_v = catenaria_v + y9;
catenaria_pos_real_w = catenaria_w + y10;
catenaria_pos_real_y = catenaria_y + y11;
catenaria_pos_real_f1 = catenaria_f1 + y12;
catenaria_pos_real_g1 = catenaria_g1 + y13;
catenaria_pos_real_h1 = catenaria_h1 + y14;
catenaria_pos_real_i1 = catenaria_i1 + y15;
catenaria_pos_real_j1 = catenaria_j1 + y16;
catenaria_pos_real_k1 = catenaria_k1 + y17;


% Plotagem da catenária entre as torres
plot(catenaria_pos_real_a, catenaria_pos_real_h,'m');
plot(catenaria_pos_real_b, catenaria_pos_real_i,'m');
plot(catenaria_pos_real_c, catenaria_pos_real_j,'m');
plot(catenaria_pos_real_d, catenaria_pos_real_k,'m');
plot(catenaria_pos_real_e, catenaria_pos_real_l,'m');
plot(catenaria_pos_real_f, catenaria_pos_real_m,'m');
plot(catenaria_pos_real_g, catenaria_pos_real_n,'m');
plot(catenaria_pos_real_o, catenaria_pos_real_t,'m');
plot(catenaria_pos_real_p, catenaria_pos_real_u,'m');
plot(catenaria_pos_real_q, catenaria_pos_real_v,'m');
plot(catenaria_pos_real_r, catenaria_pos_real_w,'m');
plot(catenaria_pos_real_s, catenaria_pos_real_y,'m');
plot(catenaria_pos_real_x, catenaria_pos_real_f1,'m');
plot(catenaria_pos_real_z, catenaria_pos_real_g1,'m');
plot(catenaria_pos_real_a1, catenaria_pos_real_h1,'m');
plot(catenaria_pos_real_b1, catenaria_pos_real_i1,'m');
plot(catenaria_pos_real_c1, catenaria_pos_real_j1,'m');
plot(catenaria_pos_real_d1, catenaria_pos_real_k1,'m');


%% Cálculo da tração em cada ponto da catenária
T = T0 + catenaria_h*ms;
Ta = T0 + catenaria_i*ms;
Tb = T0 + catenaria_j*ms;
Tc = T0 + catenaria_k*ms;
Td = T0 + catenaria_l*ms;
Te = T0 + catenaria_m*ms;
Tf = T0 + catenaria_n*ms;
Tg = T0 + catenaria_t*ms;
Th = T0 + catenaria_u*ms;
Ti = T0 + catenaria_v*ms;
Tj = T0 + catenaria_w*ms;
Tk = T0 + catenaria_y*ms;
Tl = T0 + catenaria_f1*ms;
Tm = T0 + catenaria_g1*ms;
Tn = T0 + catenaria_h1*ms;
To = T0 + catenaria_i1*ms;
Tp = T0 + catenaria_j1*ms;
Tq = T0 + catenaria_k1*ms;


figure

T_lim1 = Tlimit*ones(length(catenaria_a)); % valor limite (25% da tração)
T_lim2 = Tlimit*ones(length(catenaria_b));
T_lim3 = Tlimit*ones(length(catenaria_c));
T_lim4 = Tlimit*ones(length(catenaria_d));
T_lim5 = Tlimit*ones(length(catenaria_e));
T_lim6 = Tlimit*ones(length(catenaria_f));
T_lim7 = Tlimit*ones(length(catenaria_g));
T_lim8 = Tlimit*ones(length(catenaria_o));
T_lim9 = Tlimit*ones(length(catenaria_p));
T_lim10 = Tlimit*ones(length(catenaria_q));
T_lim11 = Tlimit*ones(length(catenaria_r));
T_lim12 = Tlimit*ones(length(catenaria_s));
T_lim13 = Tlimit*ones(length(catenaria_x));
T_lim14 = Tlimit*ones(length(catenaria_z));
T_lim15 = Tlimit*ones(length(catenaria_a1));
T_lim16 = Tlimit*ones(length(catenaria_b1));
T_lim17 = Tlimit*ones(length(catenaria_c1));
T_lim18 = Tlimit*ones(length(catenaria_d1));


plot(catenaria_a, T, 'b', catenaria_a, T_lim1,'r'); % plotagem da tração em cada ponto e do valor limite
plot(catenaria_b, Ta, 'b', catenaria_b, T_lim2,'r');
plot(catenaria_c, Tb, 'b', catenaria_c, T_lim3,'r');
plot(catenaria_d, Tc, 'b', catenaria_d, T_lim4,'r');
plot(catenaria_e, Td, 'b', catenaria_e, T_lim5,'r');
plot(catenaria_f, Te, 'b', catenaria_f, T_lim6,'r');
plot(catenaria_g, Tf, 'b', catenaria_g, T_lim7,'r');
plot(catenaria_o, Tg, 'b', catenaria_o, T_lim8,'r');
plot(catenaria_p, Th, 'b', catenaria_p, T_lim9,'r');
plot(catenaria_q, Ti, 'b', catenaria_q, T_lim10,'r');
plot(catenaria_r, Tj, 'b', catenaria_r, T_lim11,'r');
plot(catenaria_s, Tk, 'b', catenaria_s, T_lim12,'r');
plot(catenaria_x, Tl, 'b', catenaria_x, T_lim13,'r');
plot(catenaria_z, Tm, 'b', catenaria_z, T_lim14,'r');
plot(catenaria_a1, Tn, 'b', catenaria_a1, T_lim15,'r');
plot(catenaria_b1, To, 'b', catenaria_b1, T_lim16,'r');
plot(catenaria_c1, Tp, 'b', catenaria_c1, T_lim17,'r');
plot(catenaria_d1, Tq, 'b', catenaria_d1, T_lim18,'r');


A = T0/ms*(sinh(ms/T0*(xt2-x0))-sinh(ms/T0*(xt1-x0))); % equação para encontrar o comprimento da catenária.
B = T0/ms*(sinh(ms/T0*(xt3-x1))-sinh(ms/T0*(xt2-x1)));
C = T0/ms*(sinh(ms/T0*(xt4-x2))-sinh(ms/T0*(xt3-x2)));
D = T0/ms*(sinh(ms/T0*(xt5-x3))-sinh(ms/T0*(xt4-x3)));
E = T0/ms*(sinh(ms/T0*(xt6-x4))-sinh(ms/T0*(xt5-x4)));
F = T0/ms*(sinh(ms/T0*(xt7-x5))-sinh(ms/T0*(xt6-x5)));
G = T0/ms*(sinh(ms/T0*(xt8-x6))-sinh(ms/T0*(xt7-x6)));
H = T0/ms*(sinh(ms/T0*(xt9-x7))-sinh(ms/T0*(xt8-x7)));
I = T0/ms*(sinh(ms/T0*(xt10-x8))-sinh(ms/T0*(xt9-x8)));
J = T0/ms*(sinh(ms/T0*(xt11-x9))-sinh(ms/T0*(xt10-x9)));
K = T0/ms*(sinh(ms/T0*(xt12-x10))-sinh(ms/T0*(xt11-x10)));
L = T0/ms*(sinh(ms/T0*(xt13-x11))-sinh(ms/T0*(xt12-x11)));
M = T0/ms*(sinh(ms/T0*(xt14-x12))-sinh(ms/T0*(xt13-x12)));
N = T0/ms*(sinh(ms/T0*(xt15-x13))-sinh(ms/T0*(xt14-x13)));
O = T0/ms*(sinh(ms/T0*(xt16-x14))-sinh(ms/T0*(xt15-x14)));
P = T0/ms*(sinh(ms/T0*(xt17-x15))-sinh(ms/T0*(xt16-x15)));
Q = T0/ms*(sinh(ms/T0*(xt18-x16))-sinh(ms/T0*(xt17-x16)));
R = T0/ms*(sinh(ms/T0*(xt19-x17))-sinh(ms/T0*(xt18-x17)));




