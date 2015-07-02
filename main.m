clc;
close all;
clear

global gravedad
global m0 mI Fin_muestras
global Matriz_captura Baliza_detect dt_r T_parada1_s
global f_Hz Ptx_Watt Grx_Watt Gtx_Watt V_media Seg_del_tap Indice_init_Nmuestra Indice_fin_Nmuestra Num_balizas
global Bias_W_rads Bias_F_ng_ms2 Bias_Bal_EKF sigma_Bal
global C
global delta_X_anterior P P_0 k matriz_FI

T_parada1_s = 22;
gravedad = 9.80665;%[m/s^2]
m0 = zeros(3);%Matriz cuadrada de ceros orden 3
mI = eye(3);%Matriz identidad de orden 3
Matriz_captura = xlsread('DATALOG_columnaLab_164sps_4ms_1','DATALOG_columnaLab_164sps_4ms_1','A8:N17466');%Extraigo columnas de datos del excell de la captura
% HMS=datestr(xlsread('DATALOG_cuadrado_1','DATALOG_cuadrado_1','B3'),'HH:MM:SS');%Lectura hora/min/seg inicial captura
% s= second(HMS) + (minute(HMS)*60) + (hour(HMS)*3600);


%Modelado de la campana----------------------------------------------------
% f_Hz = 910;
% Ptx_Watt = 1;
% Grx_Watt = 10;
% Gtx_Watt = 1;
% V_media = 0.8;
% Seg_del_tap = Matriz_captura(:,16);
% Indice_init_Nmuestra = 4275;
% Indice_fin_Nmuestra = 4874;
% Num_balizas = 26;
%--------------------------------------------------------------------------

% 
% %Funciones-----------------------------------------------------------------
% [ W_b_rads, F_b_gs ] = Conversion_unidades ( Matriz_captura );%Funcion de escalado de los datos de la captura
% [ Ai_rad, Ai_deg, Ai_rad_2 ] = Alineacion_inicial ( F_b_gs );%Funcion de estimacion de la orientacion inicial
% [ W_b_nobias_rads, F_b_nobias_gs, sigma2_W, sigma2_F ] = Bias_det ( W_b_rads, F_b_gs, Ai_rad_2 );%Funcion que determina el bias determinista de acelerometros y giroscopos
% [ T_prima, dt_r ] = Reg_lineal_tiempos (  );
% %[ Taps, Senal_baliza, Baliza_detect ] = Modelado_baliza ( T_prima ); %A Seg_del_tap se le asignara Matriz_captura(:,16)
% [ W_b_nb_intrp_rads, F_b_nb_intrp_gs ] = Interpolacion_W_Acc ( W_b_nobias_rads, F_b_nobias_gs, T_prima );
% Fin_muestras = length(W_b_nb_intrp_rads);
% [ Step_detect, C, Abs_Acc, Desv_est, Abs_W, Step_original ] = Detec_paso ( W_b_nb_intrp_rads, F_b_nb_intrp_gs );
% [ Pos, delta_X ] = ec_navegacion ( W_b_nb_intrp_rads, F_b_nb_intrp_gs, Ai_rad_2, T_prima, Step_detect );%Funcion que contiene las ecuaciones de navegacion
% %Fin funciones-------------------------------------------------------------


%Funciones-----------------------------------------------------------------
[ W_b_rads, F_b_gs ] = Conversion_unidades ( );%Funcion de escalado de los datos de la captura
[ T_prima ] = Reg_lineal_tiempos ( );
[ Ai_rad, Ai_deg, Ai_rad_2 ] = Alineacion_inicial ( F_b_gs );%Funcion de estimacion de la orientacion inicial
[ W_b_nbdet_rads, F_b_nbdet_gs ] = Bias_det (W_b_rads, F_b_gs, Ai_rad);%Funcion que determina el bias determinista de acelerometros y giroscopos
W_b_original = W_b_nbdet_rads;
figure, plot(rad2deg(W_b_nbdet_rads(:,1)), 'b');
hold on;
plot(rad2deg(W_b_nbdet_rads(:,2)), 'k');
plot(rad2deg(W_b_nbdet_rads(:,3)), 'r');
xlabel('muestras');
ylabel('Velocidad angular [º/s]');
title('Velocidad angular BIAS determinista');
hold off;
figure, plot(F_b_nbdet_gs(:,1)*gravedad, 'b');
hold on;
plot(F_b_nbdet_gs(:,2)*gravedad, 'k');
plot(F_b_nbdet_gs(:,3)*gravedad, 'r');
xlabel('muestras');
ylabel('Aceleracion [m/s^2]');
title('Aceleracion BIAS determinista');
hold off;
Fin_muestras = length(W_b_nbdet_rads);
[ Step_detect, Abs_Acc, Desv_est, Abs_W, Step_original ] = Detec_paso ( W_b_nbdet_rads, F_b_nbdet_gs );
[ Pos, delta_X ] = ec_navegacion ( W_b_nbdet_rads, F_b_nbdet_gs, Ai_rad_2, T_prima, Step_detect );%Funcion que contiene las ecuaciones de navegacion
%Fin funciones-------------------------------------------------------------


%Representacion gráfica----------------------------------------------------
figure, plot(Pos(:,7),Pos(:,1),'b')
hold on;
plot(Pos(:,7),Pos(:,2),'k')
plot(Pos(:,7),Pos(:,3),'r')
xlabel('Tiempo [s]');
ylabel('Posicion [m]');
title('Posicion X Y Z en funcion del tiempo');
hold off;
N_muestra_parada = T_parada1_s * round(1 / dt_r);
figure, plot(Pos(1:N_muestra_parada,1),Pos(1:N_muestra_parada,2),'r')
hold on;
plot(Pos(N_muestra_parada:9450,1),Pos(N_muestra_parada:9450,2),'g')
plot(Pos(9450:11000,1),Pos(9450:11000,2),'r')
plot(Pos(11000:16600,1),Pos(11000:16600,2),'k')
plot(Pos(16600:end,1),Pos(16600:end,2),'r')
xlabel('Posicion en X [m]');
ylabel('Posicion en Y [m]');
plot(Pos(1:N_muestra_parada,1),Pos(1:N_muestra_parada,2),'r.')
hold on;
plot(Pos(N_muestra_parada:9450,1),Pos(N_muestra_parada:9450,2),'g.')
plot(Pos(9450:11000,1),Pos(9450:11000,2),'r.')
plot(Pos(11000:16600,1),Pos(11000:16600,2),'k.')
plot(Pos(16600:end,1),Pos(16600:end,2),'r.')
xlabel('Posicion en X [m]');
ylabel('Posicion en Y [m]');

plot(Pos(Step_detect == 1, 1), Pos(Step_detect == 1, 2), 'b.')
hold off;
figure, plot(Pos(:,7), rad2deg(Pos(:,4)), 'b')
hold on;
plot(Pos(:,7), rad2deg(Pos(:,5)), 'k')
plot(Pos(:,7), rad2deg(Pos(:,6)), 'r')
xlabel('Tiempo [s]');
ylabel('Giro [º]');
title('Giro X Y Z en funcion del tiempo');
hold off;
figure, plot(T_prima, Pos(:,8), 'b')
hold on;
plot(T_prima, Pos(:,9), 'k')
plot(T_prima, Pos(:,10), 'r')
title('Velocidad lineal y tiempo');
hold off;
%Fin representación gráfica------------------------------------------------