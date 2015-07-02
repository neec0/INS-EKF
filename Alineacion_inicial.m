

%-------------------Limpieza inicial de pantalla, variables y grafica-----

% format LONG;%Defino el formato de visualizacion
% load pr1_dha.txt;%Cargo el fichero
% Xacc=pr1_dha(:, 10)*0.002522;%Guardo en una matriz la columna de aceleracion (x y z) y la multiplico por el factor de escalado 2.522 mg
% Yacc=pr1_dha(:, 12)*0.002522;
% Zacc=pr1_dha(:, 14)*0.002522;
% 
% Xacc_m=mean(Xacc(1:1391, 1));%Realizo la media de los valores en estacionario de la primera parada
% Yacc_m=mean(Yacc(1:1391, 1));%para conocer la media de las componentes de aceleracion del vector gravedad
% Zacc_m=mean(Zacc(1:1391, 1));
% 
% g_module_g=sqrt(power(Xacc_m, 2)+power(Yacc_m, 2)+power(Zacc_m, 2));%Modulo del vector gravedad en g
% g_module_ms2=g_module_g*9.80665;%Modulo del vector gravedad en m/s^2
% 
% pitch_rad=asin(Xacc_m);%Calculo aproximado de pitch en radianes
% roll_rad=-asin(Yacc_m/cos(pitch_rad));%Calculo aproximado de roll en radianes
% pitch_deg=rad2deg(pitch_rad);%Resultado del pitch en grados
% roll_deg=rad2deg(roll_rad);%Resultado del roll en grados
% 
% fprintf('Roll = %.6f [rad]\nPitch = %.6f [rad]\n\nRoll = %.6f [deg]\nPitch = %.6f [deg]\n', roll_rad, pitch_rad, roll_deg, pitch_deg)


function [Ai_rad, Ai_deg, Ai_rad_2] = Alineacion_inicial (F_b_gs)

global T_parada1_s dt_r gravedad
N_muestra_parada1 = round(T_parada1_s * (1 / dt_r));%Numero de la muestra donde termina la parada inicial

Ax_m = mean(F_b_gs(1:N_muestra_parada1, 1));%Realizo la media de los valores en estacionario de la aceleracion
Ay_m = mean(F_b_gs(1:N_muestra_parada1, 2));%con el proposito de detectar el vector gravedad
Az_m = mean(F_b_gs(1:N_muestra_parada1, 3));

P_init =  asin(Ax_m);%Calculo aproximado de pitch en radianes
R_init = -asin(Ay_m/cos(P_init));%Calculo aproximado de roll en radianes
Y_init =  0;
Ai_rad = [R_init, P_init, Y_init];

% R_init = asin(Ay_m / Az_m);%Calculo aproximado de roll en radianes
% P_init = asin(Ax_m / gravedad);%Calculo aproximado de pitch en radianes
% Y_init =  0;
% Ai_rad = [R_init, P_init, Y_init];

R_init_deg = rad2deg(R_init);%Resultado del roll en grados
P_init_deg = rad2deg(P_init);%Resultado del pitch en grados
Y_init_deg = Y_init;
Ai_deg = [R_init_deg, P_init_deg, Y_init_deg];

Cbl = [ cos(P_init) (sin(P_init) * sin(R_init)) (sin(P_init) * cos(R_init));
                  0                 cos(R_init)                -sin(R_init);
       -sin(P_init) (cos(P_init) * sin(R_init)) (cos(P_init) * cos(R_init)) ];

% Tilt angles

acc1 = [Ax_m, Ay_m, Az_m];
acc2 = Cbl * acc1';

tilt_y = acc2(1);
tilt_x = acc2(2);
tilt_z = 0;

% DCM orientation matrix, with yaw=0

CLL = [      1  tilt_z -tilt_y;
       -tilt_z       1  tilt_x;
        tilt_y -tilt_x       1 ];

CBL = Cbl * CLL';

roll  = atan2(CBL(3,2), CBL(3,3));
pitch = -asin(CBL(3,1));

Ai_rad_2 = [roll, pitch, 0];
    
    
    
    
    
    
    