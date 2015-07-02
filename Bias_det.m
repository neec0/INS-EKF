
function [W_b_nbdet_rads, Acc_b_nbdet_ngra_gs] = Bias_det (W_b_rad, F_b_gs, Ai_rad)

global T_parada1_s dt_r gravedad
N_muestra_parada1 = round(T_parada1_s * (1 / dt_r));%Numero de la muestra donde termina la parada inicial


%Promedio de bias deterministas girscopos----------------------------------
W_b_nbdet_rads = [ (W_b_rad(:,1)- mean(W_b_rad(1:N_muestra_parada1 ,1))), (W_b_rad(:,2)- mean(W_b_rad(1:N_muestra_parada1 ,2))), (W_b_rad(:,3)- mean(W_b_rad(1:N_muestra_parada1 ,3))) ];
%Fin-----------------------------------------------------------------------


%Promedio de bias determinista acelerometros-------------------------------
%Declaración de matriz de rotación-----------------------------------------
r = Ai_rad(1);
p = Ai_rad(2);
y = 0;

C1 = [  cos(y) sin(y) 0;
       -sin(y) cos(y) 0;
             0      0 1 ];

C2 = [ cos(p) 0 -sin(p);
            0 1       0;
       sin(p) 0  cos(p) ];

C3 = [ 1       0      0;
       0  cos(r) sin(r);
       0 -sin(r) cos(r) ];

C1 = C1';
C2 = C2';
C3 = C3';

C = C1 * C2 * C3;

% C = [ cos(p) (sin(p) * sin(r)) (sin(p) * cos(r));
%                   0                 cos(r)                -sin(r);
%        -sin(p) (cos(p) * sin(r)) (cos(p) * cos(r)) ];

% Usamos C para rotar las aceleraciones durante la parte estática y poder
% restar la contribución de la gravedad antes de calcular los bias

Acc_b_ngra = zeros(size(F_b_gs));
C_inv = inv(C);

figure, plot(F_b_gs(:,1)*gravedad, 'b')
hold on;
plot(F_b_gs(:,2)*gravedad, 'k')
plot(F_b_gs(:,3)*gravedad, 'r')
title('Aceleracion CON BIAS con gravedad')
xlabel('Aceleracion [m/s^2]')
ylabel('Muestras')

for i=1:length(F_b_gs)
    Acc_b_ngra(i,:) = (C * F_b_gs(i,:)')' + [ 0 0 1 ];%([3X3] * [1X3]')' + [1X3] = [1X3]
    Acc_b_ngra(i,:) =  (C_inv * Acc_b_ngra(i,:)')';%deshacemos el giro
end

Ax_mean = mean(Acc_b_ngra(1:N_muestra_parada1,1));
Ay_mean = mean(Acc_b_ngra(1:N_muestra_parada1,2));
Az_mean = mean(Acc_b_ngra(1:N_muestra_parada1,3));

Acc_b_nbdet_ngra_gs = [ (F_b_gs(:, 1) - Ax_mean), (F_b_gs(:, 2) - Ay_mean), (F_b_gs(:, 3) - Az_mean) ];

plot(Acc_b_nbdet_ngra_gs(:,1)*gravedad, 'g')
plot(Acc_b_nbdet_ngra_gs(:,2)*gravedad, 'c')
plot(Acc_b_nbdet_ngra_gs(:,3)*gravedad, 'y')
hold off;
%Fin-----------------------------------------------------------------------