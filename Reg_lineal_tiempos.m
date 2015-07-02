
function [T_prima] = Reg_lineal_tiempos ()

global Matriz_captura dt_r
%Regresion lineal por minimos cuadrados para la cllumna de timepos del
%fichero Matriz_captura.xls

% derivada_t=diff(Matriz_captura(:, 2));%Derivo la columna de tiempos para convertir en ceros las muestras repetidas0
% indices_t=find(derivada_t);%Extraigo los indices de los valores diferentes de cero
% valores_indices_t=Matriz_captura(indices_t, 2);%Guardo los valores distintos de cero en un array
% x = 0:1:(size(Matriz_captura, 1)-1);
% x = x';
% dt_r = polyfit(x, valores_indices_t, 1);%Regresion lineal de los valores que existen (pendiente y ordenada)
% T=(0:dt_r(1):(size(indices_t)-1)*dt_r(1));%Con la pendiente realizo la ecuacion de la recta en funcion de los valores de muestras
% T=T';
% dt_r=dt_r(1);%Me quedo solo con el valor de la derivada y desprecio el de la ordenada en el origen


% derivada = diff(Matriz_captura(:,2));%se sacan todas la pendientes y se guardan en derivada
% dt_r = mode(derivada);%Dentro de todas las pendientes de derivada tomamos la que mas se repite (moda). Esta es la pendiente buena.
% T_prima = 0 : dt_r : max(Matriz_captura(:,2));%Defino una recta de timepos T' que es mas grande que la recta de tiempos original
% T_prima = T_prima';%Transponemos T_prima
% T = Matriz_captura(:,2);

T = Matriz_captura(:,2);
dt_r = (T(end) - T(1)) / length(T);
T_prima = 0 : dt_r : (length(T) - 1) * dt_r;