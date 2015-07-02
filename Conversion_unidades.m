
function [W_b_rad, F_b_gs] = Conversion_unidades ()

global Matriz_captura

W_b_deg = [ Matriz_captura(:, 4), Matriz_captura(:, 6), Matriz_captura(:, 8) ];
W_b_rad = deg2rad(W_b_deg);%Conversion de vel. angular a rad/s

F_b_gs=[ Matriz_captura(:, 10), Matriz_captura(:, 12), Matriz_captura(:, 14) ];%Aceleracion expresada en gs




