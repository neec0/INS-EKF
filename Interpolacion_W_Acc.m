%Interpolacion lineal

function [W_b_nb_intrp, F_b_nb_intrp_gs] = Interpolacion_W_Acc (W_b_nobias, F_b_nobias_gs, T_prima)



% W_x=interp1(T,W_b_nobias(:,1),T_prima,'linear');
% W_y=interp1(T,W_b_nobias(:,2),T_prima,'linear');
% W_z=interp1(T,W_b_nobias(:,3),T_prima,'linear');
% W_b_nb_intrp=[W_x, W_y, W_z];
W_b_nb_intrp = W_b_nobias;
% 
% F_x=interp1(T,F_b_nobias_gs(:,1),T_prima,'linear');
% F_y=interp1(T,F_b_nobias_gs(:,2),T_prima,'linear');
% F_z=interp1(T,F_b_nobias_gs(:,3),T_prima,'linear');
% F_b_nb_intrp_gs=[F_x,F_y,F_z];
F_b_nb_intrp_gs = F_b_nobias_gs;