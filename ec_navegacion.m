
function [ Pos, delta_X ] = ec_navegacion (W_b_nbdet_rads, F_b_nbdet_gs, Ai_rad, T_prima, Step_detect)

global Baliza_detect
global Bias_W_rads Bias_F_ng_ms2 C
global gravedad dt_r
global k m0 mI Fin_muestras T0 contMuestras
global yaw YawInSD prevYawToSD

Pos = zeros(Fin_muestras,10);%Defino un array de ceros para mas tarde rellenarlo
X_r_m = zeros(Fin_muestras,3);
V_r_ms = zeros(Fin_muestras,3);
% Celdas = zeros([length(F_b_nbdet_gs),5]);
delta_X = zeros(Fin_muestras, 15);
delta_Yaw = zeros(Fin_muestras, 1);
W_b_nbdet_correct_rads = zeros(Fin_muestras, 3);


for k=1:Fin_muestras
        
%         if ( k~=1 )
%             W_b_nbdet_correct_rads(k,:) = W_b_nbdet_correct_rads(k,:) - delta_X(k-1,4:6);
%             F_b_nb_intrp_ms2(k,:) = F_b_nb_intrp_ms2(k,:) - delta_X(k-1,13:15);
%         end
        
    if (k == 1)
%             C_0 = [             1 -Ai_rad(1, 3)  Ai_rad(1, 2); 
%                      Ai_rad(1, 3)             1 -Ai_rad(1, 1); 
%                     -Ai_rad(1, 2)  Ai_rad(1, 1)             1 ];

        r = Ai_rad(1);%Inicializo con aprox. allingment el roll
        p = Ai_rad(2);%Inicializo con aprox. allingment el pitch
        yaw(k,1) = Ai_rad(3);%Inicializo con aprox. allingment el yaw
%         delta_Yaw(k,1) = Ai_rad(3);
        
        C1 = [  cos(yaw(k,1)) sin(yaw(k,1)) 0;
               -sin(yaw(k,1)) cos(yaw(k,1)) 0;
                            0             0 1 ];

        C2 = [ cos(p) 0 -sin(p);
                    0 1      0;
               sin(p) 0 cos(p)];

        C3 = [ 1       0      0;
               0  cos(r) sin(r);
               0 -sin(r) cos(r) ];

        C1 = C1';
        C2 = C2';
        C3 = C3';

        C = C1 * C2 * C3;


        %             C_s(k,:)=[C(3,2),C(3,3)];
        C_0 = C;
        V_r0_ms = zeros(1,3);
        X_r0_m = zeros(1,3);
        
        contMuestras = 0;
        Step_detectAnt = 0;
        
        YawInSD = 0;
        prevYawToSD = 0;
    else
        
        W_b_nbdet_correct_rads(k,:) = W_b_nbdet_rads(k,:) - delta_X(k-1,4:6);
        
        
        B_t = [                                 0 -W_b_nbdet_correct_rads(k,3)*dt_r  W_b_nbdet_correct_rads(k,2)*dt_r; 
                 W_b_nbdet_correct_rads(k,3)*dt_r                                 0 -W_b_nbdet_correct_rads(k,1)*dt_r; 
                -W_b_nbdet_correct_rads(k,2)*dt_r  W_b_nbdet_correct_rads(k,1)*dt_r                                 0 ];
            
        %MÉTODO 1 matriz de orientación (a veces genera números imaginarios)------------------------------------
%         O = sqrt( power((W_b_nbdet_correct_rads(k,1)*dt_r),2) + power((W_b_nbdet_correct_rads(k,2)*dt_r),2) + power((W_b_nbdet_correct_rads(k,3)*dt_r),2) );
%         C = C_0 * (eye(3) + (sin(O)/O)*B_t + ((1-cos(O))/power(O,2))*power(B_t,2));
        %Fin---------------------------------------------------------------
        
        
        %MÉTODO 2 matriz orientación (sin aproximación)--------------------
%         C = C_0 * expm(B_t);
        %Fin---------------------------------------------------------------
        
        
        %MÉTODO 3 matriz orientación (aproximación Padé)-------------------
        C = C_0 * (2 * mI + B_t) / (2 * mI - B_t);
        %Fin---------------------------------------------------------------

        r = atan2(C(3,2),C(3,3));
        p = asin(-C(3,1));
        yaw(k,1) = atan2(C(2,1),C(1,1));
    end
    
    
    
    %ECUACIONES NAVEGACION-------------------------------------------------
    if(k ~= 1)
        F_b_nbdet_gs(k,:) = F_b_nbdet_gs(k,:) - (delta_X(k-1,13:15) / gravedad);%[1X3] - [1X3] = [1X3]
    end


    Acc_r_nbdet_ms2 = (C * F_b_nbdet_gs(k,:)')' * gravedad;%([3X3] * [1X3]')' * gravedad = [1X3] Pasar al marco de refrencia yaw a m/s^2
    Acc_r_nbdet_ng_ms2 = Acc_r_nbdet_ms2 + [ 0 0 gravedad ];%Sumo g para que compense la gravedad en el eje z

    
    V_r_ms(k,:) = V_r0_ms + (dt_r * Acc_r_nbdet_ng_ms2);%[1X3]
    X_r_m(k,:) = X_r0_m + (dt_r * V_r_ms(k,:));%[1X3]
    
    if( ((Step_detectAnt == 0) && (Step_detect(k,1) == 0.5)) || ((Step_detectAnt == 1) && (Step_detect(k,1) == 1)))
        contMuestras = contMuestras + 1;
%         fprintf('%d\n', contMuestras);
    elseif((Step_detectAnt == 0.5) && (Step_detect(k,1) == 0))
        contMuestras = 0;
    end
    
    if(Step_detect(k,1))
        [ delta_X(k,:) ] = EKF ( W_b_nbdet_rads(k,:), F_b_nbdet_gs, V_r_ms(k,:), Step_detect );%[1X15]
    end   
    
    if(Step_detect(k,1))
        V_r_ms(k,:) = V_r_ms(k,:) - delta_X(k,10:12);
        X_r_m(k,:) = X_r_m(k,:) - delta_X(k,7:9);
    end
    V_r0_ms = V_r_ms(k,:);
    X_r0_m = X_r_m(k,:);
    %Fin de ecuaciones de navegacion---------------------------------------
    
    
    %Update MATRIZ ORIENTACION---------------------------------------------
    delta_O = [             0 -delta_X(k,3)  delta_X(k,2);
                 delta_X(k,3)             0 -delta_X(k,1);
                -delta_X(k,2)  delta_X(k,1)             0 ];
            
    C_update = ((2 * mI + delta_O) / (2 * mI - delta_O)) * C;
    C_0 = C_update;
    %Fin-------------------------------------------------------------------
    
    Pos(k,:) = [ X_r_m(k,1), X_r_m(k,2), X_r_m(k,3), r, p, yaw(k,1), T_prima(k), V_r_ms(k,1), V_r_ms(k,2), V_r_ms(k,3)];  
    
%     if((k == 1) && (Step_detectAnt == 0) && Step_detect(k,1))
%         YawInSD = yaw(k,1);
%         prevYawToSD = yaw(k,1);
    if((Step_detectAnt == 0) && (Step_detect(k,1) == 0.5))
        YawInSD = yaw(k-1,1);
        prevYawToSD = yaw(k-2,1);
%         fprintf('Flanco ascendente YawSD: %.02f YawAntSD: %0.2f', YawInSD, prevYawToSD);
    end
    
    Step_detectAnt = Step_detect(k,1);
end
figure, plot(rad2deg(W_b_nbdet_correct_rads(:,1)), 'b');
hold on;
plot(rad2deg(W_b_nbdet_correct_rads(:,2)), 'k');
plot(rad2deg(W_b_nbdet_correct_rads(:,3)), 'r');
xlabel('muestras');
ylabel('Velocidad angular [º/s]');
title('Velocidad angular BIAS determinista & BIAS EKF');
hold off;
figure, plot(F_b_nbdet_gs(:,1)*gravedad, 'b');
hold on;
plot(F_b_nbdet_gs(:,2)*gravedad, 'k');
plot(F_b_nbdet_gs(:,3)*gravedad, 'r');
xlabel('muestras');
ylabel('Aceleracion [m/s^2]');
title('Aceleracion BIAS determinista & BIAS EKF');
hold off;