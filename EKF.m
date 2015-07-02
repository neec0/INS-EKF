
function [ delta_X ] = EKF ( W_b_nb_intrp, F_b_nbdet_gs, V_r_NoCorrect_ms, Step_detect )

global delta_X_anterior matriz_FI P P_0 K Q R H m
global dt_r Baliza_detect C
global gravedad m0 mI T0
global k Fin_muestras contMuestras
global yaw YawInSD prevYawToSD


F_r_nbdet_gs(k,:) = C * F_b_nbdet_gs(k,:)';


%DECLARACIÓN E INICIALIZACIÓN VARIABLES EKF--------------------------------
if (k == 1)%Inicialización de variables en la primera iteración
    
    %Vector de estado de error [1X15]
    delta_X = zeros(1,15);%Declaracion del tamaño, cada iteración aumenta
    Vector_X_0 = [ 10^-6 10^-6 10^-6 ];
    m0_1X3 = [ 0 0 0 ];
    delta_X_anterior = [ m0_1X3 m0_1X3 m0_1X3 m0_1X3 m0_1X3 ];%Inicialización
    %Fin declaración vector estado de error
    
    %Declaración matriz P [Fin_muestrasXFin_muestras]
    P = zeros(15);
    %Fin declaración matriz P

    %Matriz P_0, donde se guarda P de iteración anterior [15X15]   
    P_0_w0 = [ 1*10^-2       0       0;%[rad/s]
                     0 1*10^-2       0;
                     0       0 1*10^-2 ];
    P_0_a0 = [ 1*10^-2       0       0;%[m/s^2]
                     0 1*10^-2       0;
                     0       0 1*10^-2 ];

%     P_0 = [ m0     m0 m0 m0 m0;%Matriz P anterior evaluada en el TU
%             m0 P_0_w0 m0 m0 m0;
%             m0     m0 m0 m0 m0;
%             m0     m0 m0 m0 m0;
%             m0     m0 m0 m0 m0 ];

    P_0 = [ m0     m0 m0 m0     m0;%Matriz P anterior evaluada en el TU
            m0 P_0_w0 m0 m0     m0;
            m0     m0 m0 m0     m0;
            m0     m0 m0 m0     m0;
            m0     m0 m0 m0 P_0_a0 ];
    %Fin declaración matriz P
    
    %Declaracion de la matriz de covarianza de error del proceso (sensores) [15X15]
    ruido_proceso_fi = [ 1*10^-4       0       0;
                               0 1*10^-4       0;
                               0       0 1*10^-4 ];

    ruido_proceso_vel = [ 1*10^-4       0       0;
                               0 1*10^-4       0;
                               0       0 1*10^-4 ];

    Q = [ ruido_proceso_fi m0 m0                m0 m0;
                        m0 m0 m0                m0 m0;
                        m0 m0 m0                m0 m0;
                        m0 m0 m0 ruido_proceso_vel m0;
                        m0 m0 m0                m0 m0 ];
    %Fin declraracion matriz de covarianza de error del proceso
    
    %Declaracion de la matriz de covarianza de medida del error de la medida [6X6]
    R_zaru = [ 0.1   0   0;%Vector de tuning para el ZARU [rad/s]
                 0 0.1   0;
                 0   0 0.1 ];
             
    R_HDR = 0.1;%[1X1]
             
    R_zupt = [ 0.01    0    0;%Vector de tuning para el ZUPT [m/s]
                  0 0.01    0;
                  0    0 0.01 ];

%     R = [ R_zaru     m0;
%               m0 R_zupt ];
          
    R = [ R_HDR   [0 0 0] [0 0 0];
          [0;0;0]  R_zaru      m0;
          [0;0;0]      m0  R_zupt ];
    %Fin de declaración de matriz de covarianza de error de la medida
    
%     %Declaración de la matriz de medida [6X15]
%     H = [ m0 mI m0 m0 m0;
%           m0 m0 m0 mI m0 ];
%     %Fin declaracion matriz de medida
    
    %Declaración de la matriz de medida [7X15]
    H = [ [0 0 1] [0 0 0] [0 0 0] [0 0 0] [0 0 0]
               m0      mI      m0      m0      m0;
               m0      m0      m0      mI      m0 ];
    %Fin declaracion matriz de medida
end

%Inicialización de variables para todas las iteraciones
%Declaracion matriz sesgada de acelereaciones: "inclometer" [3X3]
S_A_r_nbdet_ms2 = [                0 -F_r_nbdet_gs(3)  F_r_nbdet_gs(2);
                     F_r_nbdet_gs(3)                0 -F_r_nbdet_gs(1);
                    -F_r_nbdet_gs(2)  F_r_nbdet_gs(1)                0 ];
%Fin declaracion matriz sesgada de acelereaciones: "inclometer"

%Declaración de la matriz Delta simetrica sesgada [15X15 = 5X5 vectores de 3X3]
matriz_FI = [                             mI dt_r*C m0      m0     m0;
                                          m0     mI m0      m0     m0;
                                          m0     m0 mI dt_r*mI     m0;
              -dt_r*S_A_r_nbdet_ms2*gravedad     m0 m0      mI dt_r*C;
                                          m0     m0 m0      m0     mI ]; 
%Fin declaracion matriz Delta
%FIN DECLARACIÓN E INICIALIZACIÓN------------------------------------------



%RESETEO DEL VECTOR DE ERROR DE ESTADO TRAS DETECCIÓN DE PASO--------------
if ((k ~= 1) && (Step_detect(k-1,:) == 1))%Reseteo del vector de estado de error excepto los bias para que converjan
   delta_X_anterior = delta_X_anterior .* [ 0 0 0 1 1 1 0 0 0 0 0 0 1 1 1 ];
end
%FIN RESETEO---------------------------------------------------------------



%TIME UPDATE---------------------------------------------------------------
delta_X = (matriz_FI * delta_X_anterior')';%(+ Wk)[1X15]
delta_X_anterior = delta_X;

P = matriz_FI * P_0 * matriz_FI' + Q;%[15X15]
P_0 = P;%Actualización de la matriz P
%--------------------------------------------------------------------------
    


%MEASURMENT UPDATE---------------------------------------------------------
    %ZARU [rad/s]
    if(contMuestras >= 160)%Si el paso dura mas de 180 muestras (aprox. el doble de la media de un paso)
        DELTA_W_b = W_b_nb_intrp - [ 0 0 0 ];%Viandante parado
    else
        DELTA_W_b = [ 0 0 0 ];%Viandante en movimiento
    end
    %Fin ZARU

    %ZUPT [m/s^2]
    DELTA_v = V_r_NoCorrect_ms - [ 0 0 0 ];%El pie esta en el suelo, la velocidad es 0
    %Fin ZUPT
    
    %HDR [º]
    if((YawInSD ~= 0.0) && (prevYawToSD ~= 0.0))
        DELTA_Yaw = yaw - ((YawInSD + prevYawToSD) / 2);
        Abs_DELTA_Yaw = abs(DELTA_Yaw);
%         fprintf('Abs Yaw: %.02f\n', (abs(rad2deg(DELTA_Yaw))));
        if(Abs_DELTA_Yaw >= 4)
%             fprintf('Giro\n');
%             fprintf('DELTA_Yaw > 4: %0.2f Muestra k: %d\n', (abs(rad2deg(DELTA_Yaw))), k);
        else
%             fprintf('DELTA_Yaw < 4: %.02f Muestra k: %d\n', (abs(rad2deg(DELTA_Yaw))), k);
%             fprintf('YawAnt: %.02f YawAntAnt: %.02f\n', YawInSD, prevYawToSD);
%             fprintf('Pasillo\n');
            DELTA_Yaw = 0;
        end
    else
        DELTA_Yaw = 0;
    end
    %Fin HDR

%     %Declaración del vector error de mediad m [1X6 traspuesto]
%     m = [ DELTA_W_b, DELTA_v ];%Error W y error v
%     %Fin declaración del vector de medida m
    
    %Declaración del vector error de mediad m [1X7 traspuesto]
    m = [ DELTA_Yaw, DELTA_W_b, DELTA_v ];%Error Y, error W y error v
    %Fin declaración del vector de medida m
    
    K = P * H' * pinv( H * P * H' + R );%[15X6]

    delta_X = (delta_X' + K * (m' - H * delta_X'))';%Vector de estado de error en el MU [1X15]
    delta_X_anterior = delta_X;
    
    P = ( eye(15) - K * H ) * P_0;% Vector P evaluado en el MU [15X15]
    P_0 = P;%Actualizacion del vector P con la medida del MU [15X15]
    %Fin measurement update
%--------------------------------------------------------------------------
