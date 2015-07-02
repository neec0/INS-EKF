
function [Step_detect, Modulo_Acc, Desv_est, Modulo_W, Step_original] = Detec_paso(W_b_nbdet_rads, F_b_nbdet_gs)

global Fin_muestras gravedad
Step_detect = zeros([Fin_muestras,1]);%Defino un array de ceros para mas tarde rellenarlo
Modulo_Acc = zeros([Fin_muestras,1]);%Defino un array de ceros para mas tarde rellenarlo
Modulo_W = zeros([Fin_muestras,1]);%Defino un array de ceros para mas tarde rellenarlo
C_1 = zeros([Fin_muestras,1]);
C_2 = zeros([Fin_muestras,1]);
C_3 = zeros([Fin_muestras,1]);
Desv_est = zeros([Fin_muestras,1]);
Step_original = zeros([Fin_muestras,1]);

%Condicion 1 para detectar pasos. El modulo de la aceleracion debe estar en
%un umbral

for i=1:Fin_muestras

    %Modulo de la aceleracion
    Suma_cuadrados_Acc = sum(((F_b_nbdet_gs(i,:) * gravedad).^2),2);%Suma de los cuadrados de las acerleraciones
    Modulo_Acc(i,1) = sqrt(Suma_cuadrados_Acc);%Raiz cuadrada de la suma
    
    %Condicion 1. El modulo de la aceleracion debe estar un umbral
    if((Modulo_Acc(i,:) > 9) && (Modulo_Acc(i,:) < 10.5))
        C_1(i,1)=1;
    else
        C_1(i,1)=0;
    end

    %Condicion 2. La varianza de aceleracion local debe estar por necima de en un umbral
    sample_leng = 16;%Ventana de diez muestras para que al multiplicar sean 20 (50Hz)
    if (i <= sample_leng)
        C_2(i,1) = 1;
    else
        Desv_est(i,1) = var(Modulo_Acc((i-sample_leng)+1:i));
        if (Desv_est(i,1) < 0.5)%0.5
            C_2(i,1) = 1;
        else
            C_2(i,1)=0;
        end
    end

    %Modulo de la velocidad
    Suma_cuadrados_W = sum((W_b_nbdet_rads(i,:).^2),2);%Suma de los cuadrados de las velociadades
    Modulo_W(i,1) = sqrt(Suma_cuadrados_W);%Raiz cuadrada de la suma
    
    %Condicion 3 para detectar pasos. El modulo de la velocidad angular debe
    %estar en dos umbrales.
    if((Modulo_W(i,:) >= 0) && (Modulo_W(i,:) < 0.2))%0.2
        C_3(i,1)=1;
    else
        C_3(i,1)=0;
    end
    
    if ((C_1(i,1)*C_2(i,1)*C_3(i,1)) == 1)
        Step_detect(i,1)=1;
    else
        Step_detect(i,1)=0;
    end
   
    Step_original(i,1) = Step_detect(i,1);
    
    [Step_detect] = Filtro_mediana(sample_leng, Step_detect, i);
    
end

Step_detect(7:length(Step_detect),1)=Step_detect(1:(length(Step_detect)-6),1);%La ventana de muestras de la varianza es de 5 y
Step_detect(1:6,1)=1;%la de el filtro es de 11. En el filtro se toman 5 muestras pasadas, 5 futuras y la actual del bucle (i)


figure, plot(Modulo_Acc, 'b')
hold on;
plot(Modulo_W, 'k')
% plot(sample_leng:(sample_leng + length(Desv_est) - 1),Desv_est, 'r')
plot(Desv_est, 'r')
plot(Step_detect, 'g')
title('Detección de pasos y condiciones')
grid MINOR
hold off;

%Filtro de mediana
function [Step_detect] = Filtro_mediana (sample_leng, Step_detect, i)
    ventana = (1.5 * sample_leng);
    if (i > (ventana))%La ventana la forman 11 muestras
        j = (i - ventana) + 1;
        Step_detect(j,1) = median(Step_detect(j:i,1),1);%Actualizacion despues de la mediana de la ventana
    end%Despues de dejar pasar las 11 primeras muestras (ventana) el valor filtrado se guarda empezando desde la primera muestra