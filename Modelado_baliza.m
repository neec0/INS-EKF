

function [Taps, Senal_baliza, Baliza_detect] = Modelado_baliza (T_prima) %A Seg_del_tap se le asignara pr1_dha(:,16)

global f_Hz Ptx_Watt Grx_Watt Gtx_Watt Seg_del_tap Indice_init_Nmuestra Indice_fin_Nmuestra Num_balizas dt_r

Fin_muestras = length(T_prima);
Taps = zeros(length(T_prima),1);
Senal_baliza = zeros(length(T_prima),1 );
Baliza_detect = zeros(length(T_prima),1 );
j = 1;


%Busco taps de la captura en funcion del tiempo en que han sido detectadas
%y entonces coloco la campana
for n = 1:Fin_muestras
    
    if ( T_prima(n,:)<Seg_del_tap(j,1) )
            Taps(n,1)=0;
    else
        Taps(n,1)=1;
        
        if (j==Num_balizas)%Numero total de balizas
            Taps(n,1)=0;%Para evitar que la campana se corte eliminamos intencionadamente el ultimo tap
        else
            [Senal_baliza(n-Num_balizas:n+Num_balizas,1), Max_Prx_Watt] = Campana_Gaus(dt_r);%colocamos la campana sobre el tap
            j=j+1;%Avanzamos un tap
            Seg_del_tap(j,1) = Seg_del_tap(j-1,1) + Seg_del_tap(j,1);%Sumammos en tiempo anterior a la cuenta total
        end
        
    end
end
    
%Rellenar hueco de Taps en el espacio de tiempo de parada entre la
%vuelta 1 y la 2
for i = (Indice_init_Nmuestra + 1) : 1 : (Indice_fin_Nmuestra)
    
    Tam_Nmuestras=Indice_fin_Nmuestra-Indice_init_Nmuestra;
    Taps(i,:) = 1;
    [Senal_baliza(i,1)] = Max_Prx_Watt+((rand(1,length(Tam_Nmuestras))-0.5)*0.02)';%Punto maximo de potencia de señal con ruido sumado
    
end

%Compruebo si se ha sobrepasado el umbral minimo de deteccion
for z = 1:Fin_muestras
    
    Landa = (3*power(10,2)) / f_Hz;%Velocidad de la luz entre la frecuencia (longitud de onda)
    d = ( 4*pi / Landa ) * sqrt( Senal_baliza(z,1) / Ptx_Watt * Grx_Watt * Gtx_Watt );
    
    if ( d<1 )
        Baliza_detect(z,1) = 1;
    else
        Baliza_detect(z,1) = 0;
    end
    
end





function [Prx_Watt, Max_Prx_Watt] = Campana_Gaus (dt_r)

global f_Hz Ptx_Watt Grx_Watt Gtx_Watt V_media;

%Campana de Gaus señal de balizas Poscion/Tiempo

%Parametrizacion
Prx_Watt = zeros(53,1);
Landa = (3 * power(10,2)) / f_Hz;%Velocidad de la luz entre la frecuencia (longitud de onda)
d2_m = dt_r * V_media;%Tiempo entre muestras [s] * velocidad media de la captura [m/s] =distancia [m]
d1_m = 27 * d2_m;%27 muestras que forman la mitad positiva + la parte central de la campana
d3_m = d2_m;
paso_muestras_s = d2_m / 0.8;%Espacio entre muestras en segundos

%Subida y muestra central de la campana
fin_tramo1 = 27 * paso_muestras_s;%Paso multiplicado por 26 muestras + 1 central
j = 1;

for inicio_tramo1=0:paso_muestras_s:fin_tramo1
    
    if j==27
        Prx_Watt(27,1) = Prx_Watt(26,1);
    else
        Prx_Watt(j,1) = Ptx_Watt * Grx_Watt * Gtx_Watt * power((Landa / (4*pi*d1_m)),2);
        d1_m = d1_m-d2_m;
        j = j+1;
    end
    
end

%Bajada de la campana
fin_tramo2 = 25 * paso_muestras_s;%Paso multiplicado por 25 muestras porque empieza en cero [0:25]=26
j = 1;

for inicio_tramo2=0:paso_muestras_s:fin_tramo2
    
    if j==1
        Prx_Watt(28,1) = Prx_Watt(27,1);
    else
        d2_m = d2_m + d3_m;
        Prx_Watt(27+j,1) = Ptx_Watt * Grx_Watt * Gtx_Watt * power((0.33 / (4*pi*d2_m)),2);     
    end
    j = j+1;
    
end

Prx_Watt = Prx_Watt + ((rand(1,length(Prx_Watt)) - 0.5) * 0.02)';%Suma de ruido
Max_Prx_Watt = max(Prx_Watt); 



