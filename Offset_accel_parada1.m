

function [A_noff] = Offset_accel_parada1 (A, Angles)

r=Angles(1);
p=Angles(2);
y=0;

C1=[cos(y) sin(y) 0;
    -sin(y) cos(y) 0;
    0 0 1];

C2=[cos(p) 0 -sin(p);
    0 1 0;
    sin(p) 0 cos(p)];

C3=[1 0 0;
    0 cos(r) sin(r);
    0 -sin(r) cos(r)];

C1=C1';
C2=C2';
C3=C3';

C=C1*C2*C3;

% Usamos C para rotar las aceleraciones durante la parte estática y poder
% restar la contribución de la gravedad antes de calcular los bias

N = 600;
Arot = zeros(size(A));
Cinv = inv(C);
for i=1:N
    Arot(i,:) = (C * A(i,:)')' + [0,0,1];
    Arot(i,:) = Arot(i,:) * Cinv;
end

Ax_mean = mean(Arot(1:N,1));
Ay_mean = mean(Arot(1:N,2));
Az_mean = mean(Arot(1:N,3));

Ax_noff = Arot(:, 1)-Ax_mean;
Ay_noff = Arot(:, 2)-Ay_mean;
Az_noff = Arot(:, 3)-Az_mean;
A_noff = [Ax_noff, Ay_noff, Az_noff];

figure, plot(Arot(:,1:3) * 9.8)


