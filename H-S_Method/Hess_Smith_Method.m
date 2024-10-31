%% Perfiles
% Naca 1216 

clear all;
close all;
clc;

n_nodos = 141;              % Numero de nodos
n_paneles = n_nodos - 1;    % Numero de Paneles

f = 0.01;       % Curvatura
xf = 0.2;      % Curvatura Máxima
t = 0.16;       % Espesor

angulo = 0.503;         % Grados
alfa = angulo/180*pi;   % Radianes


%% Calculo perfil

% Barrido del angulo
beta = zeros(1, n_nodos);
x = zeros(1, n_nodos);
z = zeros(1, n_nodos);
zc = zeros(1, n_nodos);
zt = zeros(1, n_nodos);

for i = 2: n_nodos
    beta(i) = beta(i - 1) + 2*pi/(n_nodos - 1);
end

for i = 1:n_nodos

    x(i) = 0.5*(1 + cos(beta(i)));

    if x(i) < xf 
        zc(i) = f*(2*xf*x(i) - x(i)^2)/(xf^2);
    else 
        zc(i) = f*((1-2*xf) + 2*xf*x(i) - x(i)^2)/((1 - xf)^2);
    end

    zt(i) = 5*t*(0.2969*sqrt(x(i)) - 0.126*x(i) - 0.3516*x(i)^2 + 0.2843*x(i)^3 - 0.1036*x(i)^4);

    if beta(i)<pi

        z(i) = zc(i) - zt(i);

    else

        z(i) = zc(i) + zt(i);

    end
end

% Corregir Naca y cerrar paneles
z(1) = 0;
z(n_nodos) = z(1);

%% Calculo puntos control y paneles

ddx = zeros(1, n_paneles);
ddz = zeros(1, n_paneles);

% Declaracion paneles 
theta = zeros(1, n_paneles);
xpc = zeros(1, n_paneles);
zpc = zeros(1, n_paneles);
c = zeros(1, n_paneles);

for i =1:n_paneles

    ddx(i) = x(i + 1) - x(i);

    ddz(i) = z(i + 1) - z(i);

    theta(i) = atan2(ddz(i),ddx(i));

    xpc(i) = (x(i) + x(i + 1))/2;

    zpc(i) = (z(i) + z(i + 1))/2;

    c(i) = abs(ddx(i)*cos(theta(i))) + abs(ddz(i)*sin(theta(i)));
    
   
end


%% Coordenadas en ejes panel

xpc_ij = zeros(n_paneles);
zpc_ij = zeros(n_paneles);
theta1 = zeros(n_paneles);
theta2 = zeros(n_paneles);
r1 = zeros(n_paneles);
r2 = zeros(n_paneles);
upcm = zeros(n_paneles); % velocidad iducida en ejes panel
wpcm = zeros(n_paneles);
upct = zeros(n_paneles);
wpct = zeros(n_paneles);

for i=1:n_paneles
    for j=1:n_paneles
        xpc_ij(i,j) = (xpc(i) - x(j))*cos(theta(j)) + (zpc(i) - z(j))*sin(theta(j));
        zpc_ij(i,j) = - (xpc(i) - x(j))*sin(theta(j)) + (zpc(i) - z(j))*cos(theta(j));

        theta1(i,j) = atan2(zpc_ij(i,j),xpc_ij(i,j));

        theta2(i,j) = atan2(zpc_ij(i,j), (xpc_ij(i,j)) - c(j));

        r1(i,j)=((xpc_ij(i,j))^2 + (zpc_ij(i,j))^2)^0.5;
        r2(i,j)=((xpc_ij(i,j) - c(j))^2 + (zpc_ij(i,j))^2)^0.5;

        upcm(i,j) = -1/(2*pi)*log(r2(i,j)/r1(i,j));
        wpcm(i,j) = (theta2(i,j) - theta1(i,j))/(2*pi);

        upct(i,j) = (theta2(i,j) - theta1(i,j))/(2*pi);
        wpct(i,j) = 1/(2*pi)*log(r2(i,j)/r1(i,j));
    end
end

%% En ejes perfil

upcm2 = zeros(n_paneles);
wpcm2 = zeros(n_paneles);
upct2 = zeros(n_paneles);
wpct2 = zeros(n_paneles);

for i=1:n_paneles
    for j=1:n_paneles
        

        upcm2(i,j) = upcm(i,j)*cos(theta(j)) - wpcm(i,j)*sin(theta(j));
        wpcm2(i,j) = upcm(i,j)*sin(theta(j)) + wpcm(i,j)*cos(theta(j));

        upct2(i,j) = upct(i,j)*cos(theta(j)) - wpct(i,j)*sin(theta(j));
        wpct2(i,j) = upct(i,j)*sin(theta(j)) + wpct(i,j)*cos(theta(j));
    end
end

%% Proyeccion de la velocidad sobre la superfice

A = zeros(n_nodos);
I = zeros(n_nodos, 1);

for i=1:n_paneles
    for j=1:n_paneles
        
        A(i,j)= -upcm2(i,j)*sin(theta(i)) + wpcm2(i,j)*cos(theta(i));
        A(i,n_nodos) = A(i,n_nodos) - upct2(i,j)*sin(theta(i)) + wpct2(i,j)*cos(theta(i));

    end
        A(i,i) = abs(A(i,i));
    I(i)= (cos(alfa)*sin(theta(i)) - sin(alfa)*cos(theta(i)));
end

%% Condición de Kutta

for j = 1:n_paneles
    A(n_nodos, j) = upcm2(1, j)*cos(theta(1)) + wpcm2(1, j)*sin(theta(1)) + upcm2(n_paneles,j)*cos(theta(n_paneles)) + wpcm2(n_paneles,j)*sin(theta(n_paneles));
    A(n_nodos, n_nodos) = A(n_nodos, n_nodos) + upct2(1, j)*cos(theta(1)) + wpct2(1, j)*sin(theta(1)) + upct2(n_paneles,j)*cos(theta(n_paneles)) + wpct2(n_paneles,j)*sin(theta(n_paneles));
end

I(n_nodos) = -(cos(alfa - theta(1)) + cos(alfa - theta(n_paneles)));

%% Resolver el sistema

B = zeros(n_paneles, n_nodos);
tan = zeros(1, n_paneles);

for i=1:n_paneles
    for j=1:n_paneles
        
        B(i,j) = upcm2(i,j)*cos(theta(i)) + wpcm2(i,j)*sin(theta(i));
        B(i, n_nodos) = B(i, n_nodos) + abs(upct2(i,j)*cos(theta(i)) + wpct2(i,j)*sin(theta(i)));
    end
    aux(i) = B(i, n_nodos);
    tan(i) = cos(alfa)*cos(theta(i)) + sin(alfa)*sin(theta(i));
end




vtan = zeros(1, n_paneles);

INV = A^(-1);


sol_q = INV*I;



for i=1:n_paneles
    for j=1:n_paneles
        vtan(i) = vtan(i) + B(i,j)*sol_q(j);
    end
     vtan(i) = vtan(i) + B(i,n_nodos)*sol_q(n_nodos) + tan(i);
end

vtan(n_paneles) = vtan(1);
Cmc4_paneles = 0;

cp = zeros(1, n_paneles);
cl_panel = zeros(1, n_paneles);
CL_panel = 0;
for i = 1:n_paneles
    cp(i)= 1 - (vtan(i))^2;
    cl_panel(i) = -cp(i)*c(i)*cos(theta(i) - alfa);
    CL_panel = CL_panel + cl_panel(i);

    Cmc4_paneles = Cmc4_paneles - cl_panel(i)*(xpc(i) - 0.25)*cos(alfa) - cl_panel(i)*zpc(i)*sin(alfa);
end

C_L = 0;
Cm0 = 0;
Cmc4 = 0;

%% Representar resultados

n_resultados = n_paneles/2;
for i=1:n_resultados
    
    xcl(i) = 0.5*(1 + cos(i*pi/n_resultados));

    if xcl<xf
        zccl(i) = f*(2*xf*x(i) - x(i)^2)/(xf^2);
    else 
        zccl(i) = f*((1-2*xf) + 2*xf*x(i) - x(i)^2)/((1 - xf)^2);
    end

    ztcl(i)= 5*t*(0.2969*x(i)^0.5-0.1260*x(i)-0.3516*x(i)^2+0.2843*x(i)^3-0.1015*x(i)^4);

    if i<n_resultados/2
        zcl(i) = zccl(i) - ztcl(i);
    else
        zcl(i) = zccl(i) + ztcl(i);
    end

    Cl(i) = cp(i) - cp(n_paneles - i );


    Cmc4 = Cmc4 - Cl(i)*(xcl(i) - 0.25)*cos(alfa)*ddx(2*i) - Cl(i)*zcl(i)*ddz(2*i)*sin(alfa);

    Cm0 = Cm0 - Cl(i)*ddx(2*i)*xcl(i)*cos(alfa) - Cl(i)*ddz(2*i)*zcl(i)*sin(alfa);

    C_L = C_L + Cl(i)*abs(ddx(2*i));
    
end



figure(1)
plot(xpc,cp)
ylabel('Cp(x) intradós y extradós')
xlabel('x')
ay =gca;
ay.YDir = 'reverse';
grid

figure(9)
plot(xpc,vtan)
ylabel('velocidad tangente')
xlabel('x')
ay =gca;
ay.YDir = 'reverse';
grid

figure(11)
plot(xpc,cl_panel)
ylabel('Cl_intrados y estrados')
xlabel('x')
ay =gca;
ay.YDir = 'reverse';
grid


figure(2)
plot(xcl, Cl);

ylabel('Cl_resta(x)')
xlabel('x')
grid

figure(3)
hold on;

plot(x,z,'Color','#000000')
plot(x,z,'o','Color','#0000FF')
plot(xpc,zpc,'s','Color','#FF0000')

%axis proportional
xlim([-0.2 0.2])
xlabel('Perfil NACA')
grid
hold off

legend('Contorno','Nodos','Puntos de Control')

%% Ayuda paneles

% theta 1_ij
for i = 1:n_paneles
    theta1_plot(i) = theta1(12,i);
    theta2_plot(i) = theta2(12,i);
    num(i) = i;
end
figure(4)
plot(num,theta1_plot)
xlabel('theta1 plot')
grid
hold on;

plot(num, theta2_plot)

%theta
figure(5)
plot(num, theta)
xlabel('theta plot')
grid

% x punto control ij
for i = 1:n_paneles
    x_plot(i) = xpc_ij(12,i);
    z_plot(i) = zpc_ij(12,i);
end

figure(6)
plot(num,x_plot)
xlabel('x z plot')
grid

hold on;

plot(num, z_plot)

% R1
for i = 1:n_paneles
    r1_plot(i) = r1(40,i);
    r2_plot(i) = r2(40,i);
end

figure(7)
plot(num,r1_plot)

grid

hold on;

plot(num, r2_plot)

% q
for i = 1:n_paneles
    q_plot(i) = sol_q(i);
end

figure(8)
plot(xpc,q_plot)
xlabel('qi plot')
grid

figure(10)
plot(xpc,aux);
ay =gca;
ay.YDir = 'reverse';
xlabel('Bi N+1 plot');

grid

%%Exportar

cp_x = zeros(2,n_paneles);

cp_x(1,:) = xpc;
cp_x(2,:) = cp;

fileID = fopen('Cp_paneles.txt', 'w');
fprintf(fileID,'%6s %12s\n','xpc','Cp');
fprintf(fileID,'%6.4f  %6.4f \r\n', cp_x);



