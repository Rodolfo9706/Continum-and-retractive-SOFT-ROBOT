clc; clear; close all;

%% ----------------- TRAYECTORIA CIRCULAR -----------------
R = 40;                  % radio del círculo (mm)
omega = 3*pi/20;         % velocidad angular (rad/s)
Xc = -32.35;             % centro X
Yc = 37.63;              % centro Y
Zc = 75.20;              % altura Z
dt = 0.01;               % paso de simulación
Tsim = 17;               % tiempo total (s)
max_iter = round(Tsim/dt);

%% ----------------- INICIALIZACIÓN -----------------
L = [60 60 60 60 35 50]'; 
Xprev = zeros(3,1);      % para derivada numérica
Kp = 50; Kd = 0.03;       % ganancias PD
max_dL = 8;             % límite de velocidad tendones
tol_pos = 0.01;

hist_X = zeros(3,max_iter);
hist_L = zeros(6,max_iter);

figure(1); hold on; grid on; axis equal;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('Soft Arm con control cartesiano PD');
view(3);



%% ----------------- LOOP DE CONTROL -----------------
hist_Xdes = zeros(3,max_iter);
hist_Xcur = zeros(3,max_iter);

for k = 1:max_iter
    t = (k-1)*dt;
    
    % Posición deseada (círculo)
    Xdes = [Xc + R*cos(omega*t);
            Yc + R*sin(omega*t);
            Zc];

    % Cinemática directa
    [Xf, Yf, Zf] = forward_end(L');
    Xcur = [Xf; Yf; Zf];
    
    % Error cartesiano
    e = Xdes - Xcur;
    
    % Velocidad cartesiana
    if k==1, Xprev = Xcur; end
    Xdot = (Xcur - Xprev)/dt;
    
    % Jacobiano numérico
    J = jacobian_numeric(L);
    
    % Ley de control PD cartesiano
    Ldot = Kp * J' * e - Kd * J' * Xdot;
    
    % Saturar velocidad
    Ldot = max(min(Ldot, max_dL), -max_dL);
    
    % Integración Euler
    L = L + Ldot*dt;
    
    % Guardar históricos
    hist_X(:,k) = Xcur;
    hist_L(:,k) = L;
    hist_Xdes(:,k) = Xdes;
    hist_Xcur(:,k) = Xcur;
    
    Xprev = Xcur;
    
    % Graficar brazo y trayectoria cada 0.05 s
    if mod(k,round(0.05/dt))==0
        softArmForward(L);

        % --- Trayectorias ---
        plot3(Xdes(1), Xdes(2), Xdes(3), 'ro', 'MarkerFaceColor','r'); % punto deseado actual
        plot3(Xcur(1), Xcur(2), Xcur(3), 'bo', 'MarkerFaceColor','b'); % efector actual
        
        plot3(hist_Xdes(1,1:k), hist_Xdes(2,1:k), hist_Xdes(3,1:k), '--r', 'LineWidth', 1.2); % trayectoria deseada (círculo)
        plot3(hist_Xcur(1,1:k), hist_Xcur(2,1:k), hist_Xcur(3,1:k), '-b', 'LineWidth', 1.5);  % trayectoria seguida
    
        legend({'','','','','','','','','','','','','','','','','','Current position','Desired trayectory'}, 'Location','best');
        drawnow;
    end
    
    % Parada por error
    if norm(e)<tol_pos
        fprintf('Error < %.2f mm en t=%.2f s\n', tol_pos, t);
        break;
    end
end


%% ----------------- FUNCIONES AUXILIARES -----------------
function [Xf,Yf,Zf] = forward_end(L)
    % Cinemática directa (usando tu función anterior)
    L1=L(1); L2=L(2); L3=L(3);
    L4=L(4); L5=L(5); L6=L(6);

    % ---- Primer segmento ----
    Lc = (L1+L2+L3)/3; r=11.4;
    DL1 = Lc-L1; DL2 = Lc-L2; DL3 = Lc-L3;
    Xdir = DL1 + DL2*cos(2*pi/3) + DL3*cos(4*pi/3);
    Ydir = DL2*sin(2*pi/3) + DL3*sin(4*pi/3);
    alfa = atan2(Ydir, Xdir);
    rp1 = r*cos(alfa);
    R1 = (Lc*rp1)/(Lc-L1+1e-6);
    B1 = Lc/R1;
    cz3d = R1*sin(B1);
    cx3d = R1*(1-cos(B1))*cos(alfa);
    cy3d = R1*(1-cos(B1))*sin(alfa);
    

    % ---- Segundo segmento ----
    Lc2 = (L4+L5+L6)/3; r=15; d=1/3;
    DL4 = Lc2-L4; DL5 = Lc2-L5; DL6 = Lc2-L6;
    X2dir = DL4*cos(d*pi) + DL5*cos(2/3*pi + d*pi) + DL6*cos(4/3*pi + d*pi);
    Y2dir = DL4*sin(d*pi) + DL5*sin(2/3*pi + d*pi) + DL6*sin(4/3*pi + d*pi);
    alfa2 = atan2(Y2dir, X2dir);
    rp4 = r*cos(alfa2 - d*pi);
    R4 = (Lc2*rp4)/(Lc2-L4+1e-6);
    B4 = Lc2/R4;
    cz2 = R4*sin(B4); cy2 = R4*(1-cos(B4)); cx2 = cy2*cos(alfa2);

    % ---- Transformación del segundo segmento ----
    Mc3d = [cx2, cy2*sin(alfa2), cz2, 1];
    theta_x = 0; theta_y = B1; theta_z = -alfa; theta_z2 = alfa;
    tx = cx3d; ty = cy3d; tz = cz3d;
    Rx = [1 0 0 0; 0 cos(theta_x) -sin(theta_x) 0; 0 sin(theta_x) cos(theta_x) 0; 0 0 0 1];
    Ry = [cos(theta_y) 0 sin(theta_y) 0; 0 1 0 0; -sin(theta_y) 0 cos(theta_y) 0; 0 0 0 1];
    Rz = [cos(theta_z) -sin(theta_z) 0 0; sin(theta_z) cos(theta_z) 0 0; 0 0 1 0; 0 0 0 1];
    Rz2 = [cos(theta_z2) -sin(theta_z2) 0 0; sin(theta_z2) cos(theta_z2) 0 0; 0 0 1 0; 0 0 0 1];
    T = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
    H = T*Rz2*Ry*Rz;
    Pf = (H*Mc3d')';
    Xf = Pf(1); Yf = Pf(2); Zf = Pf(3);
end

function J = jacobian_numeric(L)
    % Jacobiano numérico por diferencias finitas
    eps = 1e-3;
    [x0,y0,z0] = forward_end(L);
    f0 = [x0;y0;z0];
    J = zeros(3,length(L));
    for i = 1:length(L)
        Lp = L; Lp(i) = Lp(i) + eps;
        [xp,yp,zp] = forward_end(Lp);
        fp = [xp;yp;zp];
        J(:,i) = (fp - f0)/eps;
    end
end

function softArmForward(L)
    % Versión extendida del código 1 con graficado
    % Aquí se copia tu código 1, pero usando L como entrada
    L1=L(1); L2=L(2); L3=L(3);
    L4=L(4); L5=L(5); L6=L(6);
    
Lc = (L1 + L2 + L3)/3; % longitud central

% Eliminar indeterminación 1%
if L1==L2 && L2==L3
    Lc = Lc + 0.0001;
elseif (Lc-L1)==0
    L1 = L1 + 0.0001; Lc = (L1+L2+L3)/3;
elseif (Lc-L2)==0
    L2 = L2 + 0.0001; Lc = (L1+L2+L3)/3;
elseif (Lc-L3)==0
    L3 = L3 + 0.0001; Lc = (L1+L2+L3)/3;
end

% Eliminar indeterminación 2%
Lc2 = (L4 + L5 + L6)/3; % longitud central
if L4==L5 && L4==L6
    Lc2 = Lc2 + 0.0001;
elseif (Lc2-L4)==0
    L4 = L4 + 0.0001; Lc2 = (L4+L5+L6)/3;
elseif (Lc2-L5)==0
    L5 = L5 + 0.0001; Lc2 = (L4+L5+L6)/3;
elseif (Lc2-L6)==0
    L6 = L6 + 0.0001; Lc2 = (L4+L5+L6)/3;
end



%% Parámetros del brazo
r = 11.4; % radio de los tensores al centro del brazo

% Calculo de diferencias de longitudes
DL1 = Lc - L1; DL2 = Lc - L2; DL3 = Lc - L3;

% Dirección del brazo en el plano X,Y
Xdir = DL1 + DL2*cos(2*pi/3) + DL3*cos(4*pi/3);
Ydir = DL2*sin(2*pi/3) + DL3*sin(4*pi/3);
alfa = atan2(Ydir, Xdir);

% Proyección de los tensores en el plano de curvatura
rp1 = r*cos(alfa);
rp2 = r*cos(alfa + 4*pi/3);
rp3 = r*cos(alfa + 2*pi/3);

% Calculo del radio y ángulo de curvatura
R1 = (Lc*rp1)/(Lc-L1);
R2 = (Lc*rp2)/(Lc-L2);
R3 = (Lc*rp3)/(Lc-L3);

B1 = Lc/R1; B2 = Lc/R2; B3 = Lc/R3;

% Discretización del ángulo
theta1 = linspace(0, B1, 100);
theta2 = linspace(0, B2, 100);
theta3 = linspace(0, B3, 100);

% Coordenadas del arco interno
cz = R1*sin(theta1);
cy = R1*(1 - cos(theta1));
cx = zeros(size(theta1));

% Coordenadas de los arcos externos
z1 = (R1-rp1).*sin(theta1); y1 = ((R1-rp1).*(1 - cos(theta1))) + rp1; x1 = zeros(size(theta1));
z2 = (R2-rp2).*sin(theta2); y2 = ((R2-rp2).*(1 - cos(theta2))) + rp2; x2 = zeros(size(theta2));
z3 = (R3-rp3).*sin(theta3); y3 = ((R3-rp3).*(1 - cos(theta3))) + rp3; x3 = zeros(size(theta3));

% Proyección de los tensores en 3D
cz3d = cz;
cx3d = cy*cos(alfa);
cy3d = cy*sin(alfa);
z13d = z1;
x13d = y1*cos(alfa) - (rp1*cos(alfa)) + r;
y13d = y1*sin(alfa) - (rp1*sin(alfa));
z23d = z2;
x23d = y2*cos(alfa) - (rp2*cos(alfa)) + r*cos(2*pi/3);
y23d = y2*sin(alfa) - (rp2*sin(alfa)) + r*sin(2*pi/3);
z33d = z3;
x33d = y3*cos(alfa) - (rp3*cos(alfa)) + r*cos(4*pi/3);
y33d = y3*sin(alfa) - (rp3*sin(alfa)) + r*sin(4*pi/3);

%% Segundo segmento del brazo
r = 15; % radio de los tensores al centro del brazo
DL4 = Lc2 - L4; DL5 = Lc2 - L5; DL6 = Lc2 - L6;

d = 1/3; % ángulo de desfase entre tensores
X2dir = DL4*cos(d*pi) + DL5*cos(2/3*pi + d*pi) + DL6*cos(4/3*pi + d*pi);
Y2dir = DL4*sin(d*pi) + DL5*sin(2/3*pi + d*pi) + DL6*sin(4/3*pi + d*pi);
alfa2 = atan2(Y2dir, X2dir);

rp4 = r*cos(alfa2 - d*pi);
rp5 = r*cos(alfa2 + 4/3*pi - d*pi);
rp6 = r*cos(alfa2 + 2/3*pi - d*pi);

R4 = (Lc2*rp4)/(Lc2-L4);
R5 = (Lc2*rp5)/(Lc2-L5);
R6 = (Lc2*rp6)/(Lc2-L6);

B4 = Lc2/R4; B5 = Lc2/R5; B6 = Lc2/R6;

theta4 = linspace(0, B4, 100);
theta5 = linspace(0, B5, 100);
theta6 = linspace(0, B6, 100);

cz2 = R4*sin(theta4); cy2 = R4*(1-cos(theta4)); cx2 = zeros(size(theta4));
z4 = (R4-rp4).*sin(theta4); y4 = ((R4-rp4).*(1 - cos(theta4))) + rp4; x4 = zeros(size(theta4));
z5 = (R5-rp5).*sin(theta5); y5 = ((R5-rp5).*(1 - cos(theta5))) + rp5; x5 = zeros(size(theta5));
z6 = (R6-rp6).*sin(theta6); y6 = ((R6-rp6).*(1 - cos(theta6))) + rp6; x6 = zeros(size(theta6));

cz23d = cz2; cx23d = cy2*cos(alfa2); cy23d = cy2*sin(alfa2);
z43d = z4; x43d = y4*cos(alfa2) - rp4*cos(alfa2) + r*cos(d*pi); y43d = y4*sin(alfa2) - rp4*sin(alfa2) + r*sin(d*pi);
z53d = z5; x53d = y5*cos(alfa2) - rp5*cos(alfa2) + r*cos(2/3*pi + d*pi); y53d = y5*sin(alfa2) - rp5*sin(alfa2) + r*sin(2/3*pi + d*pi);
z63d = z6; x63d = y6*cos(alfa2) - rp6*cos(alfa2) + r*cos(4/3*pi + d*pi); y63d = y6*sin(alfa2) - rp6*sin(alfa2) + r*sin(4/3*pi + d*pi);

%% Transformación del segundo segmento
Mc3d = [cx23d', cy23d', cz23d', ones(length(cx23d),1)];
theta_x = 0; theta_y = B1; theta_z = -alfa; theta_z2 = alfa;
tx = cx3d(end); ty = cy3d(end); tz = cz3d(end);

Rx = [1 0 0 0; 0 cos(theta_x) -sin(theta_x) 0; 0 sin(theta_x) cos(theta_x) 0; 0 0 0 1];
Ry = [cos(theta_y) 0 sin(theta_y) 0; 0 1 0 0; -sin(theta_y) 0 cos(theta_y) 0; 0 0 0 1];
Rz = [cos(theta_z) -sin(theta_z) 0 0; sin(theta_z) cos(theta_z) 0 0; 0 0 1 0; 0 0 0 1];
Rz2 = [cos(theta_z2) -sin(theta_z2) 0 0; sin(theta_z2) cos(theta_z2) 0 0; 0 0 1 0; 0 0 0 1];

T = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];

H = T*Rz2*Ry*Rz;

Mc3d_T = (H*Mc3d')';
M43d = [x43d', y43d', z43d', ones(length(x43d),1)]; M53d = [x53d', y53d', z53d', ones(length(x53d),1)]; M63d = [x63d', y63d', z63d', ones(length(x63d),1)];
M43d_T = (H*M43d')'; M53d_T = (H*M53d')'; M63d_T = (H*M63d')';

%% Graficas de los brazos con círculos
n_circulos = 4;
idx = round(linspace(1, length(cx3d), n_circulos));
radio_circulo = r;

figure(4);
hold off;
plot3(cx3d, cy3d, cz3d, 'k', 'LineWidth', 2); % trayectoria promedio
hold on

for i = 1:n_circulos
    P = [cx3d(idx(i)), cy3d(idx(i)), cz3d(idx(i))];
    if i < n_circulos
        v = [cx3d(idx(i+1))-cx3d(idx(i)), cy3d(idx(i+1))-cy3d(idx(i)), cz3d(idx(i+1))-cz3d(idx(i))];
    else
        v = [cx3d(idx(i))-cx3d(idx(i-1)), cy3d(idx(i))-cy3d(idx(i-1)), cz3d(idx(i))-cz3d(idx(i-1))];
    end
    v = v/norm(v);
    if all(abs(v-[1 0 0]) < 1e-3)
        a = [0 1 0];
    else
        a = [1 0 0];
    end
    u = cross(v,a); u = u/norm(u);
    w = cross(v,u); w = w/norm(w);
    theta_circ = linspace(0, 2*pi, 100);
    circ = P' + radio_circulo*(u'*cos(theta_circ)+w'*sin(theta_circ));
    plot3(circ(1,:), circ(2,:), circ(3,:), 'b','LineWidth', 2.5);
    plot3(x13d, y13d, z13d, 'c', 'LineWidth', 1.5);
    plot3(x23d, y23d, z23d, 'g', 'LineWidth', 1.5);
    plot3(x33d, y33d, z33d, 'm', 'LineWidth', 1.5);
end

% Guardar primer segmento (coordenadas en el sistema mundial)
seg1X = cx3d; seg1Y = cy3d; seg1Z = cz3d;


% Brazo transformado
cx3d = Mc3d_T(:,1); cy3d = Mc3d_T(:,2); cz3d = Mc3d_T(:,3);
idx = round(linspace(1, length(Mc3d_T(:,1)), n_circulos));

plot3(Mc3d_T(:,1), Mc3d_T(:,2), Mc3d_T(:,3), 'b', 'LineWidth', 2);
hold on
for i = 1:n_circulos
    P = [cx3d(idx(i)), cy3d(idx(i)), cz3d(idx(i))];
    if i < n_circulos
        v = [cx3d(idx(i+1))-cx3d(idx(i)), cy3d(idx(i+1))-cy3d(idx(i)), cz3d(idx(i+1))-cz3d(idx(i))];
    else
        v = [cx3d(idx(i))-cx3d(idx(i-1)), cy3d(idx(i))-cy3d(idx(i-1)), cz3d(idx(i))-cz3d(idx(i-1))];
    end
    v = v/norm(v);
    if all(abs(v-[1 0 0]) < 1e-3)
        a = [0 1 0];
    else
        a = [1 0 0];
    end
    u = cross(v,a); u = u/norm(u);
    w = cross(v,u); w = w/norm(w);
    theta_circ = linspace(0,2*pi,100);
    circ = P' + radio_circulo*(u'*cos(theta_circ)+w'*sin(theta_circ));
    plot3(circ(1,:), circ(2,:), circ(3,:), 'r','LineWidth', 2.5);
    plot3(Mc3d_T(:,1), Mc3d_T(:,2), Mc3d_T(:,3), 'k', 'LineWidth', 1.5);
    plot3(M43d_T(:,1), M43d_T(:,2), M43d_T(:,3), 'c', 'LineWidth', 1.5);
    plot3(M53d_T(:,1), M53d_T(:,2), M53d_T(:,3), 'g', 'LineWidth', 1.5);
    plot3(M63d_T(:,1), M63d_T(:,2), M63d_T(:,3), 'm', 'LineWidth', 1.5);
end

Xf = Mc3d_T(end,1); Yf = Mc3d_T(end,2); Zf = Mc3d_T(end,3);
fprintf('Posición final del efector:\nX = %.2f mm\nY = %.2f mm\nZ = %.2f mm\n', Xf,Yf,Zf);


xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('soft arm');
axis equal;
grid on;
ylim([-100 100])
xlim([-100 100])
zlim([0 140])
view(3);

% ==== ROTAR LA VISTA PARA SIMULAR QUE CUELGA DEL TECHO ====
view(45, -20);      % mirar desde abajo
camup([0 0 -1]);    % invertir orientación vertical
set(gca, 'YDir', 'reverse');   % invierte el eje Y
set(gca,'ZDir','reverse');  % solo visual, no afecta dato

%% Proyecciones en las paredes (XY, XZ, YZ) para ambos segmentos pegadas a los ejes

% Segmento 1 (guardado previamente)
% seg1X, seg1Y, seg1Z

% Segmento 2 (efector transformado)
seg2X = Mc3d_T(:,1);
seg2Y = Mc3d_T(:,2);
seg2Z = Mc3d_T(:,3);

% Obtener límites actuales de la figura (asegúrate que estén definidos ya)
xl = xlim; yl = ylim; zl = zlim;

% --- Proyecciones del segmento 1 ---
% XY (Z = pared inferior zl(1))
plot3(seg1X, seg1Y, zl(1)*ones(size(seg1Z)), '--', 'LineWidth', 1.0);

% XZ (Y = pared trasera yl(1))
plot3(seg1X, yl(1)*ones(size(seg1X)), seg1Z, '--', 'LineWidth', 1.0);

% YZ (X = pared lateral xl(1))
plot3(xl(1)*ones(size(seg1Y)), seg1Y, seg1Z, '--', 'LineWidth', 1.0);

% --- Proyecciones del segmento 2 ---
% XY (Z = pared inferior zl(1))
plot3(seg2X, seg2Y, zl(1)*ones(size(seg2Z)), '--', 'LineWidth', 1.0);

% XZ (Y = pared trasera yl(1))
plot3(seg2X, yl(1)*ones(size(seg2X)), seg2Z, '--', 'LineWidth', 1.0);

% YZ (X = pared lateral xl(1))
plot3(xl(1)*ones(size(seg2Y)), seg2Y, seg2Z, '--', 'LineWidth', 1.0);


end



