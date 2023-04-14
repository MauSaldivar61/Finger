function [G_S_bits, G_H_bits] = make_bitmap_TPMS(G_greyscale, rho_e, G_S, G_H, U_xy, U_z, U_tpms, funk, sheet_q)

% Rule 1: The size of the greyscale to bitmap unit cell IS DIFFERENT than the TPMS unit cell size!!!!
% Rule 2: The code is run by getting the mesh of TPMS of the size of G_greyscale and then creating the TPMS equations
% The code follows by going through masks of every rho value and adjusting the value of hard material 

G_H_bits = repelem(G_greyscale,U_xy,U_xy,U_z);

% The Size of each unit-cell's grid;
vW = U_tpms; %Size of the grids
vH = U_tpms;
vL = U_tpms;
% The amount of unit cells:
nW = ceil(size(G_H_bits,2) / U_tpms);  %Number of units in X direction
nH = ceil(size(G_H_bits,1) / U_tpms);  %Number of units perpendicular to the direction
nL = ceil(size(G_H_bits,3) / U_tpms);  %Number of units perpendicular to the direction

% Define the structures that create our mesh:
%The linear spaces for the coordinates
W_g = nW*vW;
H_g = nH*vH;
L_g = nL*vL;
lin_W  = linspace(0,W_g,nW*vW);
lin_H  = linspace(0,H_g,nH*vH);
lin_L  = linspace(0,L_g,nL*vL);
%A mesh grid to translate the coordinates
[X,Y,Z] = meshgrid(lin_W,lin_H,lin_L);
%Calculates the total revolutions patterns the linspace and the units
kx  = 2*pi*nW/W_g;
ky  = 2*pi*nH/H_g;
kz  = 2*pi*nL/L_g;
%A mesh grid that is written in terms of radians
x = kx.*X;
y = ky.*Y;
z = kz.*Z;

% Define the Geometry:
if isequal(funk,'GYROID') || isequal(funk,'Gyroid') || isequal(funk,'gyroid')
    G = cos(x).*sin(y) + cos(y).*sin(z) + cos(z).*sin(x);
elseif isequal(funk,'DIAMOND') || isequal(funk,'Diamond') || isequal(funk,'diamond')
    G = sin(x).*sin(y).*sin(z) + sin(x).*cos(y).*cos(z) + cos(x).*sin(y).*cos(z) + cos(x).*cos(y).*cos(z);
elseif isequal(funk,'IWP') || isequal(funk,'Iwp') || isequal(funk,'iwp')
    G = cos(x).*cos(y) + cos(y).*cos(z) + cos(z).*cos(x);
elseif isequal(funk,'BCC') || isequal(funk,'Bcc') || isequal(funk,'bcc')
    G = (cos(2*x) + cos(2*y) + cos(2*z)) - 2*(cos(x).*cos(y) + cos(y).*cos(z) + cos(z).*cos(x));
elseif isequal(funk,'OCTO') || isequal(funk,'Octo') || isequal(funk,'octo')
    G = 4*(cos(x).*cos(y) + cos(y).*cos(z) + cos(z).*cos(x)) - 3*(cos(x) + cos(y) + cos(z));
elseif isequal(funk,'NEOVIUS') || isequal(funk,'Neovius') || isequal(funk,'neovius')
    G = 4*cos(x).*cos(y).*cos(z) + 3*(cos(x) + cos(y) + cos(z));
elseif isequal(funk,'SCHWARZ') || isequal(funk,'Schwarz') || isequal(funk,'schwarz')
    G = sin(x) + cos(y) + cos(z);
else
    error('I don''t know that geometry');
end

G_v = (G-min(min(min(G)))) ./ max(max(max(  G-min(min(min(G))) )))   ;
G_v = G_v(1:size(G_H_bits,1),1:size(G_H_bits,2),1:size(G_H_bits,3));

G_comb = zeros(size(G_v));

for lay = 1 : length(rho_e)
    rho_goal = rho_e(lay);
    G_aux = G_H_bits;
    G_aux(G_H_bits ~= rho_goal) = 0;
    G_aux(G_H_bits == rho_goal) = 1;
    
    G_v_aux = G_v.*G_aux;
    
    rho_test = rho_goal;
    
    %PLAY AROUND WITH THE NEXT 3
    rho_tol = 0.0005;
    rho_k = 0.05;
    end_count = 500;
    
    e_r = 1;
    countt = 0;
    clear G_comb_aux
    while abs(e_r) >  rho_tol && countt < end_count
        countt = countt + 1 ;
        if sheet_q == 1
            G_v1 = imbinarize(G_v_aux,0.5 + (rho_test)/2) .*G_aux;
            G_v2 = imbinarize(G_v_aux,0.5 - (rho_test)/2) .*G_aux;
            G_comb_aux = (~G_v1&G_v2);
        else
            G_comb_aux = double(G_v_aux<=rho_test) .*G_aux;
        end
        e_r = sum(sum(sum(G_comb_aux))) / sum(sum(sum(G_aux))) - rho_goal;
        
        rho_test = rho_test - rho_k*e_r;
    end
    G_comb = G_comb + G_comb_aux;
    
end

G_H_bits = logical(round(G_comb));
G_S_bits = ((~G_H_bits));

G_H_bits = G_H_bits&(repelem(G_H+G_S,U_xy,U_xy,U_z));
G_S_bits = G_S_bits&(repelem(G_H+G_S,U_xy,U_xy,U_z));
