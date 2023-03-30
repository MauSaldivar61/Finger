function [G_comb,X,Y,Z] = TPMS_generate_Voxels(W_G,H_G,L_G, W_U,H_U,L_U, voxX,voxY,voxZ , funk, rho_func,sheet_q)
% W_G,H_G,L_G   are the sizes of the total sample in mm
% W_U,H_U,L_U   are the dimensions of each unit cell in mm
% voxX,voxY,voxZ are the voxel dimensions
% funk is the name of the function/geometry that you wish to create
% rho_func is the profile of density that one wishes to add to the geometry
%          NOTE OF rho_func: (IF DIMENSIONS DON'T MATCH WITH W_G, IT WILL INTERPOLATE)


% The Size of each unit-cell's grid;
vW = ceil(W_U/voxX); %Size of the grids
vH = ceil(H_U/voxY);
vL = ceil(L_U/voxZ);
% The amount of unit cells:
nW = (W_G/W_U);  %Number of units in X direction
nH = (H_G/H_U);  %Number of units perpendicular to the direction
nL = (L_G/L_U);  %Number of units perpendicular to the direction

% Define the structures that create our mesh:
%The linear spaces for the coordinates
lin_W  = linspace(0,W_G,nW*vW);
lin_H  = linspace(0,H_G,ceil(H_G/voxY));
lin_L  = linspace(0,L_G,nL*vL);
%A mesh grid to translate the coordinates
[X,Y,Z] = meshgrid(lin_W,lin_H,lin_L);
%Calculates the total revolutions patterns the linspace and the units
kx  = 2*pi*nW/W_G;
ky  = 2*pi*nH/H_G;
kz  = 2*pi*nL/L_G;
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
elseif isequal(funk,'COLL') || isequal(funk,'Coll') || isequal(funk,'coll')
    
    
    Do   = [H_U/(voxY)];
    L = nW*vW;
    revs = L/Do;
    reps = 3;
    ri = sqrt( 1/3.*(Do.^2)./(reps*pi) );
    rc = Do/2 - ri;
    ri = sqrt( rho_func.*(Do.^2)./(reps*pi) );
    t = linspace(0,2*pi,L);
    clear v w
    v = rc .* sin(revs*t + 2*pi*((1:reps)')/reps)';
    v = round(v - min(v) +1 +Do/2-rc);
    w = rc .* cos(revs*t + 2*pi*((1:reps)')/reps)';
    w = round(w - min(w) +1 +Do/2-rc);
    v(v>Do) = Do;
    w(w>Do) = Do;
    [cI, rI] = meshgrid(1:Do, 1:Do);
    
    lambda_fin =  zeros(Do,Do,L);
    
    for lay = 1:L
        
        dens_lay = rho_func(lay);
        lay_aux = zeros(Do,Do);
        for j = 1:reps
            cP = (rI - v(lay,j)).^2 + (cI - w(lay,j)).^2 <= ri(lay).^2;
            lay_aux = lay_aux | cP;
        end
        aux_d = sum(sum(sum(lay_aux))) / numel(lay_aux);
        while aux_d<dens_lay
            lay_aux = imdilate(  lay_aux  ,strel('sphere',1));
            aux_d = sum(sum(sum(lay_aux))) / numel(lay_aux);
            
        end
        vox2go_n = round((aux_d - dens_lay)*numel(lay_aux)); %voxels that have to go!!
        count_lim = 25;
        count = 0;
        while vox2go_n > 0 && count < count_lim
            count = count+1;
            lay_aux_e = lay_aux & edge(lay_aux,'Canny',0.5);
            lam_e_ind = find(lay_aux_e(:));
            
            try
                kill_ind = randperm(numel(lam_e_ind), vox2go_n)';
            catch
                kill_ind = randperm(numel(lam_e_ind), numel(lam_e_ind))';
            end
            lay_aux(lam_e_ind(kill_ind)) = 0;
            aux_d = sum(sum(sum(lay_aux))) / numel(lay_aux);
            vox2go_n = round((aux_d - dens_lay)*numel(lay_aux));
            if dens_lay == 0
                lay_aux(:) = 0;
                vox2go_n = 0;
            end
        end
        lambda_fin(:,:,lay) = lay_aux;
    end
    % figure; imshow3D(lambda_fin)
    % figure; imshow3D(permute(lambda_fin,[2,3,1]))
    G_comb = repmat(  permute(lambda_fin,[1,3,2]), round(H_G/H_U), 1, round(L_G/L_U));
elseif isequal(funk,'PAR') || isequal(funk,'Par') || isequal(funk,'par')
    G_comb = PGRAD_generate_Voxels(W_G, H_G, L_G,...
                                   voxX,  voxY,  voxZ, ...
                                   rho_func);
    
else
    error('I don''t know that geometry');
end

% Change the range of -1.5 -> 1.5 to 0 -> 1
if ~(isequal(funk,'COLL') || isequal(funk,'Coll') || isequal(funk,'coll') || isequal(funk,'PAR') || isequal(funk,'Par') || isequal(funk,'par'))
    
    G_v = (G-min(min(min(G)))) ./ max(max(max(  G-min(min(min(G))) )))   ;
    max(G_v(:));
    % UNCOMMENT TO TEST THE GRAYSCALE:
    % figure
    % imshow3D(permute(G_v,[3,1,2]))
    %  figure
    %  imshow3D(G_v)
    
    
    % A SHEET GRADIENT ACROSS X THAT TRANSITIONS ACCORDING TO A GIVEN DENSITY FUNCTION
    % close all
    G_comb = zeros(size(G_v));
    nW*vW;
    
    
    for lay = 1 : nW*vW
        rho_goal = rho_func(lay);
        G_v_aux = G_v(:,lay,:);
        
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
                G_v1 = imbinarize(G_v_aux,0.5 + (rho_test)/2);
                G_v2 = imbinarize(G_v_aux,0.5 - (rho_test)/2);
                G_comb_aux = (~G_v1&G_v2);
            else
                %         G_v2 = imbinarize(G_v_aux,rho_test);
                G_comb_aux = double(G_v_aux<=rho_test);
                %         sum(sum(sum(G_comb_aux)))/numel(G_comb_aux)
            end
            e_r = sum(sum(sum(G_comb_aux))) / numel(G_comb_aux) - rho_goal;
            
            rho_test = rho_test - rho_k*e_r;
        end
        countt;
        rho_test;
        G_comb(:,lay,:) = reshape(G_comb_aux,[size(G_comb,1),1,size(G_comb,3)]);
        %     if lay < 10
        %     figure; imshow(squeeze(G_comb_aux));
        %     end
        
    end
end
% figure; imshow3D(squeeze(G_comb));