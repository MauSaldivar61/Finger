function [rho_e,X_rho_W_T] = TPMS_generate_rho_func_v3(W_T, W_G, W_h, W_s, vox_xyz, b_t,type_lab)
% function [rho_e_ste,X_rho_W] = TPMS_generate_rho_func_v2(W_G, vox_xyz, b_Ym, n_t, saw, funk, rho_t, p_t, a_t)

% W_G is the total size of the gradient in mm
% W_U is the total size of each unit cell in mm
% vox_xyz is the voxel size in the gradient direction
% b_Ym is the b parameter in E = a rho^ b
% n_t is the amount of steps existing if a step function is desired
% saw is the amount of -+rho exists in every step
% funk is a string that expects either 'Sigmoid' or 'Linear'
% rho_t is the 'transition' density of a sigmoid (ONLY SIGMOID)
% p_t is the percentage of the length of the gradient at which rho_t occurs (ONLY SIGMOID)

% Defining the vector at which the gradient will be created

W_T_v = ceil(W_T/vox_xyz);
W_G_v = ceil(W_G/vox_xyz);
W_h_v = ceil(W_h/vox_xyz);
W_s_v = ceil(W_s/vox_xyz);


[W_T_v, W_G_v, W_h_v, W_s_v ];
[W_T_v, -2*W_G_v, -2*W_h_v, -W_s_v  ,2*W_G_v + 2*W_h_v + W_s_v];

X_rho_W_T  = linspace(0,W_T,W_T_v);

a_t = 3.333;

if isequal(type_lab,'Lin')
    rho_W_G  = linspace(0,1,W_G_v).^(1./b_t);
elseif isequal(type_lab,'Par')
    rho_W_G  = linspace(0,1,W_G_v).^(1);
elseif isequal(type_lab,'Sig')
    X_rho_W = linspace(0,W_G, W_G_v);
    d_e = 4*a_t / W_G;
    z_e = 1+exp(d_e.*(X_rho_W-W_G./2));
    rho_W_G = ( 1./z_e ).^(1./b_t); % The sigmoid function
    rho_W_G = flip(rho_W_G);
    rho_W_G = rho_W_G - min(rho_W_G)/2;
    [length(X_rho_W) length(rho_W_G) W_G_v];
    
elseif isequal(type_lab,'Logit')
    
% % %         X_rho_Wp = linspace(0,W_G, W_G_v);
% % %         d_e = 4*a_t / W_G;
% % %         z_e = 1+exp(d_e.*(X_rho_Wp-W_G./2));
% % %         rho_W_Gp = ( 1./z_e ).^(1./b_t); % The sigmoid function
% % %         rho_W_Gp = flip(rho_W_Gp);
% % %         rho_W_Gp = rho_W_Gp - min(rho_W_Gp)/2;
% % %         [length(X_rho_Wp) length(rho_W_Gp) W_G_v];
% % %     
% % %         aux = rho_W_Gp.*W_G;
% % %         rho_W_Gp = X_rho_Wp./W_G;
% % %         X_rho_Wp = aux;
% % %     
% % %         X_rho_W = linspace(0,W_G, W_G_v);
% % %         rho_W_G = interp1( X_rho_Wp, rho_W_Gp, X_rho_W);
% % %         rho_W_G(1) = 0.00;
% % %         rho_W_G(end) = 1 ;
    
    
    X_rho_W = linspace(0,W_G, W_G_v);
    d_e = 4*a_t / W_G;
    z_e = 1+exp(d_e.*(X_rho_W-W_G./2));
    rho_W_G = ( 1./z_e ).^(1./b_t); % The sigmoid function
    rho_W_G = flip(rho_W_G);
    rho_W_G = rho_W_G - min(rho_W_G)/2;
    
    rho_W_G = linspace(0,1,W_G_v).^(1./b_t) - ( rho_W_G -  linspace(0,1,W_G_v).^(1./b_t) );
    rho_W_G(rho_W_G<0) = 0;
    rho_W_G(rho_W_G>1) = 1;
    
elseif isequal(type_lab,'Ste')
    
    n_divisor = 6;    %INPUT
    n_t = round(W_G/(vox_xyz)./n_divisor);
    
    rho_W_G  =( repelem(linspace(1/n_divisor,1,W_G_v/n_t),n_t)  -1/n_divisor/2).^(1./b_t);
   
elseif isequal(type_lab,'Cos')
    X_rho_W = linspace(W_G,0, W_G_v);
    rho_W_G = (1/2.*( 1 + cos(X_rho_W*pi/W_G) )).^(1/b_t);
    
end




rho_W_h  = linspace(1,1,W_h_v);
rho_W_s  = linspace(0,0,W_s_v);


rho_e = [ rho_W_h rho_W_G(end:-1:1) rho_W_s  rho_W_G  rho_W_h]';

%     figure; plot(X_rho_W,rho_W_G)
%     figure; plot(X_rho_W_T,rho_e)
    
[length(X_rho_W_T), length(rho_e)];






