function [G_d, rho_e, G_S, G_H] = create_design_layers(G_S, G_H, W_G, D_G, fun_shape)

se = strel('sphere',1); %Object to mark the interface layers

G_Hi = imdilate(G_H,se); %'Dilate' or make the hard material one voxel bigger

clock = 0; % This clock means 0 is going through direction of soft, 1 is towards hard
G_i = G_Hi & G_S; % Obtain the interface by intersecting G_Hi with G_S
%figure; imshow3D( G_i )

% Now we obtain each of the different layers of the interface:
for ii = 2:2*W_G
    clock = 1 - clock; % Invert the direction!!!
    G_a = (sum(G_i(:,:,:,1:ii-1),4)>0); %Start by the previous layer of the interface
    G_as = (imdilate(G_a,se)&~G_a)&(G_S); %Dilate that layer in soft direction
    G_ah = (imdilate(G_a,se)&~G_a)&(G_H); %Dilate that layer in hard direction
    G_i(:,:,:,ii) = G_as*(1-clock) + G_ah*clock; % Store the alternating interface layer
end

if D_G == 0 % Store all the soft direction layers, that's the odds
    G_i = G_i(:,:,:,1:2:end);
    G_i = G_i(:,:,:,1:W_G);
elseif D_G == 1 % Store all the hard direction layers, that's the odds
    G_i = flip(G_i(:,:,:,2:2:end),4);
    G_i = G_i(:,:,:,1:W_G);
elseif D_G == 2 % We need to rearrange to begin from the hard side
    G_i = cat(4, flip(G_i(:,:,:,2:2:end),4), G_i(:,:,:,1:2:end));
    G_i = G_i(:,:,:,round(size(G_i,4)/4): round(size(G_i,4)/4) + W_G - 1);
end


G_H = G_H&~(sum(G_i,4)>0); % Remove the interface layers from the hard material
G_S = G_S&~(sum(G_i,4)>0); % Remove the interface layers from the soft material
G_i = cat(4,cat(4,G_H,G_i),G_S); % Incorporate the interface layers to hard and soft


rho_e = TPMS_generate_rho_func_v3(2*W_G, W_G+2, 1, 1, 1, 1,fun_shape);
rho_e = rho_e(2:W_G+3)';
rho_e(1) = 1;
rho_e(rho_e<0.001) = 0.001;

G_d = sum(  G_i.*permute(rho_e,[1,4,3,2]) ,4); 