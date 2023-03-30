function I_rgb = grs2rgb3D(G_des, botcol, topcol)
clear II_v II II_v2
% I_rgb = zeros([size(G_des),3]);
for ii = 1:size(G_des,3)
    aux = permute(G_des(:,:,ii),[1 2 3]);
    %auxM = ceil(aux);
    auxA = repelem(aux,1,1,1);
    auxB = repelem(1-aux,1,1,1);
    red = double(auxA).*topcol(1) + double(auxB).*botcol(1);
    green = double(auxA).*topcol(2) + double(auxB).*botcol(2);
    blue = double(auxA).*topcol(3) + double(auxB).*botcol(3);
    rgb_im = cat(3,cat( 3, red.*ceil(auxA) , green.*ceil(auxA)),blue.*ceil(auxA));
    
%     if size(rgb_im,3) == 3
%         
%     else
    I_rgb(:,:,:,ii)=rgb_im;%(cat(3,aux./aux-aux,0*aux,aux));
%     end
end
I_rgb = (permute(I_rgb,[1 2 4 3]));