function TPMS_write_AbaqusInput_v2(folder_abbas,name_struct,G_comb,voxX,voxY,voxZ,W_U,H_U,L_U)

%folder_abbas = 'C:\Scratch\Functional_Grad_TPMS\abbas';

G_aux = G_comb;
aux = strrep( num2str(G_aux(:)','%1.2f,'), ' ','');
aux( length(aux): ceil(length(aux)/60)*60) = ' ';

fid = fopen([folder_abbas,'\',name_struct,'.abba'],'w');
nn = 1 ; 
for ll = 1 : length(aux)/60
    fprintf(fid, '%s\n', aux(nn:nn+59));
    nn = nn + 60;
end
fclose(fid);

fid = fopen([folder_abbas,'\',name_struct,'.hulshult'],'w');
fprintf(fid, '%s\n', num2str(voxY));
fprintf(fid, '%s\n', num2str(voxX));
fprintf(fid, '%s\n', num2str(voxZ));

fprintf(fid, '%s\n', num2str(size(G_aux,1)));
fprintf(fid, '%s\n', num2str(size(G_aux,2)));
fprintf(fid, '%s\n', num2str(size(G_aux,3)));

fprintf(fid, '%s\n', num2str(H_U/voxY));
fprintf(fid, '%s\n', num2str(W_U/voxX));
fprintf(fid, '%s\n', num2str(L_U/voxZ));
fclose(fid);

