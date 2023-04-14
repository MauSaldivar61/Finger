function make_print_files_nofluid(folder_prints, name_struct, G_A_bits, G_B_bits, G_C_bits,...
                            Mat1,Mat2,Mat3,Mat4,Mat5,Mat6)

%G_H_bits; Bits of Hard Material
%G_S_bits; Bits of Soft Material
%G_F_bits; Bits of Fluid Material

% Time to save the design
if ~exist([folder_prints,'\',name_struct], 'dir')
    mkdir([folder_prints,'\',name_struct]);
else
    rmdir([folder_prints,'\',name_struct],'s');
end

if ~exist([folder_prints,'\',name_struct,'\Mat1'], 'dir')
    mkdir([folder_prints,'\',name_struct,'\Mat1']);
    mkdir([folder_prints,'\',name_struct,'\Mat2']);
    mkdir([folder_prints,'\',name_struct,'\Mat3']);
    mkdir([folder_prints,'\',name_struct,'\Mat4']);
    mkdir([folder_prints,'\',name_struct,'\Mat5']);
    mkdir([folder_prints,'\',name_struct,'\Mat6']);
end

for lay = 1:size(G_B_bits,3)
    imwrite( logical(G_A_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat1\','II_',num2str(lay,'%04.f'),'.bmp']);
    imwrite( logical(0*G_A_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat2\','II_',num2str(lay,'%04.f'),'.bmp']);
    imwrite( logical(G_C_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat3\','II_',num2str(lay,'%04.f'),'.bmp']);
    imwrite( logical(0*G_A_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat4\','II_',num2str(lay,'%04.f'),'.bmp']);
    imwrite( logical(0*G_A_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat5\','II_',num2str(lay,'%04.f'),'.bmp']);
    imwrite( logical(G_B_bits(:,:,lay)),[folder_prints,'\',name_struct,'\Mat6\','II_',num2str(lay,'%04.f'),'.bmp']);
end

fileID = fopen([folder_prints,'\',name_struct,'.txt'],'w');
fprintf(fileID,'[Build]\n');
fprintf(fileID,'Format version = 1;\n');
fprintf(fileID,'Build Mode = 3;\n');
fprintf(fileID,'Layer thickness = 0.027;\n');
fprintf(fileID,'Number of slices= %d;\n',lay);
fprintf(fileID,'\n');
fprintf(fileID,'[Resin Type]\n');
fprintf(fileID,['Color = ',Mat1,';\n']);%Mat1 = 'M.Cleanser'
fprintf(fileID,['Resin2 = ',Mat2,';\n']);%Mat2 = 'VeroClear'
fprintf(fileID,['Resin3 = ',Mat3,';\n']);%Mat3 = 'Agilus30Clr'
fprintf(fileID,['Resin4 = ',Mat4,';\n']);%Mat4 = 'VeroYL-V'
fprintf(fileID,['Resin5 = ',Mat5,';\n']);%Mat5 = 'VeroMGT-V'
fprintf(fileID,['Resin6 = ',Mat6,';\n']);%Mat6 = 'VeroCY-V'
fprintf(fileID,'\n');
fprintf(fileID,'[Materials]\n');

% Note: FolderPath has to be changed manually in the 3D printer PC!!!!!!!
fprintf(fileID,['Material1 = C:\\FolderPath\\',name_struct,'\\Mat1\\II_xxxx.bmp\n']);
fprintf(fileID,['Material2 = C:\\FolderPath\\',name_struct,'\\Mat2\\II_xxxx.bmp\n']);
fprintf(fileID,['Material3 = C:\\FolderPath\\',name_struct,'\\Mat3\\II_xxxx.bmp\n']);
fprintf(fileID,['Material4 = C:\\FolderPath\\',name_struct,'\\Mat4\\II_xxxx.bmp\n']);
fprintf(fileID,['Material5 = C:\\FolderPath\\',name_struct,'\\Mat5\\II_xxxx.bmp\n']);
fprintf(fileID,['Material6 = C:\\FolderPath\\',name_struct,'\\Mat6\\II_xxxx.bmp\n']);
fclose(fileID);
