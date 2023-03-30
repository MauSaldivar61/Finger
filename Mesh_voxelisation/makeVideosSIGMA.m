function makeVideosSIGMA(SIGMA,path,FR)

% % y = circshift( circshift(SIGMA,ceil(size(SIGMA,2)/4),2)   ,ceil(size(SIGMA,1)/4),1);
% % SIGMA = permute(repmat( cat(3,SIGMA,y),repss(1)+1,repss(2),2 ),[2,1,3] );
% figure; imshow3D(SIGMA)
% pause
% close figure
images    = cell(size(SIGMA,3),1);
for jj = 1: size(SIGMA,3)
    pic_aux = cell(1,1);
    %try
    %[img1,map1]  =   imread([foldername,  models{i,kkk,j,o} ,'-',num2str(jj),'.png'],'png');
    %x1rgb = ind2rgb(img1, map1);
    pic_aux{1} = uint8(squeeze(255*SIGMA(:,:,jj,:)));
    %catch
    %pic_aux{1} =imread([foldername, files{i}(1:end-7),variables{kkk},num2str(jj-1),'.png'],'png');
    %pic_aux{1} =imread([foldername,  models{i,kkk,j,o} ,'-',num2str(jj),'.png'],'png');
    %end
    images{jj} = pic_aux{1};%(178:565,1:1528,:);
end
% create the video writer with 1 fps
writerObj = VideoWriter(path);
writerObj.FrameRate = FR;
writerObj.Quality = 100;
open(writerObj);

cmap = gray(256);

% write the frames to the video
for frameI = 1:jj
    if size(images{jj},3)>1
    frame = im2frame(images{frameI});
    else
    frame = im2frame(images{frameI},cmap);
    end
    writeVideo(writerObj, frame);
end
close(writerObj);
