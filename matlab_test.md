```matlab
I = imread('rice.png');
BG = imopen(I,strel('disk',15));
I2 = imsubtract(I,BG);
level = graythresh(I);
bw = im2bw(I,level); 
%%subplot(1,4,1);
%%imshow(bw);
level = graythresh(I2);
bw2 = im2bw(I2,level);    
%%subplot(1,4,2);
%%imshow(bw2);
[labeled,numObjects] = bwlabel(bw2,8);
Count = zeros(1,numObjects);
cnt = 1;
for i = 1:size(labeled,1)
    for j = 1:size(labeled,2)
        if labeled(i,j) ~= 0
            Count(labeled(i,j)) = Count(labeled(i,j)) + 1;
        end
    end
end
RGB_label = label2rgb(labeled); 
%%subplot(1,4,3);
%%imshow(RGB_label);
rgb = zeros(size(bw2,1),size(bw2,2),3);
idx = find(bw2);
rgb(:,:,1) = bw2; 
%%subplot(1,4,4);
imshow(rgb);
