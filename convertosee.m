% Convert from tif to jpg
clear all;
clc;
data = dir('*.tif');
 for i=1:20%length(data)   % Reading in different formats
 filename = data(i).name;
 h=imread(filename);
 D= uint8(h);
  j = int2str (i);
  imshow(D,[]);
 colormap jet;
 imwrite(D,jet,[j,'.jpg']);
 end
 