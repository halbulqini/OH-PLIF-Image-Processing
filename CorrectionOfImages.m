% Case 1
clear all;
clc;
 filename='Avgbpm1cut.xlsx'; %reading beam profile file 
 A=xlsread(filename);
A=double(A);
n=max(max(A));
A=A./n;
% AA=A(100:600,:);
data = dir('*.tif');
z=0;
 for t=1:length(data)
    filename = data(t).name; 
    B=imread(filename);
    B=double(B);
%     BB=B(100:600,:);
%     c=BB./AA;
    c=B./A;
%     c(700:768,:)=0;               %solution for discarding all this bar
   for i=700:768               % solution for degradation
            for j=40:175
                if c(i,j)>1000
                c(i,j)=1000;
                end
            end
   end
% % y=imagesc(c,[5 4000]);
% axis image;colormap(jet);colorbar;
% title(['Image No.',j]);
% saveas(gcf,[j,'_.jpg']);
%  imshow(c,[])
%  colormap jet;
z=z+c;
j=int2str(t); 
mx=max(max(c));
% mn=min(min(c));
 y=imagesc(c,[100 mx]);
 axis image;colormap(jet);colorbar;
 title(['Image No.',j]);
 saveas(gcf,[j,'_.jpg']);
%  imwrite(c,jet,[j,'.jpg']);
 end
 zz=z./t; %averaging
 filename='cycavg.xlsx';
 xlswrite(filename,zz);
 %imwrite(zz,'cycavg.tif');