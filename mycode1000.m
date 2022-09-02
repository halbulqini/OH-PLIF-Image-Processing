clear; clc; close all;
pix_size =0.18797;

prompt='What is the index of 1st image? ';
i_start = input(prompt);
prompt='What is the index of final image? ';
i_end= input(prompt);
prompt='What is your case? ';
case_name= input(prompt,'s');

B_t=cell([i_end,10]); %it will B_te filled with every contour as it is.

B_1=cell([i_end,6]);
B_2=cell([i_end,6]);
B_3=cell([i_end,6]);

data = dir('*.jpg');
for t=i_start:i_end 
filename = data(t).name; 

%Reading image in the current folder
H=imread(filename);

%Transfering into gray scale and cropping the reaction zone
H_gray= rgb2gray(H);
trimm(:,:)=H_gray(57:515,146:665);
%zone(:,:)=trimm(180:end,120:400);
zn=double(trimm);

%Normalizing
m1=max(max(zn));
norm_im=zn./m1;
[Lf,Wf]=size(norm_im);
%figure(1), imshow(norm_im); title('row image');

%modifying contrast of row images (Hazem 09/08/2020)
se1=strel('disk',2);
dl0=imdilate(norm_im,se1);
%figure(2), imshow(dl0); title('modified row image');

% 1-Binning by resizing
%bn=imresize(norm_im,0.5,'bilinear'); (No need after 9 August's mod.)

% 2- Smoothing and Sharpening 
sth=imgaussfilt(dl0,2);
imsh = imsharpen(sth,'Amount',1.2);

% % 3-Gradient
% grd=imgradient(sth);
% mx1=max(max(grd));
% grd1=grd./mx1;
% [Lo,Wo]=size(grd1);
% %figure(3), imshow(grd1); title('gradient image');

% 4-Threshold (Otsu's Method)
level = graythresh(imsh);
BW = im2bw(imsh,level);
[Loo,Woo]=size(BW);
%figure(4), imshow(BW); title('After Threshold');

%Closing open contours
se = strel('disk',1);
BW1 = imclose(BW,se);
ths = BW1;

%Solving the problem of salt and pepper noise in the final image
se2=strel('disk',1);
dll=imdilate(ths,se2);
%BW2 = imfill( dll ,'holes');
z0=edge(dll);
%figure(5), imshow(z); title('Taking The Edge');

if (z0(end,1)==1) && (z0(end,2)==0) && (z0(end-1,1)==0) 
    z0(end,1)=0;
end
 for i=3:Loo
     for j=3:Woo
        if (z0(i-1,j-1)==1) && (z0 (i-2,j-2)==0) && (z0(i-1,j-2)==0) && (z0(i,j-2)==0) && (z0(i-2,j-1)==0) && (z0(i,j-1)==0) && (z0(i-2,j)==0) && (z0(i-1,j)==0) && (z0(i,j)==0)
            z0(i-1,j-1)=0;
        end
    end
end

%% Taking the edges of the reaction zone (Hazem's Addition)
z=double(z0);
[L,W]=size(z);

z_t=flipud(z); %to fix the origin point to be at the bottom center

%Divide the image into 3 zones to study each one of them
z_1=z_t(1:153,:);
z_2=z_t(154:306,:);
z_3=z_t(306:end,:);

d_pix_c=2; %The number of pixels to calculate curvature on. (make sure it is EVEN number)
%% Using bwboundaries
[boundaries,~,~,~] =bwboundaries(z_t,8);
for i=1:length(boundaries)
    if length(boundaries{i})<=9 %To avoid error in filtfilt
        continue
    end
  x_t{i}=boundaries{i}(:,2);
  y_t{i}=boundaries{i}(:,1);
  ds_t = sqrt((x_t{i}(2:end)-x_t{i}(1:end-1)).^2  +  (y_t{i}(2:end)-y_t{i}(1:end-1)).^2  ) ; 
  ds_t = [ 1 ds_t' ] ;
  
   x_t{i} = x_t{i}(ds_t>0) ;       y_t{i} = y_t{i}(ds_t>0) ;       ds_t = ds_t(2:end) ;        ds_t = ds_t(ds_t>0) ;
   s_t = cumsum( [0 ds_t] ) ;      s_t = s_t/s_t(end) ;

   DS_t{i}=sqrt((x_t{i}(2:end)-x_t{i}(1:end-1)).^2  +  (y_t{i}(2:end)-y_t{i}(1:end-1)).^2);
   
%    want_to_display=true;
%      if want_to_display
%      plot (x{i},y{i},'b','LineWidth',2)
%      hold on 
%      end
end


%Rotating 
xx_t=rot90(x_t,3);

%To fix the origin point to be at the bottom center
for c0=1:length(xx_t)
    xx1_t{c0,1}=xx_t{c0,1}-259; %(half of W)
end
yy_t=rot90(y_t,3);
DS2_t=rot90(DS_t,3);

%Excluding NaN cells
xx_t = xx_t(cellfun(@(x)~isempty(x), xx_t));
xx1_t = xx1_t(cellfun(@(x)~isempty(x), xx1_t));
yy_t = yy_t(cellfun(@(x)~isempty(x), yy_t));
DS2_t = DS2_t(cellfun(@(x)~isempty(x), DS2_t));

%Smoothing the Curve to avoid sharp edges
[b,a]=butter(3,0.1);

for c00=1:length(xx_t)
    filt_x_t_2{c00,1}=filtfilt(b,a,xx_t{c00,1});
end

for c1=1:length(xx1_t)
    filt_x_t{c1,1}=filtfilt(b,a,xx1_t{c1,1});
end
for c2=1:length(yy_t)
    filt_y_t{c2,1}=filtfilt(b,a,yy_t{c2,1});
end

%Multiplying by pixel size

for j0=1:length(filt_x_t_2)
    X2_t{j0,1}=pix_size*filt_x_t_2{j0,1};
end

for j5=1:length(filt_x_t)
    X_t{j5,1}=pix_size*filt_x_t{j5,1};
end

for j6=1:length(filt_y_t)
    Y_t{j6,1}=pix_size*filt_y_t{j6,1};
end

% for j7=1:length(DS2_t)
%     D_S_t{j7,1}=pix_size*DS2_t{j7,1};
% end

for j8=1:length(X_t)
DS_filt_t{j8,1}=sqrt((X_t{j8,1}(2:end)-X_t{j8,1}(1:end-1)).^2  +  (Y_t{j8,1}(2:end)-Y_t{j8,1}(1:end-1)).^2);
end

for j1=1:length(X2_t)
dx{j1,1}=X2_t{j1,1}(2:end)-X2_t{j1,1}(1:end-1);
end

for j2=1:length(Y_t)
dy{j2,1}=Y_t{j2,1}(2:end)-Y_t{j2,1}(1:end-1);
end

B_t{t,1}=X_t; B_t{t,2}=Y_t; B_t{t,3}=DS_filt_t; B_t{t,4}=xx1_t; B_t{t,5}=yy_t; B_t{t,6}=DS2_t; B_t{t,7}=xx_t; B_t{t,8}=X2_t; B_t{t,9}=dx; B_t{t,10}=dy;

X_t=[];
Y_t=[];
DS_filt_t=[];
xx1_t=[];
yy_t=[];
DS2_t=[];
X2_t=[];
dx=[];
dy=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Zone 1
[boundaries_1,~,~,~] =bwboundaries(z_1,8);
for i=1:length(boundaries_1)
    if length(boundaries_1{i})<=9 %To avoid error in filtfilt
        continue
    end
  x_1{i}=boundaries_1{i}(:,2);
  y_1{i}=boundaries_1{i}(:,1);
  ds_1 = sqrt((x_1{i}(2:end)-x_1{i}(1:end-1)).^2  +  (y_1{i}(2:end)-y_1{i}(1:end-1)).^2  ) ; 
  ds_1 = [ 1 ds_1' ] ;
  
   x_1{i} = x_1{i}(ds_1>0) ;       y_1{i} = y_1{i}(ds_1>0) ;       ds_1 = ds_1(2:end) ;        ds_1 = ds_1(ds_1>0) ;
   s_1 = cumsum( [0 ds_1] ) ;      s_1 = s_1/s_1(end) ;

   DS_1{i}=sqrt((x_1{i}(2:end)-x_1{i}(1:end-1)).^2  +  (y_1{i}(2:end)-y_1{i}(1:end-1)).^2);
   
%    want_to_display=true;
%      if want_to_display
%      plot (x{i},y{i},'b','LineWidth',2)
%      hold on 
%      end
end


%Rotating 
xx_1=rot90(x_1,3);

%To fix the origin point to be at the bottom center
for c0=1:length(xx_1)
    xx1_1{c0,1}=xx_1{c0,1}-259; %(half of W)
end
yy_1=rot90(y_1,3);
DS2_1=rot90(DS_1,3);

%Excluding NaN cells
xx1_1 = xx1_1(cellfun(@(x)~isempty(x), xx1_1));
yy_1 = yy_1(cellfun(@(x)~isempty(x), yy_1));
DS2_1 = DS2_1(cellfun(@(x)~isempty(x), DS2_1));

%Smoothing the Curve to avoid sharp edges
[b,a]=butter(3,0.1);

for c1=1:length(xx1_1)
    filt_x_1{c1,1}=filtfilt(b,a,xx1_1{c1,1});
end
for c2=1:length(yy_1)
    filt_y_1{c2,1}=filtfilt(b,a,yy_1{c2,1});
end

%Multiplying by pixel size
for j5=1:length(filt_x_1)
    X_1{j5,1}=pix_size*filt_x_1{j5,1};
end

for j6=1:length(filt_y_1)
    Y_1{j6,1}=pix_size*filt_y_1{j6,1};
end

% for j7=1:length(DS2_1)
%     D_S_1{j7,1}=pix_size*DS2_1{j7,1};
% end

for j8=1:length(X_1)
DS_filt_1{j8,1}=sqrt((X_1{j8,1}(2:end)-X_1{j8,1}(1:end-1)).^2  +  (Y_1{j8,1}(2:end)-Y_1{j8,1}(1:end-1)).^2);
end

B_1{t,1}=X_1; B_1{t,2}=Y_1; B_1{t,3}=DS_filt_1; B_1{t,4}=xx1_1; B_1{t,5}=yy_1; B_1{t,6}=DS2_1;

X_1=[];
Y_1=[];
DS_filt_1=[];
xx1_1=[];
yy_1=[];
DS2_1=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% zone 2

[boundaries_2,~,~,~] =bwboundaries(z_2,8);
for i=1:length(boundaries_2)
    if length(boundaries_2{i})<=9 %To avoid error in filtfilt
        continue
    end
  x_2{i}=boundaries_2{i}(:,2);
  y_2{i}=boundaries_2{i}(:,1);
  ds_2 = sqrt((x_2{i}(2:end)-x_2{i}(1:end-1)).^2  +  (y_2{i}(2:end)-y_2{i}(1:end-1)).^2  ) ; 
  ds_2 = [ 1 ds_2' ] ;
  
   x_2{i} = x_2{i}(ds_2>0) ;       y_2{i} = y_2{i}(ds_2>0) ;       ds_2 = ds_2(2:end) ;        ds_2 = ds_2(ds_2>0) ;
   s_2 = cumsum( [0 ds_2] ) ;      s_2 = s_2/s_2(end) ;

   DS_2{i}=sqrt((x_2{i}(2:end)-x_2{i}(1:end-1)).^2  +  (y_2{i}(2:end)-y_2{i}(1:end-1)).^2);
   
%    want_to_display=true;
%      if want_to_display
%      plot (x{i},y{i},'b','LineWidth',2)
%      hold on 
%      end
end


%Rotating 
xx_2=rot90(x_2,3);

%To fix the origin point to be at the bottom center
for c0=1:length(xx_2)
    xx1_2{c0,1}=xx_2{c0,1}-259; %(half of W)
end
yy_2=rot90(y_2,3);
DS2_2=rot90(DS_2,3);

%Excluding NaN cells
xx1_2 = xx1_2(cellfun(@(x)~isempty(x), xx1_2));
yy_2 = yy_2(cellfun(@(x)~isempty(x), yy_2));
DS2_2 = DS2_2(cellfun(@(x)~isempty(x), DS2_2));

%Smoothing the Curve to avoid sharp edges
[b,a]=butter(3,0.1);

for c1=1:length(xx1_2)
    filt_x_2{c1,1}=filtfilt(b,a,xx1_2{c1,1});
end
for c2=1:length(yy_2)
    filt_y_2{c2,1}=filtfilt(b,a,yy_2{c2,1});
end

%Multiplying by pixel size
for j5=1:length(filt_x_2)
    X_2{j5,1}=pix_size*filt_x_2{j5,1};
end

for j6=1:length(filt_y_2)
    Y_2{j6,1}=pix_size*filt_y_2{j6,1};
end

% for j7=1:length(DS2_2)
%     D_S_2{j7,1}=pix_size*DS2_2{j7,1};
% end

for j8=1:length(X_2)
DS_filt_2{j8,1}=sqrt((X_2{j8,1}(2:end)-X_2{j8,1}(1:end-1)).^2  +  (Y_2{j8,1}(2:end)-Y_2{j8,1}(1:end-1)).^2);
end

B_2{t,1}=X_2; B_2{t,2}=Y_2; B_2{t,3}=DS_filt_2; B_2{t,4}=xx1_2; B_2{t,5}=yy_2; B_2{t,6}=DS2_2;

X_2=[];
Y_2=[];
DS_filt_2=[];
xx1_2=[];
yy_2=[];
DS2_2=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zone 3
[boundaries_3,~,~,~] =bwboundaries(z_3,8);
for i=1:length(boundaries_3)
    if length(boundaries_3{i})<=9 %To avoid error in filtfilt
        continue
    end
  x_3{i}=boundaries_3{i}(:,2);
  y_3{i}=boundaries_3{i}(:,1);
  ds_3 = sqrt((x_3{i}(2:end)-x_3{i}(1:end-1)).^2  +  (y_3{i}(2:end)-y_3{i}(1:end-1)).^2  ) ; 
  ds_3 = [ 1 ds_3' ] ;
  
   x_3{i} = x_3{i}(ds_3>0) ;       y_3{i} = y_3{i}(ds_3>0) ;       ds_3 = ds_3(2:end) ;        ds_3 = ds_3(ds_3>0) ;
   s_3 = cumsum( [0 ds_3] ) ;      s_3 = s_3/s_3(end) ;

   DS_3{i}=sqrt((x_3{i}(2:end)-x_3{i}(1:end-1)).^2  +  (y_3{i}(2:end)-y_3{i}(1:end-1)).^2);
   
%    want_to_display=true;
%      if want_to_display
%      plot (x{i},y{i},'b','LineWidth',2)
%      hold on 
%      end
end


%Rotating 
xx_3=rot90(x_3,3);

%To fix the origin point to be at the bottom center
for c0=1:length(xx_3)
    xx1_3{c0,1}=xx_3{c0,1}-259; %(half of W)
end
yy_3=rot90(y_3,3);
DS2_3=rot90(DS_3,3);

%Excluding NaN cells
xx1_3 = xx1_3(cellfun(@(x)~isempty(x), xx1_3));
yy_3 = yy_3(cellfun(@(x)~isempty(x), yy_3));
DS2_3 = DS2_3(cellfun(@(x)~isempty(x), DS2_3));

%Smoothing the Curve to avoid sharp edges
[b,a]=butter(3,0.1);

for c1=1:length(xx1_3)
    filt_x_3{c1,1}=filtfilt(b,a,xx1_3{c1,1});
end
for c2=1:length(yy_3)
    filt_y_3{c2,1}=filtfilt(b,a,yy_3{c2,1});
end

%Multiplying by pixel size
for j5=1:length(filt_x_3)
    X_3{j5,1}=pix_size*filt_x_3{j5,1};
end

for j6=1:length(filt_y_3)
    Y_3{j6,1}=pix_size*filt_y_3{j6,1};
end

% for j7=1:length(DS2_3)
%     D_S_3{j7,1}=pix_size*DS2_3{j7,1};
% end

for j8=1:length(X_3)
DS_filt_3{j8,1}=sqrt((X_3{j8,1}(2:end)-X_3{j8,1}(1:end-1)).^2  +  (Y_3{j8,1}(2:end)-Y_3{j8,1}(1:end-1)).^2);
end

B_3{t,1}=X_3; B_3{t,2}=Y_3; B_3{t,3}=DS_filt_3; B_3{t,4}=xx1_3; B_3{t,5}=yy_3; B_3{t,6}=DS2_3;

X_3=[];
Y_3=[];
DS_filt_3=[];
xx1_3=[];
yy_3=[];
DS2_3=[];
end

%% 1- Flame Curvature Calculations

 curv_limits = [-5 5] ;
 curv_limits2= [-5 5];
 number_of_bins = 500 ;
 curv_axis = linspace(curv_limits(1), curv_limits(2), number_of_bins);
 curv_axis2 = linspace(curv_limits2(1), curv_limits2(2), number_of_bins);
 d_curve = curv_axis(2) - curv_axis(1) ;

 % For the whole image (_t)

for c=1:length(B_t)
for i = 1:length(B_t{c,1})
    for k=d_pix_c+1:length(B_t{c,1}{i,1})
  
    Kapa_n_t{i,1}(k-d_pix_c,1)=8*(B_t{c,1}{i,1}(k-d_pix_c,1).*(B_t{c,2}{i,1}(k,1) - B_t{c,2}{i,1}(k-0.5*d_pix_c,1))+B_t{c,1}{i,1}(k-0.5*d_pix_c,1).*(B_t{c,2}{i,1}(k-d_pix_c,1) - B_t{c,2}{i,1}(k,1))+B_t{c,1}{i,1}(k,1).*(B_t{c,2}{i,1}(k-0.5*d_pix_c,1) - B_t{c,2}{i,1}(k-d_pix_c,1)))/(((B_t{c,1}{i,1}(k,1)-B_t{c,1}{i,1}(k-d_pix_c,1)).^2+(B_t{c,2}{i,1}(k,1) - B_t{c,2}{i,1}(k-d_pix_c,1)).^2).^1.5);
   
    
    end  
    
end
Kapa1_n_t=cell2mat(Kapa_n_t);
Kapa1_n_t = Kapa1_n_t(Kapa1_n_t>=-5 & Kapa1_n_t<=5);
kapa_all_t{c,:}=Kapa1_n_t;
Kapa_n_t=[];
NDF_curv_t=[];
NDF_curv_t (:,c) = hist(Kapa1_n_t,curv_axis);
NDF_of_curvatures_t = NDF_curv_t;
NDF_of_curvatures_t (:,c) = sum(NDF_of_curvatures_t,2);
PDF_curvature_t(:,c) = NDF_of_curvatures_t / ( d_curve * sum(NDF_of_curvatures_t));
end

pdf_n_t=sum(PDF_curvature_t,2)./i_end;

filtpdf_t=filtfilt(b,a,pdf_n_t(2:end-1));
Axes_curv_t(:,1)=curv_axis(2:end-1); Axes_curv_t(:,2)=filtpdf_t;
plot(curv_axis2(2:end-1),filtpdf_t,'-k', 'LineWidth',2); hold on;

kapa_mat_t=cell2mat(kapa_all_t);
kapa_mat2_t = kapa_mat_t(kapa_mat_t>=-1 & kapa_mat_t<=1);
radii_t=1./kapa_mat2_t;

Axes_kr_t(1,1) = mean(kapa_mat2_t);
Axes_kr_t(2,1)= median(kapa_mat2_t);
%Axes_kr_t(3,1)= mean(radii_t);
Axes_kr_t(3,1)= median(radii_t);

save kapa_mat_t_1000.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For zone 1
for c=1:length(B_1)
for i = 1:length(B_1{c,1})
    for k=d_pix_c+1:length(B_1{c,1}{i,1})
  
    Kapa_n_1{i,1}(k-d_pix_c,1)=8*(B_1{c,1}{i,1}(k-d_pix_c,1).*(B_1{c,2}{i,1}(k,1) - B_1{c,2}{i,1}(k-0.5*d_pix_c,1))+B_1{c,1}{i,1}(k-0.5*d_pix_c,1).*(B_1{c,2}{i,1}(k-d_pix_c,1) - B_1{c,2}{i,1}(k,1))+B_1{c,1}{i,1}(k,1).*(B_1{c,2}{i,1}(k-0.5*d_pix_c,1) - B_1{c,2}{i,1}(k-d_pix_c,1)))/(((B_1{c,1}{i,1}(k,1)-B_1{c,1}{i,1}(k-d_pix_c,1)).^2+(B_1{c,2}{i,1}(k,1) - B_1{c,2}{i,1}(k-d_pix_c,1)).^2).^1.5);
   
    
    end  
    
end
Kapa1_n_1=cell2mat(Kapa_n_1);
Kapa1_n_1 = Kapa1_n_1(Kapa1_n_1>=-5 & Kapa1_n_1<=5);
kapa_all_1{c,:}=Kapa1_n_1;
Kapa_n_1=[];
NDF_curv_1=[];
NDF_curv_1 (:,c) = hist(Kapa1_n_1,curv_axis);
NDF_of_curvatures_1 = NDF_curv_1;
NDF_of_curvatures_1 (:,c) = sum(NDF_of_curvatures_1,2);
PDF_curvature_1(:,c) = NDF_of_curvatures_1 / ( d_curve * sum(NDF_of_curvatures_1));
end

pdf_n_1=sum(PDF_curvature_1,2)./i_end;

filtpdf_1=filtfilt(b,a,pdf_n_1(2:end-1));
Axes_curv_1(:,1)=curv_axis(2:end-1); Axes_curv_1(:,2)=filtpdf_1;
plot(curv_axis2(2:end-1),filtpdf_1,'-b', 'LineWidth',2); hold on;

kapa_mat_1=cell2mat(kapa_all_1);
kapa_mat2_1 = kapa_mat_1(kapa_mat_1>=-1 & kapa_mat_1<=1);
radii_1=1./kapa_mat2_1;

kr_1(1,1) = mean(kapa_mat2_1);
kr_1(2,1)= median(kapa_mat2_1);
%kr_1(3,1)= mean(radii_1);
kr_1(3,1)= median(radii_1);

save kapa_mat_1_1000.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For zone 2

for c=1:length(B_2)
for i = 1:length(B_2{c,1})
    for k=d_pix_c+1:length(B_2{c,1}{i,1})
  
    Kapa_n_2{i,1}(k-d_pix_c,1)=8*(B_2{c,1}{i,1}(k-d_pix_c,1).*(B_2{c,2}{i,1}(k,1) - B_2{c,2}{i,1}(k-0.5*d_pix_c,1))+B_2{c,1}{i,1}(k-0.5*d_pix_c,1).*(B_2{c,2}{i,1}(k-d_pix_c,1) - B_2{c,2}{i,1}(k,1))+B_2{c,1}{i,1}(k,1).*(B_2{c,2}{i,1}(k-0.5*d_pix_c,1) - B_2{c,2}{i,1}(k-d_pix_c,1)))/(((B_2{c,1}{i,1}(k,1)-B_2{c,1}{i,1}(k-d_pix_c,1)).^2+(B_2{c,2}{i,1}(k,1) - B_2{c,2}{i,1}(k-d_pix_c,1)).^2).^1.5);
   
    
    end  
    
end
Kapa1_n_2=cell2mat(Kapa_n_2);
Kapa1_n_2 = Kapa1_n_2(Kapa1_n_2>=-5 & Kapa1_n_2<=5);
kapa_all_2{c,:}=Kapa1_n_2;
Kapa_n_2=[];
NDF_curv_2=[];
NDF_curv_2 (:,c) = hist(Kapa1_n_2,curv_axis);
NDF_of_curvatures_2 = NDF_curv_2;
NDF_of_curvatures_2 (:,c) = sum(NDF_of_curvatures_2,2);
PDF_curvature_2(:,c) = NDF_of_curvatures_2 / ( d_curve * sum(NDF_of_curvatures_2));
end

pdf_n_2=sum(PDF_curvature_2,2)./i_end;

filtpdf_2=filtfilt(b,a,pdf_n_2(2:end-1));
Axes_curv_2(:,1)=curv_axis(2:end-1); Axes_curv_2(:,2)=filtpdf_2;
plot(curv_axis2(2:end-1),filtpdf_2,'-g', 'LineWidth',2); hold on;

kapa_mat_2=cell2mat(kapa_all_2);
kapa_mat2_2 = kapa_mat_2(kapa_mat_2>=-1 & kapa_mat_2<=1);
radii_2=1./kapa_mat2_2;

kr_2(1,1) = mean(kapa_mat2_2);
kr_2(2,1)= median(kapa_mat2_2);
%kr_2(3,1)= mean(radii_2);
kr_2(3,1)= median(radii_2);

save kapa_mat_2_1000.mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For zone 3

for c=1:length(B_3)
for i = 1:length(B_3{c,1})
    for k=d_pix_c+1:length(B_3{c,1}{i,1})
  
    Kapa_n_3{i,1}(k-d_pix_c,1)=8*(B_3{c,1}{i,1}(k-d_pix_c,1).*(B_3{c,2}{i,1}(k,1) - B_3{c,2}{i,1}(k-0.5*d_pix_c,1))+B_3{c,1}{i,1}(k-0.5*d_pix_c,1).*(B_3{c,2}{i,1}(k-d_pix_c,1) - B_3{c,2}{i,1}(k,1))+B_3{c,1}{i,1}(k,1).*(B_3{c,2}{i,1}(k-0.5*d_pix_c,1) - B_3{c,2}{i,1}(k-d_pix_c,1)))/(((B_3{c,1}{i,1}(k,1)-B_3{c,1}{i,1}(k-d_pix_c,1)).^2+(B_3{c,2}{i,1}(k,1) - B_3{c,2}{i,1}(k-d_pix_c,1)).^2).^1.5);
   
    
    end  
    
end
Kapa1_n_3=cell2mat(Kapa_n_3);
Kapa1_n_3 = Kapa1_n_3(Kapa1_n_3>=-5 & Kapa1_n_3<=5);
kapa_all_3{c,:}=Kapa1_n_3;
Kapa_n_3=[];
NDF_curv_3=[];
NDF_curv_3 (:,c) = hist(Kapa1_n_3,curv_axis);
NDF_of_curvatures_3 = NDF_curv_3;
NDF_of_curvatures_3 (:,c) = sum(NDF_of_curvatures_3,2);
PDF_curvature_3(:,c) = NDF_of_curvatures_3 / ( d_curve * sum(NDF_of_curvatures_3));
end

pdf_n_3=sum(PDF_curvature_3,2)./i_end;

filtpdf_3=filtfilt(b,a,pdf_n_3(2:end-1));
Axes_curv_3(:,1)=curv_axis(2:end-1); Axes_curv_3(:,2)=filtpdf_3;
plot(curv_axis2(2:end-1),filtpdf_3,'-r', 'LineWidth',2); hold on;

kapa_mat_3=cell2mat(kapa_all_3);
kapa_mat2_3 = kapa_mat_3(kapa_mat_3>=-1 & kapa_mat_3<=1);
radii_3=1./kapa_mat2_3;

kr_3(1,1) = mean(kapa_mat2_3);
kr_3(2,1)= median(kapa_mat2_3);
%kr_3(3,1)= mean(radii_3);
kr_3(3,1)= median(radii_3);

save kapa_mat_3_1000.mat;

kkrr(:,1)=kr_1; kkrr(:,2)=kr_2; kkrr(:,3)=kr_3; kkrr(:,4)=Axes_kr_t;
pdfs(:,1)=curv_axis2(2:end-1); pdfs(:,2)=filtpdf_1; pdfs(:,3)=filtpdf_2; pdfs(:,4)=filtpdf_3; pdfs(:,5)=filtpdf_t;

xlswrite('Avg k and r_1000',kkrr);

filename1=[case_name,'_PDFs.xlsx'];
xlswrite(filename1,pdfs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2- Flame Surface Density  Calculations

d_pix=8; %1.5 mm
d=d_pix*pix_size;
maap1=zeros(L,W);

for c=1:length(B_t)
for i = 1:length(B_t{c,3})
    for k=d_pix+1:length(B_t{c,3}{i,1})
    
    Lf(k-d_pix,1)=B_t{c,3}{i,1}(k-d_pix,1)+B_t{c,3}{i,1}(k-d_pix+1,1)+B_t{c,3}{i,1}(k-d_pix+2,1)+B_t{c,3}{i,1}(k-d_pix+3,1)+B_t{c,3}{i,1}(k-d_pix+4,1)+B_t{c,3}{i,1}(k-d_pix+5,1)+B_t{c,3}{i,1}(k-d_pix+6,1)+B_t{c,3}{i,1}(k-d_pix+7,1)+B_t{c,3}{i,1}(k-d_pix+8,1);
  
    sig{i,1}(k-d_pix,1)=Lf(k-d_pix,1)./(d*d);
    
    maap1(B_t{c,5}{i,1}(k-0.5*d_pix,1),B_t{c,7}{i,1}(k-0.5*d_pix,1))=maap1(B_t{c,5}{i,1}(k-0.5*d_pix,1),B_t{c,7}{i,1}(k-0.5*d_pix,1))+sig{i,1}(k-d_pix,1);
    end  
    
end 
Lf=[];
sig=[];
 end
maap=flipud(maap1);
maap2=maap./c;
% 
figure (2); imshow(maap2,[]); colormap jet; caxis([0 0.5]);
axis image;colormap(jet);colorbar;
title([case_name]);
saveas(gcf,[case_name,'_Map_fsd.png']);

Flame_Surface_Density_Map = maap2;
rad_pos=[Flame_Surface_Density_Map(400,:);Flame_Surface_Density_Map(266,:);Flame_Surface_Density_Map(133,:)];
rad_pos2=rot90(rad_pos,3);

filename2=[case_name,'_Map_FSD.xlsx'];
xlswrite (filename2,Flame_Surface_Density_Map);


filename3=[case_name,'_Diff_Lvls_FSD.xlsx'];
xlswrite (filename3,rad_pos2);

meanfsd=mean(Flame_Surface_Density_Map);
mean_fsd=rot90(meanfsd,3);
xax=[1:520];
xfsd=rot90(pix_size*(xax-259),3);

figure (3); plot(xfsd,rad_pos(1,:)); hold on; plot(xfsd,rad_pos(2,:));hold on;plot(xfsd,rad_pos(3,:));
title([case_name]); 
legend({'Z/D=0.25','Z/D=0.5','Z/D=0.75'},'Location','northwest');
saveas(gcf,[case_name,'_Diff_Lvls_FSD.png']);

Axes_fsd(:,1)=xfsd; Axes_fsd(:,2)=mean_fsd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3- Flame Front Angle Calculations

curv_limits_theta = [-90 90] ;
curv_axis_theta = linspace(curv_limits_theta(1), curv_limits_theta(2), number_of_bins);
d_curve_theta = curv_axis_theta(2) - curv_axis_theta(1) ;
 
for c=1:length(B_t)
for i = 1:length(B_t{c,9})
    for k=d_pix_c+1:length(B_t{c,9}{i,1})
  
    theta{i,1}(k-d_pix_c,1)=atand(B_t{c,10}{i,1}(k-d_pix_c,1)./B_t{c,9}{i,1}(k,1)); 
   % theta is the angle with the +ve horizontal axis
    
    end  
    
end
theta1=cell2mat(theta);
theta_all{c,:}=theta1;
theta=[];
NDF_theta=[];
NDF_theta (:,c) = hist(theta1,curv_axis_theta);
NDF_of_theta = NDF_theta;
NDF_of_theta (:,c) = sum(NDF_of_theta,2);
PDF_theta(:,c) = NDF_of_theta / ( d_curve_theta * sum(NDF_of_theta));
end
pdf_n_theta=sum(PDF_theta,2)./i_end;
figure (4); plot(curv_axis_theta, pdf_n_theta); hold on;

filtpdf_theta=filtfilt(b,a,pdf_n_theta(1:end));
Axes_theta(:,1)=curv_axis_theta(1:end); Axes_theta(:,2)=filtpdf_theta;
plot(curv_axis_theta(1:end),filtpdf_theta,'-r', 'LineWidth',2); 

theta_mat=cell2mat(theta_all);
theta_avg = mean(theta_mat);

save theta_mat_1000.mat;