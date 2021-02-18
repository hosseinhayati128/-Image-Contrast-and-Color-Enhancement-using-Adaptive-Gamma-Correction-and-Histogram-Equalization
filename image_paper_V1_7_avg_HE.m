% Name : hossein hayati 


clc;
clear;
close all;

imagefiles = dir('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\distorted_images\*.bmp');      
nfiles = length(imagefiles);  
Sat_sum = 0;
Contrast_sum = 0;
Color_sum = 0;
Entropy_sum = 0;

for ii=1:nfiles
   disp((ii/nfiles)*100)
   currentfilename = imagefiles(ii).name;
   img_folder = imagefiles(ii).folder;
   sample_img = imread(fullfile(img_folder,currentfilename));
   sample_img = im2double(sample_img);
   HE_tst = histeq(sample_img);

% 4.1 Entropy:
Entropy = entropy(HE_tst);
Entropy_sum = Entropy_sum + Entropy;

end
disp("Entropy")
disp(Entropy_sum/nfiles)
