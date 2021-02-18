clc 
clear 
close all 

imagefiles = dir('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\distorted_images\I19*.bmp');      
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles
%    results_dir_folder = 'C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\frh_result\';
   currentfilename = imagefiles(ii).name;
%    result_dir = fullfile(results_dir_folder,currentfilename);
   img_folder = imagefiles(ii).folder;
   currentimage = imread(fullfile(img_folder,currentfilename));
   f = figure('Name',currentfilename);
   imhist(currentimage)
   ylim([0,3000]);
   saveas(f, sprintf('%d.png',ii));
end
