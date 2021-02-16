% Name : hossein hayati 

clc;
clear;
close all;

%%%%% Section 2 : Conventional image enhancement method %%%%%

% just load and plot an image for test
sample_img = imread('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\distorted_images\I19_17_4.bmp');
% sample_img = imread('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\reference_images\I19.bmp');

sample_img = im2double(sample_img);
figure('Name','input image double')
imshow(sample_img)
figure('Name','input image histogram')
imhist(sample_img)

% constants
L = 256; % The magnitude of each and every color channel is confined within the range [0 , L-1]
M = 256;% total intensity levels ???? ( wrong ) 


%%%%% 3. Proposed work %%%%%


% Extract color channels. REPEATED
redChannel = sample_img(:,:,1); % Red channel
greenChannel = sample_img(:,:,2); % Green channel
blueChannel = sample_img(:,:,3); % Blue channel


%%%%% 3.1 color channel stretching %%%%%

max_blue = max(max(blueChannel));
max_green = max(max(greenChannel));
max_red = max(max(redChannel));

min_blue = min(min(blueChannel));
min_green = min(min(greenChannel));
min_red = min(min(redChannel));

bn = blueChannel - min_blue;
rn = redChannel - min_red;
gn = greenChannel - min_green;

max_bn = max(max(bn));
max_rn = max(max(rn));
max_gn = max(max(gn));

b_stretched = bn/max_bn;
r_stretched = rn/max_rn;
g_stretched = gn/max_gn;
before_stretched = (blueChannel+greenChannel+redChannel)/3;
after_stretched = (b_stretched+r_stretched+g_stretched)/3;
% close all
% figure('Name','imhist before_stretched')
% imhist(before_stretched);
% figure('Name','imhist after_stretched')
% imhist(after_stretched);
% Recombine separate color channels into an RGB image.
rgb_stretched_Image = cat(3, r_stretched, g_stretched, b_stretched);


%%%%% Convert RGB to HSI %%%%%
% The image can be enhanced by conserving H and S channels
% components while improving only intensity channel

hsi_image = rgb2hsi(rgb_stretched_Image);

intensity = hsi_image(:, :, 3);

% figure('Name','imhist intensity')
[counts , binLocations] = imhist(intensity);
% imhist(intensity);
hist = counts;


% (17) the clipping limit is computed based on the mean value of the
% intensity 

Tc = mean(hist);

% (16) h (i) c is histogram clipping which is used to control the level of excessive
% contrast enhancement and it is constructed as
length_hist = length(hist);
clipped_hist = zeros(length_hist,1);
for hist_id = 1:length_hist
    if hist(hist_id)<Tc
        clipped_hist(hist_id) = hist(hist_id);
        continue
    end
    clipped_hist(hist_id) = Tc;
end
% 
% figure('Name','imhist hist')
% plot(hist)
% figure('Name','imhist clipped_hist')
% plot(clipped_hist)

%%

% (15) the corresponding PDF (p(i)) is calculated as:
% M is the total intensity levels   
P = clipped_hist / sum(clipped_hist);  %         QQQ
% P = clipped_hist / M;

% (14) the CDF, c(i) is formulated as:
C = cumsum(P);

% (13) the corresponding Weighted Histogram Distribution function (WHD) is
% constructed as:
% [ Pmin,Pmax : min and max PDF of clipped histogram ]
Pmin = min(P);
Pmax = max(P);
Pw = Pmax * power((P-Pmin)/(Pmax-Pmin),C);

Cw = zeros(1,L);
gamma = zeros(1,L);

% (12) the weighted PDF sum is defined as follows:
Sum_Pw = sum(Pw(1:L));

% (11) the weighted CDF is defined as:
Cw = cumsum(Pw)/Sum_Pw;

% (10) the gamma parameter with weighted CDF is computed as:

gamma = 1 - Cw;
[Length , Width] = size(intensity) ;
result = zeros(Length,Width);
intensity_max = max(intensity(:));
% intensity_max = max(max(intensity(:))); 
% intensity = uint8(intensity*255);


% (9) Transformed pixel intensity 
for i=1:L
    disp(i)
    for l = 1:Length
        for w = 1:Width
            if uint8(intensity(l,w)*255+1) == i
               result(l,w) =  power((intensity(l,w)/intensity_max),gamma(i));
            end
        end
    end
end
% test =  power(round(intensity/intensity_max),gamma) ;

H1 = hsi_image(:, :, 1);
S1 = hsi_image(:, :, 2);
I1 = result;
result_rgb = hsi2rgb(H1,S1,I1);
result_rgb = im2double(result_rgb);
figure('Name','imhist rgb output')

imhist(result_rgb)

figure('Name','rgb output')
imshow(result_rgb)

% 4.1 Entropy:
% result_rgb = imread('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\reference_images\I19.bmp');
% result_rgb = im2double(result_rgb);

Entropy = entropy(result_rgb);
disp("Entropy")
disp(Entropy)
