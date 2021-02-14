% Name : hossein hayati 
% Student Number : 980042875


clc;
clear;
close all;

%%%%% Section 2 : Conventional image enhancement method %%%%%

% just load and plot an image for test
sample_img = imread('C:\Users\h256.DESKTOP-NJDT39C\Documents\projects\azad\image\DS\reference_images\I08.BMP');

sample_img = im2double(sample_img);

% constants
L = 256; % The magnitude of each and every color channel is confined within the range [0 , L-1]
M = 256;% total intensity levels ???? ( wrong ) 

% Extract color channels. REPEATED
redChannel = sample_img(:,:,1); % Red channel
greenChannel = sample_img(:,:,2); % Green channel
blueChannel = sample_img(:,:,3); % Blue channel


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
figure()
imhist(before_stretched);
figure()
imhist(after_stretched);
% Recombine separate color channels into an RGB image.
rgb_stretched_Image = cat(3, r_stretched, g_stretched, b_stretched);


%%%%% Convert RGB to HSI %%%%%
% The image can be enhanced by conserving H and S channels
% components while improving only intensity channel

hsi_image = rgb2hsi(rgb_stretched_Image);

intensity = hsi_image(:, :, 3);

figure()
[counts , binLocations] = imhist(intensity);
imhist(intensity);
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

figure()
plot(hist)
figure()
plot(clipped_hist)

%%

% (15) the corresponding PDF (p(i)) is calculated as:
% M is the total intensity levels   

% Pi = clipped_hist / M;
tpi = zeros(1,length_hist);
tpi_LUT = 1:L;

for i = 1:L
    P(i) = clipped_hist(i) / sum(clipped_hist);
end
for i = 1:L
    
    % (14) the CDF, c(i) is formulated as:
    C(i) = sum(P(1:i));         
end
for i = 1:L
    Alpha = C(i);

    % (13) the corresponding Weighted Histogram Distribution function (WHD) is
    % constructed as:
    % [ Pmin,Pmax : min and max PDF of clipped histogram ]
    
    Pmin = min(P);
    Pmax = max(P);

    Pw(i) = Pmax * power((P(i)-Pmin)/(Pmax-Pmin),Alpha);
end
for i = 1:L
    Alpha = C(i);

%     plot(Pw(1:i))

    % (12) the weighted PDF sum is defined as follows:
    %  [xL: because intensity is double between 0 and 1]  
    % intensity_max = max(max(intensity*L)); 
    intensity_max = max(max(intensity(1:i))); 
%     Sum_Pw(i) = sum(Pw(2:intensity_max));
%     Sum_Pw(i) = sum(Pw(1:intensity_max));
    Sum_Pw = sum(Pw(1:L));
   
    % (11) the weighted CDF is defined as:
    Cw(i) = sum(Pw(1:intensity_max)/Sum_Pw); 

    % (10) the gamma parameter with weighted CDF is computed as:
    gamma(i) = 1 - Cw(i);
    tpi_LUT(i) =  power(round(tpi_LUT(i)/i),gamma(i)) ;

    % (9) Transformed pixel intensity 
end
% tpi = intlut(intensity,tpi_LUT);
[length , width] = size(intensity) ;

result = zeros(length,width);

for i=1:L
    disp(i)
    for l = 1:length
        for w = 1:width
            if uint8(intensity(l,w)*255) == i
               result(l,w) =  power(round(intensity(l,w)/intensity_max),gamma(i));
            end
        end
    end
end
% test =  power(round(intensity/intensity_max),gamma) ;

H1 = hsi_image(:, :, 1);
S1 = hsi_image(:, :, 2);
I1 = result;
result_rgb = hsi2rgb(H1,S1,I1);


figure()
imshow(result_rgb)

