function RGB1 = hsi2rgb(H1,S1,I1)
 %Multiply Hue by 360 to represent in the range [0 360]  
 H1=H1*360;                                               
    


 %Preallocate the R,G and B components  
 R1=zeros(size(H1));  
 G1=zeros(size(H1));  
 B1=zeros(size(H1));  
 RGB1=zeros([size(H1),3]);  
    


 %RG Sector(0<=H<120)  
 %When H is in the above sector, the RGB components equations are  


    
 B1(H1<120)=I1(H1<120).*(1-S1(H1<120));  
 R1(H1<120)=I1(H1<120).*(1+((S1(H1<120).*cosd(H1(H1<120)))./cosd(60-H1(H1<120))));  
 G1(H1<120)=3.*I1(H1<120)-(R1(H1<120)+B1(H1<120));  


    
 %GB Sector(120<=H<240)  
 %When H is in the above sector, the RGB components equations are  


    
 %Subtract 120 from Hue  
 H2=H1-120;  


    
 R1(H1>=120&H1<240)=I1(H1>=120&H1<240).*(1-S1(H1>=120&H1<240));  
 G1(H1>=120&H1<240)=I1(H1>=120&H1<240).*(1+((S1(H1>=120&H1<240).*cosd(H2(H1>=120&H1<240)))./cosd(60-H2(H1>=120&H1<240))));  
 B1(H1>=120&H1<240)=3.*I1(H1>=120&H1<240)-(R1(H1>=120&H1<240)+G1(H1>=120&H1<240));  


    
 %BR Sector(240<=H<=360)  
 %When H is in the above sector, the RGB components equations are  


    
 %Subtract 240 from Hue  
 H2=H1-240;  


    
 G1(H1>=240&H1<=360)=I1(H1>=240&H1<=360).*(1-S1(H1>=240&H1<=360));  
 B1(H1>=240&H1<=360)=I1(H1>=240&H1<=360).*(1+((S1(H1>=240&H1<=360).*cosd(H2(H1>=240&H1<=360)))./cosd(60-H2(H1>=240&H1<=360))));  
 R1(H1>=240&H1<=360)=3.*I1(H1>=240&H1<=360)-(G1(H1>=240&H1<=360)+B1(H1>=240&H1<=360));  


    
 %Form RGB Image  
 RGB1(:,:,1)=R1;  
 RGB1(:,:,2)=G1;  
 RGB1(:,:,3)=B1;  


    
 %Represent the image in the range [0 255]  
 RGB1=im2uint8(RGB1);  
