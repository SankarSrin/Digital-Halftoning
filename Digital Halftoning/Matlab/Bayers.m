%%%%%%%Function to generate Bayer's Halftone Patterns%%%%%%%%%%

%%   Author: 
%%   Mr. Sankarasrinivasan
%%   Research Fellow, Multimedia Signal Processing Lab, NTUST, Taiwan
%%   Advisor: Prof. Jing Ming Guo
%%   Dated: Apr, 2018

%%Details
%%Bayers(im,met) im==> input image;
%% met=1; Bayers Dispersed Dot && met=2; Bayers Clustered Dot

function [HOD]=Bayers(im,met)
[s1 s2]=size(im);
im=im2double(im);
if (met==1)
%Bayers Dispersed Dot

DA=[00 48 12 60 03 51 15 63;
    32 16 44 28 35 19 47 31;
    08 56 04 52 11 59 07 55;
    40 24 36 20 43 27 39 23;
    02 50 14 62 01 49 13 61;
    34 18 46 30 33 17 45 29;
    10 58 06 54 09 57 05 53;
    42 26 38 22 41 25 37 21]/64;

else
%Bayers Clustered Dot

DA=[24 10 12 26 35 47 49 37;
    08 00 02 14 45 59 61 51;
    22 06 04 16 43 57 63 53;
    30 20 18 28 33 41 55 39;
    34 46 48 36 25 11 13 27;
    44 58 60 50 09 01 03 15;
    42 56 62 52 23 07 05 17;
    32 40 54 38 31 21 19 29;]/64;

end

mask=repmat(DA,s1/8,s2/8);
H=zeros(s1,s2);
HOD=im>mask;

end
