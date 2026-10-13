%Programmed By: Sankarasrinivasan S
%Multimedia Signal Processing Lab, Dept. of Elec Engg. NTUST
%Oct 2016

function [y,q,k] = floyd_serp(in)

in=im2double(in);
if size(in,3)==3
    in=rgb2gray(in);
end
[y,q,k] = errdiff(in,[0 -99*16 7; 3 5 1]/16,0,-1,0); %FS but param -1 for Serpentine Scan, 0 = quiet
end
