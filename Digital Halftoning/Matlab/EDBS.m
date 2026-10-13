%Programmed By: Sankarasrinivasan S
%Multimedia Signal Processing Lab, Dept. of Elec Engg. NTUST
%Oct 2016

%Efficient Direct Binary search halftoning
%
%   dst = EDBS(im)              im: grayscale (uint8 0..255 or double 0..1) or RGB
%   dst = EDBS(im, init)        start from a given binary halftone (0/1)
%   dst = EDBS(im, init, maxit) stop after maxit full passes (default Inf)
%
% Revision (2026):
%  - Default start pattern is rand(rows,cols)>0.5 (a 50% random pattern).
%    The old randn(rows,cols)>0.5 gave a pattern that was only ~31% ones,
%    probably a typo for rand.
%  - Input converted with im2double (old code did double(im)/255, which
%    silently gave an all-black result for a double image in 0..1).
%  - Removed the "double a0=0; uint8 Cpx=0;" lines: in MATLAB these are
%    command syntax calls (double('a0=0')) and declared nothing.
%  - Optional start pattern and pass limit; RGB input accepted.
%  The DBS cost update itself is unchanged.

function [dst]=EDBS(im,init,maxit)
if nargin<3, maxit=Inf; end
imr=im2double(im);
if size(imr,3)==3
    imr=rgb2gray(imr);
end
%Initiation
[rows,cols]=size(imr);
fs=7;    %Gaussian filter size
if nargin<2 || isempty(init)
    dst=double(rand(rows,cols)>0.5);   %Random start pattern
else
    dst=double(init~=0);
end
%Gaussian Filter
%Generation of Gaussian Filter
d=(fs-1)/6;
gaulen=(fs-1)/2;
GF=zeros(fs,fs);
for k=-gaulen:gaulen
    for l=-gaulen:gaulen
        c=(k*k + l*l)/(2*d*d);
        GF(k+gaulen+1,l+gaulen+1)=exp(-c)/(2*3.14*d*d);
    end
end
HalfCPPSize=6;
CPP=conv2(GF,GF); %AutoCorrelation of gaussian filter (13x13)
Err=dst-imr;
CEP=conv2(Err,rot90(CPP,2));   %Cross Correlation between Error and Gaussian (= xcorr2(Err,CPP), no toolbox needed)

iter=0;
while(iter<maxit)
    iter=iter+1;
    CountB=0;
    for i=1:1:rows
        for j=1:1:cols
            a0c=0; a1c=0; Cpx=0; Cpy=0; EPS_MIN=0;
            for y=-1:1:1
                if (i+y < 1 || i+y > rows)
                    continue;
                end
                for x=-1:1:1
                    if (j+x < 1 || j+x > cols)
                        continue;
                    end
                    if(y==0 && x==0)            %% Toggle
                        if(dst(i,j)==1)
                            a0=-1; a1=0;
                        else
                            a0=1; a1=0;
                        end
                    else                        %% Swap
                        if (dst(i+y,j+x) ~= dst(i,j))
                            if (dst(i,j) == 1)
                                a0 = -1;
                            else
                                a0 = 1;
                            end
                            a1 = -a0;
                        else
                            a0 = 0;
                            a1 = 0;
                        end
                    end
                    EPS =(a0*a0+a1*a1)*CPP(HalfCPPSize+1,HalfCPPSize+1)+2*a0*a1*CPP(HalfCPPSize+y+1,HalfCPPSize+x+1)+2*a0*CEP(i+HalfCPPSize,j+HalfCPPSize)+2*a1*CEP(i+y+HalfCPPSize,j+x+HalfCPPSize);
                    if (EPS_MIN > EPS)
                        EPS_MIN = EPS; a0c = a0; a1c = a1; Cpx = x; Cpy = y;
                    end
                end
            end
            if(EPS_MIN<0)
                for y=-HalfCPPSize:1:HalfCPPSize
                    for x=-HalfCPPSize:1:HalfCPPSize
                        CEP(i+y+HalfCPPSize,j+x+HalfCPPSize) = CEP(i+y+HalfCPPSize,j+x+HalfCPPSize)+ a0c*CPP(y+HalfCPPSize+1,x+HalfCPPSize+1);
                    end
                end
                for y=-HalfCPPSize:1:HalfCPPSize
                    for x=-HalfCPPSize:1:HalfCPPSize
                        CEP(i+y+Cpy+HalfCPPSize,j+x+Cpx+HalfCPPSize) = CEP(i+y+Cpy+HalfCPPSize,j+x+Cpx+HalfCPPSize)+a1c*CPP(y+HalfCPPSize+1,x+HalfCPPSize+1);
                    end
                end
                dst(i,j)=dst(i,j)+ (a0c);
                dst(i+Cpy,j+Cpx)=dst(i+Cpy,j+Cpx)+ (a1c);
                CountB=CountB+1;
            end
        end
    end
    if(CountB==0)
        break;
    end
end
end
