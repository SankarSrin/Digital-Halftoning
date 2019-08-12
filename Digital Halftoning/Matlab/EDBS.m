%Programmed By: Sankarasrinivasan S
%Multimedia Signal Processing Lab, Dept. of Elec Engg. NTUST
%Oct 2016

%Efficient Direct Binary search halftoning

function [dst]=EDBS(im)
im=double(im);
%Initiation
[rows cols]=size(im);
fs=7;    %Gaussian filter size
%Random MAtrix
dst=double(randn(rows,cols)>0.5);
%Gaussian Filter
%Generation of Gaussian Filter
d=(fs-1)/6;
gaulen=(fs-1)/2;
for k=-gaulen:gaulen
    for l=-gaulen:gaulen
        c=(k*k + l*l)/(2*d*d);
        GF(k+gaulen+1,l+gaulen+1)=exp(-c)/(2*3.14*d*d);
    end
end
CPP=zeros(13,13);
HalfCPPSize=6;
CPP=CPP+(conv2(GF,GF)); %AutoCorrelation of gaussian filter
imr=double(im/255);
Err=dst-imr;
CEP=xcorr2(Err,CPP);   %Cross Correlation between Error and Gaussian

double EPS=0; double EPS_MIN=0;  CountB=0;

while(1)
    CountB=0;
   double a0=0; double a1=0; double a0c=0; double a1c=0;
   uint8 Cpx=0; uint8 Cpy=0;
    
    for i=1:1:rows
        for j=1:1:cols
         a0c=0; a1c=0; Cpx=0; Cpy=0; EPS_MIN=0;
         for y=-1:1:1
             if (i+y < 1|| i+y >rows)
						continue;			
             end
             for x=-1:1:1
                 if (j+x < 1 || j+x >cols)
							continue;
                 end
                 if(y==0 && x==0)
                     if(dst(i,j)==1)
                         a0=-1; a1=0;
                     else
                         a0=1; a1=0;
                     end
                     else	%% Swap
							if (dst((i+y), (j+x))~= dst(i,j))
								if (dst(i ,j) == 1)	
									a0 = -1;	
									a1 = -a0;
                                else	
									a0 = 1;
									a1 = -a0;
                                end
                            else
								a0 = 0;
								a1 = 0;
                            end
                 end                 
EPS =(a0*a0+a1*a1)*CPP(HalfCPPSize+1,HalfCPPSize+1)+2*a0*a1*CPP(HalfCPPSize+y+1,HalfCPPSize+x+1)+2*a0*CEP(i+ HalfCPPSize,j+HalfCPPSize)+2*a1*CEP(i+y+HalfCPPSize,j+x+HalfCPPSize);     
if (EPS_MIN > EPS)
    EPS_MIN = EPS; 	a0c = a0;	a1c = a1;	Cpx = x;	Cpy = y;
end

 end             
         end
            if(EPS_MIN<0)
                 for y=-HalfCPPSize:1:HalfCPPSize
                     for x=-HalfCPPSize:1:HalfCPPSize
                         CEP(i+y+HalfCPPSize,j+x+HalfCPPSize) =  CEP(i+y+HalfCPPSize,j+x+HalfCPPSize)+ a0c*CPP(y+HalfCPPSize+1,x+HalfCPPSize+1);
                     end
                 end
                 for y=-HalfCPPSize:1:HalfCPPSize
                     for x=-HalfCPPSize:1:HalfCPPSize
                         CEP(i+y+Cpy+HalfCPPSize,j+x+Cpx+HalfCPPSize) = CEP(i+y+Cpy+HalfCPPSize,j+x+Cpx+HalfCPPSize)+a1c*CPP(y+HalfCPPSize+1,x+HalfCPPSize+1);
                     end
                 end
                 dst(i,j)=dst(i,j)+ (a0c);
                 dst(Cpy+i,j+Cpx)=dst(i+Cpy,j+Cpx)+ (a1c);
                 CountB=CountB+1;
end
        end
    end
    if(CountB==0)
        break;
    end

end
dst=dst;
end



                 
                     
                     
                     
