

%Modified code from The University of Texas at Austin

function [out, qn, k] = errdiff(in, fc, l, dir, v)

%   [OUT, QN, K] = ERRDIFF(IN, FC, L, DIR, V) performs error
%	diffusion on image IN using error filter FC, modified parameter L,
%	and direction DIR, where 1 is raster scan and -1 is serpentine.
%	When V (verbose) is non-zero, progress is printed to the output.
%	Input range is 0 to 1, as is output range.  



if nargin<5					% default to verbose
  v=1; end
if nargin<4					% default to raster
  dir=1; end
if nargin<3					% default to unmodified
  l=0; end
if nargin<2					% default to Floyd-Steinberg
  fc=[0 -99*16 7; 3 5 1]/16; end

[ri,ci]=size(in);
[rm,cm]=size(fc);
[r0,c0]=find(fc==-99);				% find origin of error filter
fc(r0,c0)=0;

rm=rm-1; cm=cm-1;
inpad=zeros(ri+rm,ci+cm);			% modified input image
inpad(r0:r0+ri-1,c0:c0+ci-1)=in;
out=zeros(ri,ci); qn=out;
sp=1; ep=ci; step=1;				% for direction changing
r0=r0-1; c0=c0-1;

for y=1:ri
  for x=sp:step:ep
    inpix=inpad(y+r0,x+c0);
    outpix=(inpix+l*in(y,x))>=0.5;
    out(y,x)=outpix;
    qerr=outpix-inpix;
    qn(y,x)=qerr;
    inpad(y:y+rm,x:x+cm)=inpad(y:y+rm,x:x+cm)-qerr*fc;
  end
  if dir==-1
    step=-step; temp=ep; ep=sp; sp=temp;
    fc=fc(:,cm+1:-1:1); end
  if v
   % fprintf('\rDithering... %3d%% done',round(y/ri*100)), 
   end
end

if v
  fprintf('\n')
end

if nargout==3
  xp=out(:)-0.5-qn(:);
  k=sum(abs(xp))/(2*sum(xp.^2));
end
