%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
% Point Process                                                     %
%		This function converts a binary halftone to a sample of a   %
%		point process. If the input image is a color halftone, pp   %
%		returns a point process of each channel.                    %
%                                                                   %
% Synopsis:                                                         %
%       Y=pp(X)                                                     %
%			X=Binary halftone                                       %
%           Y=Point process sample of type logical.                  %
%                                                                   %
% Dr. Daniel Leo Lau                                                %
% Copyright June 13, 2000                                           %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function Y=pp(X)

if (any(X(:)~=0 & X(:)~=1))
	error('X must contain only zeros and ones!');
end;

% Revision (2026): returns a LOGICAL array. uint8(logical(...)) is plain
% uint8 since MATLAB 6.5, which pc/rsmm/ddf reject ("must be LOGICAL").
% mean2 (Image Processing Toolbox) replaced by mean(...(:)).
[M,N,C]=size(X);
Y=false(size(X));
for c=1:C
	Xc=X(:,:,c)~=0;
	if (mean(Xc(:))<=0.5)
		Y(:,:,c)=Xc;
	else
		Y(:,:,c)=~Xc;
	end;
end;
