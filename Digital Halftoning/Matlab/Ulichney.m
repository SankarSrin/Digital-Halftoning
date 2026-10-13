%%%%%%%Function to generate Ulichney Halftone Patterns%%%%%%%%%%

%%   Author:
%%   Mr. Sankarasrinivasan
%%   Research Fellow, Multimedia Signal Processing Lab, NTUST, Taiwan
%%   Advisor: Prof. Jing Ming Guo
%%   Dated: Apr, 2018

%met=1 (clustered dot 4x4); met=2 (clustered dot 8x8) met=3 dispersed dot
%(14x14 void-and-cluster)
%
% Revision (2026):
%  - works for any image size: the threshold array is tiled and cropped
%    (met=3 no longer resizes the image with imresize, so the output has the
%    same size as the input).
%  - met=1 thresholds are now 0..15/16 (were 0..15/15, so a pure-white input
%    still printed one black dot per 4x4 cell).
%  - accepts RGB input; errors on an unknown met.

function [HOD]=Ulichney(im,met)
if nargin<2, met=1; end
im=im2double(im);
if size(im,3)==3
    im=rgb2gray(im);
end
[s1,s2]=size(im);


if(met==1)

od=[12 05 06 13;
    04 00 01 07;
    11 03 02 08;
    15 10 09 14;]/16;

bs=4;

elseif(met==2)

 od=[0.567	0.635	0.608	0.514	0.424	0.365	0.392	0.486
    0.847	0.878	0.910	0.698	0.153	0.122	0.090	0.302
    0.820	0.969	0.941	0.667	0.180	0.031	0.059	0.333
    0.725	0.788	0.757	0.545	0.275	0.212	0.243	0.455
    0.424	0.365	0.392	0.486	0.567	0.635	0.608	0.514
    0.153	0.122	0.090	0.302	0.847	0.878	0.910	0.698
    0.180	0.031	0.059	0.333	0.820	0.969	0.941	0.667
    0.275	0.212	0.243	0.455	0.725	0.788	0.757	0.545];
bs=8;

elseif(met==3)
    od=[131 187 8 78 50 18 134 89 155 102 29 95 184 73;
        22 86 113 171 142 105 34 166 9 60 151 128 40 110;
        168 137 45 28 64 188 82 54 124 189 80 13 156 56;
        7 61 186 121 154 6 108 177 24 100 38 176 93 123;
        83 148 96 17 88 133 44 145 69 161 139 72 30 181;
        115 27 163 47 178 65 164 14 120 48 5 127 153 52;
        190 58 126 81 116 21 106 77 173 92 191 63 99 12;
        76 144 4 185 37 149 192 39 135 23 117 31 170 132;
        35 172 103 66 129 79 3 97 57 159 70 141 53 94;
        114 20 49 158 19 146 169 122 183 11 104 180 2 165;
        152 87 182 118 91 42 67 25 84 147 43 85 125 68;
        16 136 71 10 193 112 160 138 51 111 162 26 194 46;
        174 107 41 143 33 74 1 101 195 15 75 140 109 90;
        32 62 157 98 167 119 179 59 36 130 175 55 0 150;]/196;
    bs=14;

else
    error('Ulichney: met must be 1, 2 or 3');
end

mask=repmat(od,ceil(s1/bs),ceil(s2/bs));
mask=mask(1:s1,1:s2);
HOD=im>mask;

end
