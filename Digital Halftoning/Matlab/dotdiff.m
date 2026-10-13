%Function to generate halftoning by Dot Diffusion MEthod
%Ref- Digital Halftoning by Donald E Knuth 1987

%im-input image (grayscale or RGB, any size)
%met 1-MeseVaidhayanath CM 8*8 2-Knuth CM 8*8 else Mese CM 16*16

%Programmed By: Sankarasrinivasan S
%Multimedia Signal Processing Lab, Dept. of Elec Engg. NTUST
%Oct 2016
%
% Revision (2026):
%  - The quantisation error is now divided by the total weight of the
%    neighbours that still receive error (class > current class). The old
%    code divided by sum(DM(:)) = 1, so part of the error was thrown away and
%    the output tone was wrong (e.g. 25% gray -> 12-15% ink).
%    If no neighbour has a higher class, the error is dropped (Knuth's
%    "baron" case) instead of dividing by zero.
%  - Knuth class matrix: the entry at row 7, column 8 is 44 (was a duplicate 42).
%  - Fixed CM1(2:s1+1,2:s1+1) -> CM1(2:s1+1,2:s2+1) (non-square images crashed).
%  - Any image size (class matrix tiled and cropped); no Image Processing
%    Toolbox needed (padarray replaced by explicit zero padding).
%  - Simplified the pixel loop (the old nested while loops visited each pixel
%    of a class exactly once; this one does the same).


function [H]=dotdiff(im,met)

if nargin<2, met=1; end
im=im2double(im);
if size(im,3)==3
    im=rgb2gray(im);
end

[s1,s2]=size(im);
impad=zeros(s1+2,s2+2);
impad(2:s1+1,2:s2+1)=im;
im=impad;

if(met==1)
    
              CM = [47 31	51	24	27	45	5	21;
                    37	63	53	11	22	4	1	33;
                    61	0	57	16	26	29	46	8;
                    20	14	9	62	18	41	38	6;
                    17	13	25	15	55	48	52	58;
                    3	7	2	32	30	34	56	60;
                    28	40	36	39	49	43	35	10;
                    54	23	50	12	42	59	44	19;];  %Messe-Vaidhya
                bs=8; nn=63;
elseif(met==2)

               CM=  [34 48 40 32 29 15 23 31;    
				     42 58 56 53 21 5 7 10;
				     50 62 61 45 13 1 2 18;
				     38 46 54 37 25 17 9 26;
				     28 14 22 30 35 49 41 33;
				     20 4 6 11 43 59 57 52;
				     12 0 3 19 51 63 60 44;
				     24 16 8 27 39 47 55 36;]; %Knuth
                 bs=8; nn=63; 
else

CM= [207	0	13	17	28	55	18	102	81	97	74	144	149	169	170	172
3	6	23	36	56	50	65	87	145	130	137	158	182	184	195	221
7	14	24	37	67	69	86	5	106	152	150	165	183	192	224	1
15	26	43	53	51	101	115	131	139	136	166	119	208	223	226	4
22	39	52	71	84	103	164	135	157	173	113	190	222	225	227	16
40	85	72	83	104	117	167	133	168	180	200	219	231	228	12	21
47	120	54	105	123	132	146	176	179	202	220	230	245	2	20	41
76	73	127	109	138	134	178	181	206	196	229	244	246	19	42	49
80	99	112	147	142	171	177	203	218	232	243	248	247	33	48	68
108	107	140	143	185	163	204	217	233	242	249	255	44	45	70	79
110	141	88	75	175	205	214	234	241	250	254	38	46	77	116	100
111	148	160	174	201	215	235	240	251	252	253	61	62	93	94	125
151	159	189	199	197	216	236	239	25	31	60	82	92	95	124	114
156	188	191	209	213	237	238	29	32	59	64	91	118	78	128	155
187	194	198	212	9	10	30	35	58	63	90	96	122	129	154	161
193	210	211	8	11	27	34	57	66	89	98	121	126	153	162	186];
bs=16; nn=255;

end


%Diffusion Matrix

DM=[1 2 1;
    2 0 2;
    1 2 1;];

CM=repmat(CM,ceil(s1/bs),ceil(s2/bs));
CM=CM(1:s1,1:s2);
CM1=zeros(s1+2,s2+2)-1;          % border gets class -1: never receives error
CM1(2:s1+1,2:s2+1)=CM;
CM=CM1;
H=zeros(s1+2,s2+2);

for ii=0:1:nn
    [p1,p2]=find(CM==ii);
    for k=1:numel(p1)
        r=p1(k); c=p2(k);
        if(im(r,c)>0.5)
            H(r,c)=1;
        end
        qerr=im(r,c)-H(r,c);
        w=DM.*(CM(r-1:r+1,c-1:c+1)>ii);   % only not-yet-processed neighbours
        sw=sum(w(:));
        if(sw>0)
            im(r-1:r+1,c-1:c+1)=im(r-1:r+1,c-1:c+1)+(qerr/sw)*w;
        end
    end
end
H=H(2:(s1+1),2:(s2+1));

end
