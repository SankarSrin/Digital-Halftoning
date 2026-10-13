% Computes reference outputs of the MATLAB/C code for tests/test_parity.py.
% octave-cli ref.m <matlabdir> <mexdir> <ppdir> <out.mat>   (see README.md)
a=argv(); addpath(a{1}); addpath(a{2}); addpath(a{3}); pkg load image; pkg load signal;
load(fullfile(fileparts(mfilename('fullpath')),'inputs.mat'));
X=logical(X); Xc=logical(Xc); Xr=logical(Xr); X3=logical(X3);
R=struct(); 
function R=tryset(R,name,f)
  try, R.(name)=f(); catch e, R.(name)=['ERR: ' e.message]; end
end
imgs={'pep','crop','ramp','rgb'}; 
for k=1:numel(imgs)
  I=eval(imgs{k}); n=imgs{k};
  R=tryset(R,['bayers1_' n],@() double(Bayers(I,1)));
  R=tryset(R,['bayers2_' n],@() double(Bayers(I,2)));
  for m=1:3, R=tryset(R,sprintf('ulichney%d_%s',m,n),@() double(Ulichney(I,m))); end
  for m=1:3, R=tryset(R,sprintf('dotdiff%d_%s',m,n),@() double(dotdiff(I,m))); end
  R=tryset(R,['floyd_' n],@() double(floyd_serp(I)));
end
J=[0 0 -99*48 7 5; 3 5 7 5 3; 1 3 5 3 1]/48; A=[0 0 -99*16 7; 1 3 5 0]/16;
c=im2double(crop);
[o,q,kk]=errdiff(c,J,0,1,0); R.ed_jarvis_raster=o; R.ed_jarvis_raster_q=q; R.ed_jarvis_raster_k=kk;
[o,q,kk]=errdiff(c,J,0,-1,0); R.ed_jarvis_serp=o; R.ed_jarvis_serp_k=kk;
R=tryset(R,'ed_asym_serp',@() errdiff(c,A,0,-1,0));
R=tryset(R,'ed_asym_raster',@() errdiff(c,A,0,1,0));
[o,q,kk]=errdiff(c,[0 -99*16 7; 3 5 1]/16,0.5,-1,0); R.ed_fs_l05=o; R.ed_fs_l05_k=kk;
R=tryset(R,'edbs_a',@() EDBS(dbs_a,init_a));
R=tryset(R,'edbs_b',@() EDBS(dbs_b,init_b));
% point process
R=tryset(R,'pp_X3',@() double(pp(double(X3))));
R=tryset(R,'pp_inv',@() double(pp(double(~X))));
R=tryset(R,'int_X3',@() intnsty(double(X3)));
R=tryset(R,'lau_bayer',@() double(bayer(dimg)));
D=double(X);
S={ {32}, {64,128,[16 32]}, {33}, {16,20}, {64,64,8}, {48,64} };
for k=1:numel(S)
  args=S{k};
  try, [p,f,cn]=rapsd(D,args{:}); R.(sprintf('rapsd%d',k))=p; R.(sprintf('rapsd%d_f',k))=f; R.(sprintf('rapsd%d_c',k))=cn; catch e, R.(sprintf('rapsd%d',k))=['ERR: ' e.message]; end
  try, [p,f,cn]=anis(D,args{:}); R.(sprintf('anis%d',k))=p; R.(sprintf('anis%d_c',k))=cn; catch e, R.(sprintf('anis%d',k))=['ERR: ' e.message]; end
end
try, [p,f]=rapsd(double(Xr),30); R.rapsd_rect=p; catch e, R.rapsd_rect=['ERR: ' e.message]; end
try, [g,r]=pc(X,1); R.pc1=g; R.pc1_r=r; catch e, R.pc1=['ERR: ' e.message]; end
R=tryset(R,'pc2',@() pc(X,0.5,30));
R=tryset(R,'pc3',@() pc(Xc,1,40,'clus'));
R=tryset(R,'pc4',@() pc(Xr,0.7));
R=tryset(R,'pc5',@() pc(Xr',1,100));
R=tryset(R,'pc6',@() pc(X,1,32,'xyz'));
R=tryset(R,'pc7',@() pc(Xc,2,60,'cluster'));
try, [g,r]=vpc(X,[1 1],0.5); R.vpc1=g; R.vpc1_r=r; catch e, R.vpc1=['ERR: ' e.message]; end
R=tryset(R,'vpc2',@() vpc(X3,[1 2],1));
R=tryset(R,'vpc3',@() vpc(X3,[3 1],1,20));
R=tryset(R,'vpc4',@() vpc(Xr,[1 1],0.75));
R=tryset(R,'rsmm1',@() rsmm(X));
R=tryset(R,'rsmm2',@() rsmm(Xc,20));
R=tryset(R,'rsmm3',@() rsmm(Xr));
try, [T,an]=ddf(X,[1 3 5],[3 5 8]); R.ddf1=T; R.ddf1_a=an; catch e, R.ddf1=['ERR: ' e.message]; end
R=tryset(R,'ddf2',@() ddf(Xc,[2 6],[6 12],36));
R=tryset(R,'ddf3',@() ddf(Xr,[0.5],[4],12));
save('-v7',a{4},'R');
