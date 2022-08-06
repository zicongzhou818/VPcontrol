clc
clear;
close all;
format long
%% 
N=129;h=1;h1=h*2;h2=h*4;
m=N;n=N;
x=1:n;
y=1:m;
[X,Y]=ndgrid(x,y);
White_bg = zeros(N, N)+255;
%% load problem
T_original=double(rgb2gray(imread('I_little.png')));
T=imresize(T_original,[m,n]);
R_original=double(rgb2gray(imread('I_large.png')));
R=imresize(R_original,[m,n]);
%% fwd registration
h_step=8;h_0=1;
%% these are the parameters tuned for J to V and V to J example
jminbnd = 0.6;
jmaxbnd = 5.8;
ratio_bnd = 0.02;
%%
tic
[phix_T2R,phiy_T2R,Tn_T2R,U1n_T2R,U2n_T2R]=Img_Reg2D_diffeo_regriddingPluspostf(T,R,N);
toc
%  r:0.06053 maxJ:2.7622 minJ:0.25401
% Elapsed time is 13.706418 seconds.
[jsc, dice, Jmax, Jmin] = evaluations2(Tn_T2R, R, phix_T2R,phiy_T2R)
% jsc =   0.984992965452556
% dice =   0.992439754292015
% Jmax =   2.762243264183202
% Jmin =   0.254008128721549

figure(101)
imshow(imresize(T, [512 512]),[]), title('I_m'); 
figure(102)
imshow(imresize(R, [512 512]),[]), title('I_f'); 
figure(103)
imshow(ones(512, 512)*255,[]), title('\phi'), hold on
for i = 1:h:N
    plot(phiy_T2R(i,1:+h:end),phix_T2R(i,1:+h:end),'r-'), hold on
    plot(phiy_T2R(1:+h:end,i),phix_T2R(1:+h:end,i),'r-'), hold on
end
axis ([0,N+1,0,N+1])

figure(104)
imshow(ones(512, 512)*255,[]), title('\phi'), hold on
quiver(Y(1:+h2:end,1:+h2:end),X(1:+h2:end,1:+h2:end),U2n_T2R(1:+h2:end,1:+h2:end),U1n_T2R(1:+h2:end,1:+h2:end),0);
axis ([0,N+1,0,N+1])
figure(105)
imshow(imresize(Tn_T2R, [512 512]),[]), title('I_m(\phi)');

%% guess inverse
tic
[phi1_comp,phi2_comp,U1_comp,U2_comp,phi_m1, phi_m2,U1_m,U2_m,ratio1]=PJDC_on_given_mesh2fast(ones(N, N),zeros(N, N),N,phix_T2R,phiy_T2R);
toc
% ts: 1.16 r: 9.9756e-05 ei: 132 ti: 250
% Elapsed time is 3.449949 seconds.
guess_inv = interpn(X, Y, R, phi_m1, phi_m2, 'makima'); 

figure(106)
imshow(ones(512, 512)*255,[]), title('\phi^{-1}'), hold on
for i = 1:h:N
    plot(phi_m2(i,1:+h:end),phi_m1(i,1:+h:end),'r-'), hold on
    plot(phi_m2(1:+h:end,i),phi_m1(1:+h:end,i),'r-'), hold on
end
axis ([0,N+1,0,N+1])
figure(107)
imshow(ones(512, 512)*255,[]), title('\phi^{-1} o \phi vs id'), hold on
for i = 1:h:N
    plot(Y(i,1:+h:end),X(i,1:+h:end),'k-'), hold on
    plot(Y(1:+h:end,i),X(1:+h:end,i),'k-'), hold on
end
for i = 1:h:N
    plot(phi2_comp(i,1:+h:end),phi1_comp(i,1:+h:end),'r-'), hold on
    plot(phi2_comp(1:+h:end,i),phi1_comp(1:+h:end,i),'r-'), hold on
end
axis ([0,N+1,0,N+1])
figure(108)
imshow(ones(512, 512)*255,[]), title('v'), hold on
quiver(Y(1:+h2:end,1:+h2:end),X(1:+h2:end,1:+h2:end),U2_m(1:+h2:end,1:+h2:end),U1_m(1:+h2:end,1:+h2:end),0);
axis ([0,N+1,0,N+1])
figure(109)
imshow(imresize(guess_inv, [512 512]),[]), title('I_f(\phi^{-1})');
figure(110)
imshow(abs(imresize(guess_inv-T, [512 512])),[]);

[jsc, dice, Jmax, Jmin] = evaluations2(guess_inv, T, phi_m1, phi_m2)
ssd_inv=sum(sum((guess_inv-T).^2))
r_inv=ssd_inv/sum(sum((R-T).^2))
% jsc =   0.983072100313480
% dice =   0.991463800189693
% Jmax =   2.843619819814246
% Jmin =   0.321114198044457
% ssd_inv =     5.264603839790616e+05
% r_inv =   0.065039903539434

