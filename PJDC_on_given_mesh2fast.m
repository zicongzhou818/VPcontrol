function  [phi1_comp,phi2_comp,U1,U2,phi_m1, phi_m2,U1_m,U2_m,Rot]=PJDC_on_given_mesh2fast(JD_prescribed,CV_prescribed,N,given_x,given_y)
%% initialization
better = 1;
ei = 0;
r = 1;
ti = 0;
ts = 1;
Npts = N-2;
imax = 5.0e4;
ts_r = 1e-6;
rTol = 1.0e-4;
Rot=[imax,1];
%% preset control functions corresponsive to the identity map
f = ones(Npts);
g = zeros(Npts);
%% preset new deformation based on control functions that are corresponsive to the identity map
[X,Y] = ndgrid(1:N,1:N);
phi_m1 = X;
phi_m2 = Y;
%% reading prescriptions
JD_prescribed=JD_prescribed*(N^2/(sum(sum(JD_prescribed))));
f0=JD_prescribed(2:Npts+1,2:Npts+1);
g0=CV_prescribed(2:Npts+1,2:Npts+1);
%% reading given mesh 
h=1;
phi_o1=given_x;
phi_o2=given_y;
[phi_o1_y,phi_o1_x]=gradient(phi_o1,h);
[phi_o2_y,phi_o2_x]=gradient(phi_o2,h);
%% (ps: the final phi initially is the same as given phi here because the current phi is id)
phi1_comp=given_x;
phi2_comp=given_y;
%% compute ssd
[phi_o_jd, phi_o_cv]=compute_JD_and_Curl(phi_o1,phi_o2,h);
Pn=(phi_o_jd(2:Npts+1,2:Npts+1)-f0);
Qn=(phi_o_cv(2:Npts+1,2:Npts+1)-g0);

ssd_initial=sum(sum(Pn.^2+Qn.^2));
ssd_old=ssd_initial;

%% main loop
while ei<imax && ts>ts_r && r>rTol
    ti=ti+1;
    if better % form gradients wrt F when ssd decreases
        ei=ei+1;
        [phi_jd, phi_cv]=compute_JD_and_Curl(phi1_comp, phi2_comp,h);
        [phi_cur1_y,phi_cur1_x]=gradient(phi_m1,h);
        [phi_cur2_y,phi_cur2_x]=gradient(phi_m2,h);
        
        P=(phi_jd - JD_prescribed).*phi_o_jd;
        Q=phi_cv - CV_prescribed;
        
        P_11=-P.*phi_cur2_y;
        Q_phi_pre1_y=Q.*phi_o1_y;
        b11=P_11+Q_phi_pre1_y;
                
        P_12=P.*phi_cur2_x;
        Q_phi_pre2_y=Q.*phi_o2_y;
        b12=P_12+Q_phi_pre2_y;
        
        P_21=P.*phi_cur1_y;
        Q_phi_pre1_x=-Q.*phi_o1_x;
        b21=P_21+Q_phi_pre1_x;

        P_22=-P.*phi_cur1_x;
        Q_phi_pre2_x=-Q.*phi_o2_x;
        b22=P_22+Q_phi_pre2_x;
        
        %% form B
        [~,b11_x]=gradient(b11,h);
        [b12_y,~]=gradient(b12,h);
        lap_B1=b11_x+b12_y;
        
        [~,b21_x]=gradient(b21,h);
        [b22_y,~]=gradient(b22,h);
        lap_B2=b21_x+b22_y;

        B1=pois2fft2(lap_B1);
        B2=pois2fft2(lap_B2);

        [a1_x2, a1_x1] = gradient(B1); 
        [a2_x2, a2_x1] = gradient(B2); 
        
        del_f = - (a1_x1 + a2_x2);
        del_g = - a2_x1 + a1_x2;
        
    end
    %update f1n,f2n to calculate u1,u2
    fn=f-ts*del_f(2:Npts+1,2:Npts+1);
    gn=g -ts*del_g(2:Npts+1,2:Npts+1);

    [f_x2, f_x1] = gradient(fn); 
    [g3_x2, g3_x1] = gradient(gn); 
    
    f1n = f_x1 - g3_x2;
    f2n = f_x2 + g3_x1;

    %update u1,u2
    u1=pois2fft2(f1n);
    u2=pois2fft2(f2n);
    U1_m = matrixpad(u1,0);
    U2_m = matrixpad(u2,0);
    %update T_new to compute ssd

    phi_m1=X+U1_m;
    phi_m2=Y+U2_m;
    phi1_comp = interpn(X, Y, phi_o1, phi_m1, phi_m2, 'makima'); 
    phi2_comp = interpn(X, Y, phi_o2, phi_m1, phi_m2, 'makima'); 

    %update JT,curlT
    [phi_jd, phi_cv]=compute_JD_and_Curl(phi1_comp,phi2_comp,h);
    Pn=(phi_jd(2:Npts+1,2:Npts+1)-f0);
    Qn=(phi_cv(2:Npts+1,2:Npts+1)-g0);
    ssd_new=sum(sum(Pn.^2+Qn.^2));
    r=ssd_new/ssd_initial;
    % if ssd decrease, then update f,T, increase tstep
    % if not, decrese tstep
%     display([' ts:',num2str(ts),' r:',num2str(r), ' ei:',num2str(ei), ' ti:',num2str(ti)]);
    Rot(ti,1)=r;
    if (ssd_new<ssd_old)
        ts=ts*1.1;
        f=fn;
        g=gn;
        ssd_old=ssd_new;
        better=1;
    else
        better=0;
        ts=ts*0.9;
    end 
end
display(['ts: ',num2str(ts),' r: ',num2str(r), ' ei: ',num2str(ei), ' ti: ',num2str(ti)]);
U1 = phi1_comp-X;
U2 = phi2_comp-Y;
Rot=Rot(1:ti,1);


end
