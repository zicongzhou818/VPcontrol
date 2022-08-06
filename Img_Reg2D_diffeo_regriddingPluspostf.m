function [pos_xkF,pos_ykF,T_tempk,U1F,U2F]=Img_Reg2D_diffeo_regriddingPluspostf(T,R,N)
T_test=T;
m=N;n=N;
% White_bg = zeros(m,n)+255;
ts = 1e-8;% tstep: initial tstep 
ts_r = 1e-10;% tstep_ratio: stop when tstep less than tstep_ratio
ts_up=1.05;% tstep_up: magnification factor when ssd decreases
ts_dn=0.95;% tstep_down: reduction factor when doesn't ssd decrease
ite_max=300;% ite_max: preset maximum number of iterations
rTol=1e-5;% ratiotol: preset ratio to stop 
r=1;
%% Step-1
grid_size=[m-2,n-2];% remove zero boundaries

F1=zeros(grid_size);
F2=zeros(grid_size);
% MM=(m-2)*(n-2);

% Compute initial ssd on uniform grid
ssd_old=sum(sum(((T_test-R).^2)));
ssd_initial=ssd_old;

% iteration flag
Better = 1;% Better: indicate whether the ssd decreases(Better=1) or not 
ti=0;% ii: total iteration steps

% set up initial values for Phi(x,y)=(x+u1(x,y),y+u2(x,y))
[xI,yI] = ndgrid(1:m,1:n);% (i,j)--->(pos_x(i,j),pos_y(i,j)) on (1 to m,1 to n)

pos_x0=xI;
pos_y0=yI;
pos_xk1=xI;
pos_yk1=yI;
%% new criterias for regridding
bd_J_min=0.75;
k=0;
ei=0;
while  (ts>ts_r) && (ei<ite_max) && (r>rTol)
    ti=ti+1; % record the number of iterations
    %% Step-2
    if Better  % --(PS: Better holds for the first step in "while" loop)
        ei=ei+1;
        
        % compute the difference from R to T (ps: not the same as from T to R)
        T_temp = interpn(xI, yI, T_test, pos_x0,pos_y0,'makima'); 
        TR_diff = T_temp-R; 
        [T_x2, T_x1] = gradient(T_test); 
        T_x1 = interpn(xI, yI, T_x1, pos_x0,pos_y0,'makima'); 
        T_x2 = interpn(xI, yI, T_x2, pos_x0,pos_y0,'makima'); 
        G1 = (TR_diff.*T_x1);
        G2 = (TR_diff.*T_x2);
        % remove boundaries of a, since a=0 on boundaries 
        lap_a1_interior=G1(2:m-1, 2:n-1);
        lap_a2_interior=G2(2:m-1, 2:n-1);
        % w1,w2 
        a1=pois2fft2(lap_a1_interior);
        a2=pois2fft2(lap_a2_interior);
    end
    % add boundaries back to get U, where U=0 on boundaries
    F1_new = F1 - a1 * ts;
    F2_new = F2 - a2 * ts;
    
    intU1=pois2fft2(F1_new);
    intU2=pois2fft2(F2_new);  
    U1k=matrixpad(intU1,0);
    U2k=matrixpad(intU2,0); 
    pos_xk=xI+U1k;
    pos_yk=yI+U2k;

    [JD, ~]=compute_JD_and_Curl(pos_xk,pos_yk,1);
    JD_min=min(min(JD));

    if JD_min<=bd_J_min %|| JD_max>=bd_J_max %|| ii>ite_max
        k=k+1;
        display([' ts:',num2str(ts),' r:',num2str(r),' k:',num2str(k),' ei:',num2str(ei),' ti:',num2str(ti)]);
        pos_xk1_temp=pos_xk1;
        pos_yk1_temp=pos_yk1;
        
        pos_xk1=interpn(xI, yI, pos_xk1_temp, pos_x0, pos_y0,'makima');
        pos_yk1=interpn(xI, yI, pos_yk1_temp, pos_x0, pos_y0,'makima');
        T_tempk = interpn(xI, yI, T, pos_xk1, pos_yk1,'makima');
        
        T_test = T_tempk;
        pos_x0=xI;
        pos_y0=yI;
        Better = 1;
        ts = 1e-8;
        F1=zeros(grid_size);
        F2=zeros(grid_size);

    else
        pos_x0=pos_xk;
        pos_y0=pos_yk;    
        T_temp = interpn(xI, yI, T_test, pos_x0,pos_y0,'makima');
        ssd=sum(sum(((T_temp-R).^2))); % use sums to approximate integrals
        r=ssd/ssd_initial; 
        if ssd>ssd_old
            ts=ts*ts_dn;% reduction factor when ssd doesn't decrease
            Better = 0;
        else
            ts=ts*ts_up;% magnification factor when ssd decreases
            ssd_old=ssd;
            Better = 1;
            F1=F1_new;
            F2=F2_new;
        end 
    end
    
end
pos_xk1_temp=pos_xk1;
pos_yk1_temp=pos_yk1;
pos_xk1=interpn(xI, yI, pos_xk1_temp, pos_x0, pos_y0,'makima');
pos_yk1=interpn(xI, yI, pos_yk1_temp, pos_x0, pos_y0,'makima');
T_tempk = interpn(xI, yI, T, pos_xk1,pos_yk1,'makima');

display(['ts:',num2str(ts),' r:',num2str(r),' ei:',num2str(ei),' ti:',num2str(ti)]);

T_2=T_tempk;
T_test=T_2;
Jmin=0.2;
f = ones(grid_size);
f1 = sqrt(f - Jmin);
g = zeros(grid_size);
ti=0;
ei=0;
ts = 1e-6;
pos_x0=xI;
pos_y0=yI;
pos_xk2=xI;
pos_yk2=yI;
Better = 1;
bd_J_min=0.6;
ite_max=5000;
ts_r = 1e-7;
while  (ts>ts_r)  && (r>rTol)  && (ei<ite_max)
    ti=ti+1; % record the number of iterations
    %% Step-2
    if Better  % --(PS: Better holds for the first step in "while" loop)
        ei=ei+1;
        % compute the difference from R to T (ps: not the same as from T to R)
        T_temp = interpn(xI, yI, T_test, pos_x0, pos_y0,'makima'); 
        TR_diff = T_temp-R; 
        [T_x2, T_x1] = gradient(T_temp); 
        % construct  lap(a)=[R(phi(X))-T(X)]*grad(R(phi(X)))
        lap_a1_interior = (TR_diff(2:m-1, 2:n-1).*T_x1(2:m-1, 2:n-1));
        lap_a2_interior = (TR_diff(2:m-1, 2:n-1).*T_x2(2:m-1, 2:n-1));
        % remove boundaries of a, since a=0 on boundaries 
        % w1,w2 
        a1=pois2fft2(lap_a1_interior);
        a2=pois2fft2(lap_a2_interior);
        
        [a1_x2, a1_x1] = gradient(a1); 
        [a2_x2, a2_x1] = gradient(a2); 
        
        del_f1 = -2* f1 .* (a1_x1 + a2_x2);
        del_g = -a2_x1+a1_x2;
    end
    
    f1_new = f1 - del_f1 * ts;
    f_new = Jmin + f1_new.^2;
    g_new = g - del_g * ts;
    
    %% Step-3 Update F
    [f_x2, f_x1] = gradient(f_new); 
    [g_x2, g_x1] = gradient(g_new); 
    
    F1 = f_x1 - g_x2;
    F2 = f_x2 + g_x1;
    
    %% Step-4 Solve for U, where lap(U)=F
    intU1=pois2fft2(F1);
    intU2=pois2fft2(F2);
    % add boundaries back to get U, where U=0 on boundaries
    U1k=matrixpad(intU1,0);
    U2k=matrixpad(intU2,0); 

    %% Step-6  % phi(x,y)=X+U
    pos_xk=xI+U1k;
    pos_yk=yI+U2k;

    [JD, ~]=compute_JD_and_Curl(pos_xk,pos_yk,1);
    JD_min=min(min(JD));

    if JD_min<=bd_J_min %|| JD_max>=bd_J_max %|| ii>ite_max
        k=k+1;
        display([' ts:',num2str(ts),' r:',num2str(r),' k:',num2str(k),' ei:',num2str(ei),' ti:',num2str(ti)]);
        pos_xk2_temp=pos_xk2;
        pos_yk2_temp=pos_yk2;
        
        pos_xk2=interpn(xI, yI, pos_xk2_temp, pos_x0,pos_y0,'makima');
        pos_yk2=interpn(xI, yI, pos_yk2_temp, pos_x0,pos_y0,'makima');
        T_tempk = interpn(xI, yI, T_2, pos_xk2, pos_yk2,'makima');
        
        T_test = T_tempk;
        pos_x0=xI;
        pos_y0=yI;
        Better = 1;
        f=ones(grid_size);
        f1 = sqrt(f - Jmin);
        g=zeros(grid_size);

    else
        pos_x0=pos_xk;
        pos_y0=pos_yk;    
        T_temp = interpn(xI, yI, T_test, pos_x0,pos_y0,'makima');
        ssd=sum(sum(((T_temp-R).^2))); % use sums to approximate integrals
        r=ssd/ssd_initial; 
        if ssd>ssd_old
            ts=ts*ts_dn;% reduction factor when ssd doesn't decrease
            Better = 0;
        else
            ts=ts*ts_up;% magnification factor when ssd decreases
            ssd_old=ssd;
            Better = 1;
            f1 = f1_new;
            g = g_new;
        end 
    end

end
pos_xk2_temp=pos_xk2;
pos_yk2_temp=pos_yk2;
pos_xk2=interpn(xI, yI, pos_xk2_temp, pos_x0, pos_y0,'makima');
pos_yk2=interpn(xI, yI, pos_yk2_temp, pos_x0, pos_y0,'makima');
pos_xkF=interpn(xI, yI, pos_xk1, pos_xk2, pos_yk2,'makima');
pos_ykF=interpn(xI, yI, pos_yk1, pos_xk2, pos_yk2,'makima');
T_tempk = interpn(xI, yI, T, pos_xkF, pos_ykF,'makima');


ssd=sum(sum(((T_tempk-R).^2))); % use sums to approximate integrals
r=ssd/ssd_initial;  

display(['ts:',num2str(ts),' r:',num2str(r),' ei:',num2str(ei),' ti:',num2str(ti)]);
U1F=pos_xkF-xI;
U2F=pos_ykF-yI;

[ff, ~] = compute_JD_and_Curl(pos_xkF,pos_ykF,1);
ha1 = max(max(ff));
ha2 = min(min(ff));

display([' r:',num2str(r),' maxJ:',num2str(ha1),' minJ:',num2str(ha2)]);
