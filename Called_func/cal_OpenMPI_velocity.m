function [velocity, velocity_norm] = cal_OpenMPI_velocity


u0 = 4*pi*1e-7;           
Fs=2.5e6;     Fn=0.02154*Fs;     t=0.5/Fs:1/Fs:Fn/Fs-0.5/Fs;
Gx=1/u0;fx=2.5e6/102;   Ax=12e-3/u0;                                
Gy=1/u0;fy=2.5e6/99;  Ay=12e-3/u0;  
Gz=2/u0;fz=2.5e6/96;   Az=12e-3/u0; 

G=[Gx;Gy;Gz];

theta=pi/2;
Hx_t=Ax*cos(2*pi*fx*t+theta);            
Hy_t=Ay*cos(2*pi*fy*t+theta);   
Hz_t=Az*cos(2*pi*fz*t+theta);   
H_t=[Hx_t;Hy_t;Hz_t]; 

dHx_t=-Ax*2*pi*fx*sin(2*pi*fx*t+theta);   
dHy_t=-Ay*2*pi*fy*sin(2*pi*fy*t+theta);
dHz_t=-Az*2*pi*fz*sin(2*pi*fz*t+theta);

dH_t=[dHx_t;dHy_t;dHz_t];                          

dXffp=dH_t./repmat(G,1,size(dH_t,2)); 
dXffp_Amp=sqrt(sum(dXffp.^2));  

Xffp=H_t./repmat(G,1,size(H_t,2));    
xffp=Xffp(1,:);	
yffp=Xffp(2,:);
zffp=Xffp(3,:);



%% grid 
velocity = zeros(19,19,19);
velocity_norm = zeros(19,19,19);

x = -0.019 : 0.002 : 0.019;
y = -0.019 : 0.002 : 0.019;
z = -0.019/2 : 0.001 : 0.019/2;

for i=1:length(x)-1
    for j=1:length(y)-1
        for k=1:length(z)-1
            temp_x_start = x(i);  temp_x_end   = x(i+1);
            temp_y_start = y(j);  temp_y_end   = y(j+1);
            temp_z_start = z(k);  temp_z_end   = z(k+1);

            in_temp_voxel = (xffp >= temp_x_start) & (xffp <= temp_x_end) &...
                            (yffp >= temp_y_start) & (yffp <= temp_y_end) &...
                            (zffp >= temp_z_start) & (zffp <= temp_z_end);

            velocity(i,j,k) = sum(dXffp_Amp(in_temp_voxel)) / sum(in_temp_voxel) / max(dXffp_Amp);
            velocity_norm(i,j,k) = velocity(i,j,k) / max(dXffp_Amp);
        end
    end
end
velocity(ismissing(velocity))=0;


end