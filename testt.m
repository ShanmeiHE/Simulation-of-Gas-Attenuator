clear
r=0.01;
z=0.01;
e=1.6*10^(-19);
T=300;
k_B=1.38*10^(-23);
m=4.65*10^(-26);
me=9.1096*10^(-31);
epsilon=8.854187817*10^(-12);
mu=4*pi*10^(-7);
delta = 1*10^(-4);
gridx = -r:delta:r;
gridy = gridx;
gridz = 0:delta:z;
delta_t=0.5*10^(-13);
sigma=10^(-15);
n=10^5;
W=1;
c=3*10^8;
phi = 0.5;
delta_r = delta;
delta_z = delta;

H_theta_gridr = delta_r/2:delta_r:r-delta_r/2;
H_theta_gridz = delta_z/2:delta_z:z-delta_z/2;

H_theta_lenr = length(H_theta_gridr);
H_theta_lenz = length(H_theta_gridz);
H_theta_len = length(H_theta_gridr)*length(H_theta_gridz);

[H_theta_R,H_theta_Z] = meshgrid(H_theta_gridr,H_theta_gridz);
H_theta_R = H_theta_R';
H_theta_Z = H_theta_Z';



E_r_gridr = delta_r/2:delta_r:r-delta_r/2;
E_r_gridz = 0:delta_z:z;

E_r_lenr = length(E_r_gridr);
E_r_lenz = length(E_r_gridz);
E_r_len = length(E_r_gridr)*length(E_r_gridz);

[E_r_R,E_r_Z] = meshgrid(E_r_gridr,E_r_gridz);
E_r_R = E_r_R';
E_r_Z = E_r_Z';



E_z_gridr = 0:delta_r:r;
E_z_gridz = delta_z/2:delta_z:z-delta_z/2;

E_z_lenr = length(E_z_gridr);
E_z_lenz = length(E_z_gridz);
E_z_len = length(E_z_gridr)*length(E_z_gridz);

[E_z_R,E_z_Z] = meshgrid(E_z_gridr,E_z_gridz);
E_z_R = E_z_R';
E_z_Z = E_z_Z';

%sqrt(electron_velocity(:,1).^2+electron_velocity(:,2).^2+electron_velocity(:,3).^2);



%{
electron_position(1,1) = 0;
electron_position(1,2) = 0;
electron_position(1,3) = z/2;
electron_velocity(1,1) = 0;
electron_velocity(1,2) = 0;
electron_velocity(1,3) = 10^6;
electron_position(2,1) = 0;
electron_position(2,2) = 0;
electron_position(2,3) = z/2;
electron_velocity(2,1) = 0;
electron_velocity(2,2) = 0;
electron_velocity(2,3) = 10^6;
%}


gridE_r = zeros(E_r_lenr,E_r_lenz);
gridj_r = zeros(E_r_lenr,E_r_lenz);
gridE_z = zeros(E_z_lenr,E_z_lenz);
gridj_z = zeros(E_z_lenr,E_z_lenz);
gridH_theta = zeros(H_theta_lenr,H_theta_lenz,1);
gridj_r(1,3) = 10^3;
gridj_r(2,9) = 10^3;
delta_t_repeat = delta_t;
repeat = 1;
for time=1:10000
    delta_t_repeat = delta_t/repeat;
	gridH_theta = gridH_theta - delta_t_repeat/2/(mu*delta_r)*( gridE_z(2:E_z_lenr,:)-gridE_z(1:(E_z_lenr-1),:) ) + delta_t_repeat/2/(mu*delta_z)*( gridE_r(:,2:E_r_lenz)-gridE_r(:,1:(E_r_lenz-1)) );
    for timee = 1:repeat
        gridH_theta = gridH_theta + delta_t_repeat/(mu*delta_r)*( gridE_z(2:E_z_lenr,:)-gridE_z(1:(E_z_lenr-1),:) ) - delta_t_repeat/(mu*delta_z)*( gridE_r(:,2:E_r_lenz)-gridE_r(:,1:(E_r_lenz-1)) );
        %sum(sum(delta_t_repeat/(mu*delta_r)*( gridE_z(2:E_z_lenr,:)-gridE_z(1:(E_z_lenr-1),:) ) - delta_t_repeat/(mu*delta_z)*( gridE_r(:,2:E_r_lenz)-gridE_r(:,1:(E_r_lenz-1)) )))
        gridE_r(:,2:(E_r_lenz-1)) = gridE_r(:,2:(E_r_lenz-1)) - delta_t_repeat/epsilon*gridj_r(:,2:(E_r_lenz-1)) - delta_t_repeat/(epsilon*delta_z)*( gridH_theta(:,2:H_theta_lenz)-gridH_theta(:,1:(H_theta_lenz-1)) );
        %sum(sum(delta_t_repeat/(epsilon*delta_z)*( gridH_theta(:,2:H_theta_lenz)-gridH_theta(:,1:(H_theta_lenz-1)) )))
        gridE_z(2:(E_z_lenr-1),:) = gridE_z(2:(E_z_lenr-1),:) - delta_t_repeat/epsilon*gridj_z(2:(E_z_lenr-1),:) + delta_t_repeat/epsilon*(1./(2*E_z_R(2:(E_z_lenr-1),:)+1/delta_r)).*gridH_theta(2:H_theta_lenr,:) + delta_t_repeat/epsilon*(1./(2*E_z_R(2:(E_z_lenr-1),:)-1/delta_r)).*gridH_theta(1:(H_theta_lenr-1),:) ;
        %sum(sum(delta_t_repeat/epsilon*(1./(2*E_z_R(2:(E_z_lenr-1),:)+1/delta_r)).*gridH_theta(2:H_theta_lenr,:) + delta_t_repeat/epsilon*(1./(2*E_z_R(2:(E_z_lenr-1),:)-1/delta_r)).*gridH_theta(1:(H_theta_lenr-1),:) ))
        gridE_z(1,:) = gridE_z(1,:) + 4*delta_t_repeat/(epsilon*delta_r)*gridH_theta(1,:);
    end
    gridH_theta = gridH_theta + delta_t_repeat/2/(mu*delta_r)*( gridE_z(2:E_z_lenr,:)-gridE_z(1:(E_z_lenr-1),:) ) - delta_t_repeat/2/(mu*delta_z)*( gridE_r(:,2:E_r_lenz)-gridE_r(:,1:(E_r_lenz-1)) );
        
    
    energy_of_electric_field = trapz(E_r_gridz,trapz(E_r_gridr,gridE_r.^2.*E_r_R))+trapz(E_z_gridz,trapz(E_z_gridr,gridE_z.^2.*E_z_R))
    energy_of_magnetic_field = trapz(H_theta_gridz,trapz(H_theta_gridr,gridH_theta.^2.*H_theta_R))    

    %energy_of_electric_field = 1/2*epsilon*energy_of_electric_field
    %energy_of_magnetic_field = 1/2*mu*energy_of_magnetic_field
    energy_total = energy_of_electric_field+energy_of_magnetic_field

    quiver(H_theta_R,H_theta_Z,gridE_r(:,2:101),gridE_z(1:100,:));
    fmat = getframe;
    gridj_r = zeros(E_r_lenr,E_r_lenz);
    gridj_z = zeros(E_z_lenr,E_z_lenz);

end
%{
        
        curl_r = curlB_r;
        curl_z = curlB_z;
        %zeros(lenr,lenz,lentheta)
        theta = 0:0.1:6.28;
        lentheta = length(theta);
        [Z3,R3,theta3] = meshgrid(gridr,gridz,theta);
        
        gridB3 = kron(reshape(curl_r,[len,1]),ones(1,lentheta));
        gridB_z = kron(reshape(curl_z,[len,1]),ones(1,lentheta));
        gridB3 = reshape(gridB3,[lenr,lenz,lentheta]);
        gridB_z = reshape(gridB_z,[lenr,lenz,lentheta]);
        sum(sum(abs(gridB3(:,:,1)-gridB_theta)));
        X = R3.*cos(theta3);
        Y = R3.*sin(theta3);
        gridB_x = -gridB3.*R3.*sin(theta3);
        gridB_y = gridB3.*R3.*cos(theta3);
        
        quiver3(X,Y,Z3,gridB_x,gridB_y,gridB_z,0.1);
        %}




