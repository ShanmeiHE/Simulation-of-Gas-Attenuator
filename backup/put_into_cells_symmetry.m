function [INDEX] = put_into_cells_symmetryadjust(r_divide,phi_divide,z_divide,dphi,r,phi,z,N_total,Rp,Lp,cellnumxy,sum_cell)
%put N_total mol. into cells, and sign the mol. with INDEX according to the
%cell where it locates
dz=1.0*Lp/z_divide;
dr=1.0*Rp/r_divide;



INDEX=zeros(N_total,1);
index_z = round(z/dz+0.5);
index_r = round(r/dr+0.5);
index_phi = round(phi./dphi(index_r)+0.5);
a = find(index_r~=1);
INDEX(a) = (index_z(a)-1)*cellnumxy+sum_cell(index_r(a)-1)+index_phi(a);
a = find(index_r==1);
INDEX(a) = (index_z(a)-1)*cellnumxy+0+index_phi(a);


end

