function [VC] = calculate_VC( Rp,Lp,r_divide,phi_divide,z_divide,dphi,cellnumxy,sum_cell,cellnum )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
VC=zeros(cellnum,1);
dz=1.0*Lp/z_divide;
dr=1.0*Rp/r_divide;
for i=1:r_divide
    for j=1:phi_divide(i)
        if i~=1
            cellorder=sum_cell(i-1)+j;
        else
            cellorder=j;
        end
        VC(cellorder)=0.5*((i*dr)^2-((i-1)*dr)^2)*dphi(i)*dz;
    end
end
for i=1:cellnumxy
    for j=1:z_divide-1
    VC(i+j*cellnumxy)=VC(i);
    end
end
end

