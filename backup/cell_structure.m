function [cellnum,dphi,phi_divide,sum_cell] = cell_structure( phimax,Rp,r_divide,lambda)
%calculate the num of cells of one plane sector, and dphi of each r
%sum_cell is the sum of cell for r<ri

n=0;
dphi=zeros(r_divide,1);
phi_divide=zeros(r_divide,1);
sum_cell=zeros(r_divide,1);
for i=1:r_divide
    l=i/r_divide*Rp*phimax;
    ratio=mod(l,lambda)/lambda;  %to decide which is closer to l, nfloor*lambda or (nfloor+1)*lambda
    phi_divide(i)=floor(l/lambda);
    if ratio>0.5
        phi_divide(i)=phi_divide(i)+1;
    end
    if phi_divide(i)==0
        phi_divide(i)=1;
    end
    n=n+phi_divide(i);
    dphi(i)=phimax/phi_divide(i);
      if i>1
        sum_cell(i)=sum_cell(i-1)+phi_divide(i);
      else
          sum_cell(i)=phi_divide(i);
      end
    
end
cellnum=n;
end

