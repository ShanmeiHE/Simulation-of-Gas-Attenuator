function remainder =mod_y (a,z,m)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
q=floor(m/a);r=mod(m,a);
s=a*mod(z,q)-r*floor(z/q);
	if(s>=0) 
        remainder=s;
	else 
        remainder=s+m;
    end

end

