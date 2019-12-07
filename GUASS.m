function [p] = GUASS(delta_r,laser_r)

sigma = laser_r/3;
p = [];
fun = @(x) x.*exp(-x.^2/(2*sigma^2));
parfor i = 1:round(laser_r/delta_r)
    p(i) = integral(fun,(i-1)*delta_r,i*delta_r);
end
p = p/sum(p);
