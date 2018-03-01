function [ x0,y0 ] = affine_coordinates(x,y,w,h,t,theta)
%AFFINE_COORDINATES Summary of this function goes here
%   Detailed explanation goes here

[x0,y0] = affine_coordinates_inner(x, y, w, h, 1/t, 1, -theta);

end


function [x0,y0] = affine_coordinates_inner(x1, y1, w1, h1, t1, t2, theta)

  Rtheta = theta*pi/180;
  Rot = [cos(Rtheta) -sin(Rtheta);sin(Rtheta) cos(Rtheta)];
  
  % Simulate rotation -> [x1;y1] = Rot*[x1;y1]
  A = Rot*([x1;y1]);
  x1 = A(1,:);
  y1 = A(2,:);
  
  % where is the new origin
  A = Rot*[ [0;h1] [w1;0] [w1;h1] [0;0] ];
  x_ori = min(A(1,:));
  y_ori = min(A(2,:));

  x1 = x1 - x_ori;
  y1 = y1 - y_ori;
  
  % Simulate tilt 
  x0 = x1 * t2;
  y0 = y1 * t1;
  
end
