function R = euler2rot(euler)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

roll = euler(1);
pitch = euler(2);
yaw = euler(3);

Rx = [1,         0,          0;
      0, cos(roll), -sin(roll);
      0, sin(roll),  cos(roll)];
Ry = [ cos(pitch), 0, sin(pitch);
                0, 1,          0;
      -sin(pitch), 0, cos(pitch)];
Rz = [cos(yaw), -sin(yaw), 0;
      sin(yaw),  cos(yaw), 0;
             0,         0, 1];
R = Rz*Ry*Rx;
end

