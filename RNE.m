%% Recursive Newton-Euler Algorithm for dynamic model of serial robots with modified DH parameters.
%% Take iiwa7 as example.
function tau_list = RNE(q,qd,qdd,G)
%% 输入参数：
% q：广义关节坐标，此处为关节转角，1×7矩阵，每一行向量对应一组关节转角，单位：rad；
% qd： 广义关节坐标一阶导数，此处为关节角速度，1×7矩阵，每一行向量对应一组关节角速度，单位：rad/s；
% qdd： 广义关节坐标二阶导数，此处为关节角加速度，1×7矩阵，每一行向量对应一组关节角加速度，单位：rad/s^2；
% G：重力项，当G = 1,有重力影响，当G = 0,无重力影响；
% note:：三个输入参数的长度需保持一致。

%% 输出参数：
% tau_list ：关节力矩，1×7矩阵，每一行向量对应一组关节力矩，单位：Nm

% 判断输入是否符合规则
rows = size(q,1);
if rows ~= size(qd,1) || rows ~= size(qdd,1)
    error("输入参数长度不一致");
end

% 参数初始化
% DH_list：机器人DH参数，4×7矩阵
%          alpha      a     d     theta
dh_list = [0          0    0.34     0;
           -pi/2      0    0        0;
           pi/2       0    0.4      0;
           -pi/2      0    0        0;
           pi/2       0    0.4      0;
           -pi/2      0    0        0;
           pi/2       0    0.1266   0;
           0          0    0        0];
       
% mass_list: 连杆的质量，1×7矩阵，单位：kg       
mass_list = [2702.4 2725.8 3175.01 2725.80 1693.85 1836.74 269.17]/1000;

% mass_center_list：连杆质心在连杆坐标系下的位置，3×7矩阵，单位：m                             
%                   x         y           z
mass_center_list = [0        -34.73       -69.48;
                    0        -67.33       34.41;
                    -0.03    29.56        -89.00;
                    0.03     -67.33       -34.41;
                    0        -21.39       -140.03;
                    0        -2.12        0.29;
                    0.01      0           -25.22]/1000;
                
% inertia_tensor_list：连杆关于质心坐标系的惯性张量，质心坐标系与连杆坐标系方位一致，7个3×3矩阵，单位kg*m^2
%         I                =      Ixx            -Ixy           -Ixz
%                                 -Ixy           Iyy            -Iyz
%                                 -Iyz           -Iyz           Izz
inertia_tensor_list(:,:,1) = [17085955.96       29.13           520.90;
                              29.13             16299848.40     3041655.32;
                              520.90            3041655.32      6028717.92]/1e9;
                          
inertia_tensor_list(:,:,2) = [17049081.74       -379.65         -251.02;
                              -379.65           6095392.93      2836636.56;
                              -251.02           2836636.56      16245722.28]/1e9;
                          
inertia_tensor_list(:,:,3) = [25077403.67       1373.03         4792.34;
                              1373.03           23806776.16     -4872887.52;
                              4792.34           -4872887.52      7607337.29]/1e9;
                          
inertia_tensor_list(:,:,4) = [17049008.27       -2836.83        605.34;
                              -2836.83           6095457.66      2836684.49;
                              605.34            2836684.49      16245750.37]/1e9; 
                          
inertia_tensor_list(:,:,5) = [10079214.34       -74.51          17.98;
                              -74.51             8702598.01      3090329.67;
                              17.98             3090329.67      4469563.66]/1e9;
                          
inertia_tensor_list(:,:,6) = [5094485.23        -96.90           -67.21;
                              -96.90            3542620.51      -249580.52;
                              -67.21            -249580.52      4899002.53]/1e9;
                          
inertia_tensor_list(:,:,7) = [198764.02         0.88            -27.71;
                              0.88              195312.12       0.20;
                              -27.71            0.20            322516.91]/1e9;
                          
% f_external：施加在末端连杆的外力和外力矩
f_external = zeros(2,3);

number_of_links = 7;
z = [0,0,1]';  % 关节轴
%%判断是否施加重力
if G == 1
    g = 9.81;     % 重力加速度，单位m/s^2
else
    g = 0;
end

% 位姿变换矩阵参数设置
for i = 1:number_of_links+1
    dh = dh_list(i,:);
    alpha(i) = dh(1);
    a(i) = dh(2);
    d(i) = dh(3);
    theta(i) = dh(4);
    if i == number_of_links+1
        q(i) = 0;
    end
    T(:,:,i) = [cos(q(i)),            -sin(q(i)),           0,           a(i);
            sin(q(i))*cos(alpha(i)), cos(q(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i);
            sin(q(i))*sin(alpha(i)), cos(q(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i);
            0,                     0,                     0,          1];
    T = T(:,:,i);
    % 提取旋转矩阵并求逆
    R(:,:,i) = inv(T(1:3,1:3));
    P(:,:,i) = T(1:3,4:4);
end

for k = 1: rows
    % 外推 --->
    for i = 0:number_of_links-1
        if i == 0
            wi = [0,0,0]';      % 初始角速度为0
            dwi = [0,0,0]';     % 初始角加速度为0
            dvi = [0, 0, g]';   % 初始加速度，根据坐标系0对重力加速度方向进行设置
        else
            wi = w(:,i);
            dwi = dw(:,i);
            dvi = dv(:,i);
        end
        w(:,:,i+1) = R(:,:,i+1)*wi + qd(k,i+1)*z;
        dw(:,:,i+1) = R(:,:,i+1)*dwi + cross(R(:,:,i+1)*wi,qd(k,i+1)*z) + qdd(k,i+1)*z;
        dv(:,:,i+1) = R(:,:,i+1)*(cross(dwi,P(:,:,i+1)) + cross(wi,cross(wi,P(:,:,i+1))) + dvi);
        dvc(:,:,i+1) = cross(dw(:,:,i+1),mass_center_list(i+1,:)')...
                        + cross(w(:,:,i+1),cross(w(:,:,i+1),mass_center_list(i+1,:)'))...
                        + dv(:,:,i+1);
        F(:,:,i+1) = mass_list(i+1)*dvc(:,:,i+1);
        N(:,:,i+1) = inertia_tensor_list(:,:,i+1)*dw(:,:,i+1) + cross(w(:,:,i+1),inertia_tensor_list(:,:,i+1)*w(:,:,i+1));
    end
    % 内推 <---
    for i = number_of_links:-1:1
        if i == number_of_links
            f(:,:,i+1) = f_external(1,:)';
            n(:,:,i+1) = f_external(2,:)';
        end
        f(:,:,i) = R(:,:,i+1)\f(:,:,i+1) + F(:,:,i);
        n(:,:,i) = N(:,:,i) + R(:,:,i+1)\n(:,:,i+1) + cross(mass_center_list(i,:)',F(:,:,i))...
                    + cross(P(:,:,i+1),R(:,:,i+1)\f(:,:,i+1));
        tau_list (k,i) = dot(n(:,:,i),z);
    end
end
