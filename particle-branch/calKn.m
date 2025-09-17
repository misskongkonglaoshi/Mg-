% Mg 颗粒在不同粒径、压强、温度条件下的 Kn 数变化
% 考虑气体可能为多组分的情况
clear;
clc;

% 定义气体常数和阿伏伽德罗常数
R_u = 8.314; % 通用气体常数，单位：J/(mol·K)
N_A = 6.022e23; % 阿伏伽德罗常数，单位：mol^-1

% 定义气体组分的平均直径（可以根据实际气体组分进行调整）
% 单位：米
D_a = [3.3e-10,3.3e-10,3.3e-10]; % 示例：气体组分的平均直径

% 调用固定温度下的计算函数
T_fixed = 2000;  % 固定温度值（K）
calculateKnWithFixedTemperature(T_fixed, D_a, R_u, N_A);

% 调用固定压强下的计算函数
P_fixed = 1e5;   % 固定压强值（Pa）
calculateKnWithFixedPressure(P_fixed, D_a, R_u, N_A);

function calculateKnWithFixedTemperature(T_fixed, D_a, R_u, N_A)
    % 在固定温度下计算不同压强和粒径的 Kn 数
    % 输入参数：
    % T_fixed: 固定温度值（K）
    % D_a: 气体分子平均直径数组（m）
    % R_u: 通用气体常数
    % N_A: 阿伏伽德罗常数
    
    % 定义压强和粒径范围
    P_a = 0.5e5:1000:2e5; % 压强范围（Pa）
    D_p = 2e-5:1e-7:1e-4; % 颗粒直径范围（m）
    
    % 初始化 Kn 数矩阵
    Kn = zeros(length(D_p), length(P_a));
    
    % 计算 Kn 数
    for i = 1:length(D_p)
        for j = 1:length(P_a)
            % 计算多组分气体的平均直径
            D_a_avg = mean(D_a);
            % 计算 Kn 数
            Kn(i,j) = (R_u * T_fixed) / (sqrt(2) * pi * D_a_avg^2 * N_A * P_a(j) * D_p(i));
        end
    end
    
    % 创建网格数据
    [DP, PA] = meshgrid(D_p * 1e6, P_a);
    
    % 绘制三维曲面图
    figure('Name', ['固定温度 T = ' num2str(T_fixed) ' K']);
    subplot(2,1,1);
    surf(DP, PA, Kn');
    xlabel('颗粒直径（μm）');
    ylabel('环境压强（Pa）');
    zlabel('Kn 数');
    title(['Kn 数在温度 ' num2str(T_fixed) ' K 下的变化']);
    colorbar;
    grid on;
    
    % 绘制等高线图
    subplot(2,1,2);
    contourf(DP, PA, Kn');
    xlabel('颗粒直径（μm）');
    ylabel('环境压强（Pa）');
    title(['Kn 数在温度 ' num2str(T_fixed) ' K 下的等高线图']);
    colorbar;
    grid on;
    
    % 输出特征值
    fprintf('\n固定温度 T = %d K 下的 Kn 数特征值：\n', T_fixed);
    fprintf('最小值：%.4f\n', min(Kn(:)));
    fprintf('最大值：%.4f\n', max(Kn(:)));
    fprintf('平均值：%.4f\n', mean(Kn(:)));
end

function calculateKnWithFixedPressure(P_fixed, D_a, R_u, N_A)
    % 在固定压强下计算不同温度和粒径的 Kn 数
    % 输入参数：
    % P_fixed: 固定压强值（Pa）
    % D_a: 气体分子平均直径数组（m）
    % R_u: 通用气体常数
    % N_A: 阿伏伽德罗常数
    
    % 定义温度和粒径范围
    T_a = 0:200:3000; % 温度范围（K）
    D_p = 2e-5:1e-7:1e-4; % 颗粒直径范围（m）
    
    
    % 初始化 Kn 数矩阵
    Kn = zeros(length(D_p), length(T_a));
    
    % 计算 Kn 数
    for i = 1:length(D_p)
        for j = 1:length(T_a)
            % 计算多组分气体的平均直径
            D_a_avg = mean(D_a);
            % 计算 Kn 数
            Kn(i,j) = (R_u * T_a(j)) / (sqrt(2) * pi * D_a_avg^2 * N_A * P_fixed * D_p(i));
        end
    end
    
    % 创建网格数据
    [DP, TA] = meshgrid(D_p * 1e6, T_a);
    
    % 绘制三维曲面图
    figure('Name', ['固定压强 P = ' num2str(P_fixed) ' Pa']);
    subplot(2,1,1);
    surf(DP, TA, Kn');
    xlabel('颗粒直径（μm）');
    ylabel('温度（K）');
    zlabel('Kn 数');
    title(['Kn 数在压强 ' num2str(P_fixed) ' Pa 下的变化']);
    colorbar;
    grid on;
    
    % 绘制等高线图
    subplot(2,1,2);
    contourf(DP, TA, Kn');
    xlabel('颗粒直径（μm）');
    ylabel('温度（K）');
    title(['Kn 数在压强 ' num2str(P_fixed) ' Pa 下的等高线图']);
    colorbar;
    grid on;
    
    % 输出特征值
    fprintf('\n固定压强 P = %.0f Pa 下的 Kn 数特征值：\n', P_fixed);
    fprintf('最小值：%.4f\n', min(Kn(:)));
    fprintf('最大值：%.4f\n', max(Kn(:)));
    fprintf('平均值：%.4f\n', mean(Kn(:)));
end