% 运行Mg颗粒燃烧仿真的主脚本

% 清理工作空间
clear;
clc;

% 初始化参数
params = Parameters();

% 创建物理模型
model = PhysicalModel(params);

% 运行仿真并获取完整结果
fprintf('开始仿真 (求解器: %s)...\n', params.solver_config.type);
[t, y, rates_history] = model.solve();
fprintf('仿真完成！\n');

% 初始化结果存储
num_steps = length(t);
results = struct();
results.time = t;
results.temperature = y(1,:);
results.mass_mg = y(2,:);
results.mass_mgo = y(3,:);
results.mass_c = y(4,:);
results.melted_fraction = y(5,:);
results.diameter = 2 * y(6,:);  

% 从历史记录中提取速率
results.heat_rate = [rates_history.heat_rate];
% 将反应速率从 kg/s 转换为 mmol/s 用于绘图
reaction_rate_kg_s = [rates_history.reaction_rate];
results.reaction_rate = (reaction_rate_kg_s / params.materials.Mg.molar_mass) * 1000; % (kg/s -> mol/s -> mmol/s)

% 绘制结果
fprintf('正在生成图表...\n');
figure('Name', ['Mg颗粒燃烧仿真结果 (求解器: ' params.solver_config.type ')']);

% 温度变化
subplot(2, 2, 1);
plot(results.time * 1000, results.temperature, 'LineWidth', 1.5);
xlabel('时间 (ms)');
ylabel('温度 (K)');
title('颗粒温度演化');
grid on;

% 熔化分数
subplot(2, 2, 2);
plot(results.time * 1000, results.melted_fraction, 'LineWidth', 1.5);
xlabel('时间 (ms)');
ylabel('熔化分数');
title('熔化进度');
grid on;

% 热量变化率
subplot(2, 2, 3);
plot(results.time * 1000, results.heat_rate / 1000, 'LineWidth', 1.5);
xlabel('时间 (ms)');
ylabel('热量变化率 (kW)');
title('净热流');
grid on;

% 反应速率
subplot(2, 2, 4);
plot(results.time * 1000, results.reaction_rate, 'LineWidth', 1.5);
xlabel('时间 (ms)');
ylabel('反应速率 (mmol/s)');
title('Mg消耗速率');
grid on;

% 调整图形布局
set(gcf, 'Position', [100, 100, 1000, 800]);
fprintf('完成！\n'); 