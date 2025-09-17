% Mg颗粒燃烧仿真主程序
% 分阶段模拟Mg颗粒在CO2环境中的加热、熔融和反应过程

% 清理工作空间
clear;
clc;

% 初始化参数
params = Parameters();
% 创建物理模型
model = ParticleModel(params);
 
% 创建阶段管理器
manager = StageManager(params);

% 运行仿真
fprintf('开始仿真...\n');
fprintf('初始氧化层厚度: %.2e m\n', params.initial_oxide_thickness);
results = manager.solve();
fprintf('仿真完成！\n');

% 调用外部可视化函数
%Visualization.plotResults(results);

% 保存结果
%if ~exist('results', 'dir')
%    mkdir('results');
%end
%save('results/simulation_results.mat', 'results');
% 使用Visualization类的专用方法导出CSV，而不是直接转换结构体
%Visualization.exportToCSV(results, 'results/simulation_results.csv');
% writetable(struct2table(results), 'results/simulation_results.csv');


% 打印仿真信息
fprintf('\n仿真信息:\n');
fprintf('总仿真时间: %.6f 秒\n', results.time(end));
fprintf('最终温度: %.2f K\n', results.temperature(end));
fprintf('最终熔融分数: %.4f\n', results.melted_fraction(end));
fprintf('最终氧化层厚度: %.4e m (%.2f µm)\n', results.oxide_thickness(end), results.oxide_thickness(end)*1e6);

% 如果存在质量数据，打印质量信息
if isfield(results, 'mass_mg') && isfield(results, 'mass_mgo')
    fprintf('初始Mg粒径: %.6f um\n', results.radius(1)*1000000*2);
    fprintf('初始Mg质量: %.4e kg\n', results.mass_mg(1));
    fprintf('最终Mg质量: %.4e kg\n', results.mass_mg(end));
    fprintf('Mg消耗量: %.2f%%\n', (1 - results.mass_mg(end)/results.mass_mg(1)) * 100);
    
    fprintf('初始MgO质量: %.8e kg\n', results.mass_mgo(1));
    fprintf('最终MgO质量: %.4e kg\n', results.mass_mgo(end));
    fprintf('MgO增长量: %.2f%%\n', (results.mass_mgo(end)/results.mass_mgo(1) - 1) * 100);

    fprintf('火焰温度: %.6f K\n', results.flame_temperature(1));
    fprintf('初始Mg质量: %.4e kg\n', results.mass_mg(1));
end

% 打印各阶段时间
unique_stages = unique(results.stage);
fprintf('\n各阶段时间:\n');
for i = 1:length(unique_stages)
    stage_indices = strcmp(results.stage, unique_stages{i});
    stage_time = results.time(stage_indices);
    if ~isempty(stage_time)
        fprintf('%s: %.6f - %.6f ms (持续 %.6f ms)\n', ...
            unique_stages{i}, stage_time(1)*1000, stage_time(end)*1000, ...
            stage_time(end)*1000 - stage_time(1)*1000);
    end
end 