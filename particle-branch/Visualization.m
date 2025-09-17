classdef Visualization
    % Visualization 结果可视化类
    
    methods (Static)
        function plotResults(results)
            % 绘制仿真结果
            %
            % 输入:
            %   results: 仿真结果结构体, 包含火焰面数据
            
            % 检查结果是否为空
            if isempty(results) || ~isfield(results, 'time') || isempty(results.time)
                fprintf('结果为空, 无法绘图。\n');
                return;
            end
            
            figure('Name', 'Mg颗粒燃烧仿真结果', 'Position', [100, 100, 1800, 500]);

            % --- 子图1: 火焰面和颗粒温度 ---
            subplot(1, 3, 1);
            hold on;
            plot(results.time * 1000, results.temperature, 'b-', 'LineWidth', 1.5, 'DisplayName', '颗粒温度');
            
            % 仅在有火焰数据时绘制
            if isfield(results, 'flame_temperature')
                valid_flame_indices = ~isnan(results.flame_temperature);
                if any(valid_flame_indices)
                    plot(results.time(valid_flame_indices) * 1000, results.flame_temperature(valid_flame_indices), 'r--', 'LineWidth', 1.5, 'DisplayName', '火焰面温度');
                end
            end
            
            hold off;
            xlabel('时间 (ms)');
            ylabel('温度 (K)');
            title('温度演化');
            legend('Location', 'best');
            grid on;

            % --- 子图2: 半径变化 ---
            subplot(1, 3, 2);
            hold on;
            
            % 计算核心半径历史
            % 注意: 此处硬编码了Mg的密度, 因为它不在results结构中
            density_mg = 1740; % kg/m^3
            % 确保质量不为负值
            safe_mass_mg = max(results.mass_mg, 0);
            r_c_history = (3 .* safe_mass_mg ./ (4 * pi * density_mg)).^(1/3);
            
            % 确保所有数据都是实数
            plot(results.time * 1000, real(results.radius * 1e6), 'k-', 'LineWidth', 1.5, 'DisplayName', '颗粒外径 (r_p)');
            plot(results.time * 1000, real(r_c_history * 1e6), 'g-.', 'LineWidth', 1.5, 'DisplayName', 'Mg核心半径 (r_c)');
            
            if isfield(results, 'flame_radius')
                valid_flame_indices = ~isnan(results.flame_radius);
                if any(valid_flame_indices)
                    plot(results.time(valid_flame_indices) * 1000, real(results.flame_radius(valid_flame_indices) * 1e6), 'r--', 'LineWidth', 1.5, 'DisplayName', '火焰面半径');
                end
            end
            
            hold off;
            xlabel('时间 (ms)');
            ylabel('半径 (µm)');
            title('几何尺寸演化');
            legend('Location', 'best');
            grid on;

            % --- 子图3: 质量变化 ---
            subplot(1, 3, 3);
            hold on ;
            
            % 使用正确的系数且确保数据是实数
            plot(results.time * 1000, real(results.mass_mg * 1e8), 'm--', 'LineWidth',1.5, 'DisplayName','mg质量');
            plot(results.time * 1000, real(results.mass_mgo * 1e8), 'k-', 'LineWidth',1.5, 'DisplayName','mgo质量');
            plot(results.time * 1000, real(results.mass_c * 1e8), 'g-', 'LineWidth',1.5, 'DisplayName','c质量');


            hold off ;
            xlabel('时间 (ms)');
            ylabel('质量');
            title('质量变化');
            % ylim([-0.05, 1.05]);
            grid on;

            % 调整布局
            set(gcf, 'Color', 'w');
        end
        
        function saveResults(results, filename)
            % 保存仿真结果到MAT文件
            %
            % 输入:
            %   results: 仿真结果结构体
            %   filename: 保存文件名
            
            save(filename, 'results');
            fprintf('结果已保存到文件: %s\n', filename);
        end
        
        function exportToCSV(results, filename)
            % 导出仿真结果到CSV文件
            %
            % 输入:
            %   results: 仿真结果结构体
            %   filename: 保存文件名
            
            % 准备数据表
            data_table = table(results.time, results.temperature, results.melted_fraction, ...
                'VariableNames', {'Time', 'Temperature', 'MeltedFraction'});
            
            % 添加阶段信息
            data_table.Stage = categorical(results.stage);
            
            % 添加氧化层厚度（如果存在）
            if isfield(results, 'oxide_thickness')
                data_table.OxideThickness = results.oxide_thickness;
            end
            
            % 添加质量信息（如果存在）
            if isfield(results, 'mass_mg')
                data_table.MassMg = results.mass_mg;
            end
            
            if isfield(results, 'mass_mgo')
                data_table.MassMgO = results.mass_mgo;
            end
            
            if isfield(results, 'mass_c')
                data_table.MassC = results.mass_c;
            end
            
            % 写入CSV文件
            writetable(data_table, filename);
            fprintf('结果已导出到CSV文件: %s\n', filename);
        end
    end
end 