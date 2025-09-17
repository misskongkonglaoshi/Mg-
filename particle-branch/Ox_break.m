classdef Ox_break < handle
    % Ox_break 氧化层破裂模型
    % 计算热膨胀和反应产物体积变化导致的接触应力，判断氧化层破裂
    
    properties
        params              % 参数对象
        break_history       % 破裂历史记录
        last_break_time     % 上次破裂的时间
        accumulated_stress  % 累积应力
    end
    
    methods
        function obj = Ox_break(params)
            % 构造函数
            obj.params = params;
            obj.break_history = struct('time', [], 'temperature', [], 'stress', [], 'break_index', []);
            obj.last_break_time = -1;  % 初始化为负值表示未发生过破裂
            obj.accumulated_stress = 0;
        end
        
        function [is_broken, break_factor, stress_info] = calculate_oxide_break(obj, particleState, current_time, dt)
            % 计算氧化层是否发生破裂
            % 输入:
            %   particleState - 当前颗粒状态
            %   current_time - 当前时间点
            %   dt - 时间步长
            % 输出:
            %   is_broken - 是否发生破裂
            %   break_factor - 破裂程度(0-1)，影响有效反应面积
            %   stress_info - 应力相关信息的结构体
            
            % 默认结果
            fprintf(' 开始氧化层破裂 计算\n');
            is_broken = false;
            break_factor = 0;
            stress_info = struct(...
                'thermal_stress', 0, ...
                'reaction_stress', 0, ...
                'total_stress', 0, ...
                'accumulated_stress', 0, ...
                'limit_stress', 1, ...
                'stress_ratio', 0);
            
            % 检查阶段和参数设置
            if ~isfield(obj.params, 'current_stage') || ...
               ~strcmp(obj.params.current_stage, 'heating_and_melting') || ...
               ~obj.params.consider_oxide_break
                return;
            end
            
            % 检查温度是否超过沸点
            if particleState.T_p >= obj.params.materials.Mg.boiling_point
                % 温度超过沸点，不考虑破裂
                return;
            end
            
            % 计算应力
            thermal_stress = obj.calculate_thermal_expansion_stress(particleState);
            reaction_stress = obj.calculate_reaction_volume_stress(particleState);
            
            % 计算净应力（考虑相互抵消）
            if reaction_stress < 0 && thermal_stress > 0 && abs(reaction_stress) >= thermal_stress
                % 如果反应产生的压缩应力能抵消热膨胀引起的张力
                net_stress = 0;
                fprintf('  t=%.6f s: 化学反应应力已抵消热膨胀应力，氧化层稳定\n', current_time);
            else
                % 正常计算总应力
                net_stress = thermal_stress + reaction_stress;
            end
            
            % 累积应力（考虑应力松弛和疲劳效应）
            relaxation_factor = 0.95;
            obj.accumulated_stress = obj.accumulated_stress * relaxation_factor + net_stress * (1 - relaxation_factor);
            
            % 判断是否超过极限应力
            limit_stress = obj.get_oxide_limit_stress(particleState.T_p);
            stress_ratio = obj.accumulated_stress / limit_stress;
            
            % 应力信息结构体
            stress_info = struct(...
                'thermal_stress', thermal_stress, ...
                'reaction_stress', reaction_stress, ...
                'total_stress', thermal_stress + reaction_stress, ...
                'accumulated_stress', obj.accumulated_stress, ...
                'limit_stress', limit_stress, ...
                'stress_ratio', stress_ratio);
            
            % 移除恢复期检查，只保留破裂判断逻辑
            % 破裂判断
            if stress_ratio >= 1.0
                is_broken = true;
                
                % 根据超过极限的程度确定破裂因子
                break_factor = min(1.0, (stress_ratio - 1.0) * 2 + 0.3);
                
                % 记录破裂事件
                obj.record_break_event(current_time, particleState.T_p, obj.accumulated_stress, break_factor);
                
                % 应力完全释放为0
                obj.accumulated_stress = 0;
                
                % 更新上次破裂时间（仅用于记录）
                obj.last_break_time = current_time;
                
                fprintf('  t=%.6f s: 氧化层破裂! 温度=%.1f K, 应力比=%.2f, 破裂因子(break_factor)=%.2f\n', ...
                       current_time, particleState.T_p, stress_ratio, break_factor);
            end
            
            % 定期记录应力历史，无论是否发生破裂
            if mod(round(current_time/1e-6), 10) == 0  % 每10步记录一次
                obj.record_stress_history(current_time, particleState.T_p, obj.accumulated_stress);
            end
        end
        
        function stress = calculate_thermal_expansion_stress(obj, particleState)
            % 计算热膨胀引起的接触应力
            
            % 获取材料参数
            T_p = particleState.T_p;
            T_ref = obj.params.initial_temperature;  % 参考温度
            r_c = particleState.r_c;  % 核心半径
            thickness = particleState.oxide_thickness;  % 氧化层厚度
            
            % 热膨胀系数 (1/K)
            alpha_mg = 25e-6;  % 镁热膨胀系数
            alpha_oxide = 8e-6;  % 氧化镁热膨胀系数
            
            % 弹性模量 (Pa)
            E_mg = 45e9;  % 镁弹性模量
            E_oxide = 250e9;  % 氧化镁弹性模量
            
            % 泊松比
            nu_mg = 0.35;  % 镁泊松比
            nu_oxide = 0.25;  % 氧化镁泊松比
            
            % 使用参数中的值（如果有）
            if isfield(obj.params, 'material_properties')
                props = obj.params.material_properties;
                if isfield(props, 'alpha_mg')
                    alpha_mg = props.alpha_mg;
                end
                if isfield(props, 'alpha_oxide')
                    alpha_oxide = props.alpha_oxide;
                end
                if isfield(props, 'E_mg')
                    E_mg = props.E_mg;
                end
                if isfield(props, 'E_oxide')
                    E_oxide = props.E_oxide;
                end
                if isfield(props, 'nu_mg')
                    nu_mg = props.nu_mg;
                end
                if isfield(props, 'nu_oxide')
                    nu_oxide = props.nu_oxide;
                end
            end
            
            % 温度变化
            delta_T = T_p - T_ref;
            
            % 热膨胀引起的径向位移差异
            delta_r_mg = r_c * alpha_mg * delta_T;
            delta_r_oxide_inner = r_c * alpha_oxide * delta_T;
            
            % 径向位移不匹配量
            mismatch = delta_r_mg - delta_r_oxide_inner;
            
            % 复合圆筒理论计算接触应力
            % 考虑球形几何的简化模型
            k_mg = (1 - nu_mg) / E_mg;
            k_oxide = (1 - nu_oxide) / E_oxide;
            
            % 接触压力估计
            if mismatch > 0 && thickness > 0
                % 正接触压力（膨胀导致挤压）
                contact_pressure = mismatch / (r_c * (k_mg + k_oxide * (1 + r_c/(r_c+thickness))));
                
                % 氧化层的径向应力
                stress = contact_pressure;
            else
                % 如果是收缩，应力为零或者拉应力（但拉应力容易导致开裂，这里简化为零）
                stress = 0;
            end
            
            % 考虑温度对氧化层强度的影响
            if T_p > 700  % 高温会降低材料强度
                temp_factor = 1 + (T_p - 700) / 300;  % 温度校正因子
                stress = stress * temp_factor;
            end
        end
        
        function stress = calculate_reaction_volume_stress(obj, particleState)
            % 计算反应产物体积变化引起的应力
            
            % 如果未发生反应或反应初期，返回零应力
            if particleState.m_c < 1e-12
                stress = 0;
                return;
            end
            
            % 获取材料参数
            r_c = particleState.r_c;  % 核心半径
            thickness = particleState.oxide_thickness;  % 氧化层厚度
            m_mg_reacted = particleState.m_c * (24.3/12);  % 已反应的Mg质量（基于碳的质量和化学计量比）
            
            % 材料密度 (kg/m³)
            rho_mg = obj.params.materials.Mg.density;
            rho_mgo = obj.params.materials.MgO.density;
            rho_c = obj.params.materials.C.density;
            
            % 摩尔质量 (kg/mol)
            mw_mg = obj.params.materials.Mg.molar_mass;
            mw_mgo = obj.params.materials.MgO.molar_mass;
            mw_c = obj.params.materials.C.molar_mass;
            
            % 计算反应前后的体积变化
            V_mg_reacted = m_mg_reacted / rho_mg;  % 反应的Mg体积
            
            % 根据反应方程式：Mg + CO2 → MgO + C
            n_mol = m_mg_reacted / mw_mg;  % 反应的摩尔数
            m_mgo_produced = n_mol * mw_mgo;  % 产生的MgO质量
            m_c_produced = n_mol * mw_c;  % 产生的C质量
            
            V_mgo_produced = m_mgo_produced / rho_mgo;  % 产生的MgO体积
            V_c_produced = m_c_produced / rho_c;  % 产生的C体积
            
            % 反应产物总体积
            V_products = V_mgo_produced + V_c_produced;
            
            % 体积变化比例
            volume_ratio = V_products / V_mg_reacted;
            
            % 弹性模量 (Pa)
            E_oxide = 250e9;  % 氧化镁弹性模量
            
            % 体积应变转换为应力
            % 简化模型：假设体积变化在球壳内均匀分布，并转换为径向应力
            if volume_ratio > 1.0
                % 膨胀情况
                volumetric_strain = volume_ratio - 1.0;
                % 径向应变约等于体积应变的1/3
                radial_strain = volumetric_strain / 3;
                % 应力估计（胡克定律简化版）
                stress = E_oxide * radial_strain;
            else
                % 收缩情况，产生压应力（负值）
                volumetric_strain = 1.0 - volume_ratio;
                radial_strain = volumetric_strain / 3;
                stress = -E_oxide * radial_strain * 0.5;
            end
            
            % 应力衰减因子 - 考虑氧化层厚度的影响
            if thickness > 0
                thickness_factor = 1.0 / (1.0 + 5.0 * thickness / r_c);
                stress = stress * thickness_factor;
            end
        end
        
        function limit_stress = get_oxide_limit_stress(obj, temperature)
            % 获取氧化层在当前温度下的极限应力
            
            % 基础极限应力值（室温）
            base_limit_stress = 300e6;  % 300 MPa
            
            % 使用参数中的值（如果有）
            if isfield(obj.params, 'material_properties') && ...
               isfield(obj.params.material_properties, 'oxide_limit_stress')
                base_limit_stress = obj.params.material_properties.oxide_limit_stress;
            end
            
            % 温度对强度的影响（高温强度降低）
            T_ref = 300;  % 参考温度(K)
            if temperature <= T_ref
                temp_factor = 1.0;
            else
                % 温度升高导致强度下降
                temp_factor = exp(-(temperature - T_ref) / 500);
            end
            
            limit_stress = base_limit_stress * temp_factor;
        end
        
        function update_reaction_area_factor(obj, particleState, break_factor)
            % 根据破裂情况更新颗粒的有效反应面积因子
           
            % 初始反应面积因子（未破裂情况下由氧化层厚度决定）
            base_factor = 1.0;
            if particleState.oxide_thickness > particleState.initial_oxide_thickness
                base_factor = particleState.initial_oxide_thickness / particleState.oxide_thickness;
            end
            
            % 考虑破裂的影响 - 增加有效反应面积
            if break_factor > 0
                % 破裂使得内部金属更容易与外部气体接触
                effective_factor = base_factor + break_factor * (1.0 - base_factor);
                
                % 随机变异，模拟破裂的随机性
                randomness = 0.9 + 0.2 * rand();
                effective_factor = effective_factor * randomness;
                
                % 更新颗粒状态中的反应面积因子
                particleState.reaction_area_factor = min(1.0, effective_factor);
            else
                % 未破裂，使用基础因子
                particleState.reaction_area_factor = base_factor;
            end
            fprintf('更新有效反应面积因子\n');
        end
        
        function record_break_event(obj, time, temperature, stress, break_factor)
            % 检查阶段和温度
            if ~isfield(obj.params, 'current_stage') || ...
               ~strcmp(obj.params.current_stage, 'heating_and_melting') || ...
               temperature >= obj.params.materials.Mg.boiling_point
                return; % 不符合条件，不记录
            end
            
            % 检查温度是否为负或不合理值
            if temperature <= 0 || ~isfinite(temperature)
                fprintf('警告: 记录到不合理的温度值: %f K, 使用默认值\n', temperature);
                temperature = 1000; % 使用一个合理的默认温度
            end
            
            % 记录破裂事件
            idx = length(obj.break_history.time) + 1;
            obj.break_history.time(idx) = time;
            obj.break_history.temperature(idx) = temperature;
            obj.break_history.stress(idx) = stress;
            obj.break_history.break_index(idx) = break_factor;
            
            % 最多记录100个事件
            if idx > 100
                obj.break_history.time = obj.break_history.time(end-99:end);
                obj.break_history.temperature = obj.break_history.temperature(end-99:end);
                obj.break_history.stress = obj.break_history.stress(end-99:end);
                obj.break_history.break_index = obj.break_history.break_index(end-99:end);
            end
        end
        
        function visualize_break_history(obj)
            % 可视化破裂历史
            if isempty(obj.break_history.time)
                fprintf('没有记录到氧化层破裂事件\n');
                return;
            end
            
            figure('Name', '氧化层破裂历史', 'Position', [100, 100, 800, 600]);
            
            % 时间-温度曲线，标记破裂点
            subplot(2, 1, 1);
            plot(obj.break_history.time, obj.break_history.temperature, 'b-', 'LineWidth', 1.5);
            hold on;
            scatter(obj.break_history.time, obj.break_history.temperature, 50, ...
                    obj.break_history.break_index, 'filled');
            colorbar;
            xlabel('时间 (s)');
            ylabel('温度 (K)');
            title('氧化层破裂时刻的温度');
            grid on;
            
            % 总应力与极限应力对比
            subplot(2, 1, 2);
            plot(obj.break_history.time, obj.break_history.stress, 'r-', 'LineWidth', 2, 'DisplayName', '累积应力');
            hold on;
            
            % 获取每个时间点对应的极限应力
            limit_stresses = zeros(size(obj.break_history.time));
            for i = 1:length(obj.break_history.time)
                limit_stresses(i) = obj.get_oxide_limit_stress(obj.break_history.temperature(i));
            end
            
            plot(obj.break_history.time, limit_stresses, 'k--', 'LineWidth', 1.5, 'DisplayName', '极限应力');
            
            % 标记破裂点
            break_times = obj.break_history.time;
            break_stresses = obj.break_history.stress;
            scatter(break_times, break_stresses, 80, 'r', 'filled', 'DisplayName', '破裂点');
            
            xlabel('时间 (s)');
            ylabel('应力 (Pa)');
            title('应力演变与极限应力对比');
            legend('Location', 'best');
            grid on;
            
            sgtitle('氧化层破裂历史分析');
        end
    end
end
