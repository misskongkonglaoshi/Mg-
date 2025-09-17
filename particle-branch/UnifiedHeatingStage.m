classdef UnifiedHeatingStage < handle
    % UnifiedHeatingStage 统一的加热/熔化阶段求解器 (基于焓方法)
    %%%  修改这部分关于热量的处理  计算热通量时，温度大923K时，加入化学反应的考虑
    properties
        params
        physicalModel
        initialStateTemplate % 用于在事件函数中作为模板
        oxidation_started = false % 标记氧化反应是否开始
        oxidation_threshold_temp % 氧化反应开始的温度阈值
        oxide_break_model % 氧化层破裂模型
    end
    
    methods
        function obj = UnifiedHeatingStage(params, physicalModel)
            obj.params = params;
            obj.physicalModel = physicalModel;
            obj.oxidation_threshold_temp = params.T_reaction_begin; % 使用参数中定义的反应温度阈值
            obj.oxide_break_model = Ox_break(params); % 初始化氧化层破裂模型
        end
        
        % 注释掉原始的solve方法，改用扩展版本
            % function [t_history, state_history] = solve(obj, particleState, t_start, target_temp)
            %     % 求解统一的加热/熔化阶段, 直到达到目标温度
            %     
            %     % 将初始状态保存为模板，供事件函数使用
            %     obj.initialStateTemplate = particleState;
            %     
            %     % 1. 根据初始状态计算初始总焓
            %     H0 = obj.physicalModel.get_enthalpy_from_state(particleState);
            %     
            %     % 2. 定义求解总焓H的ODE
            %     ode_func = @(t, H) obj.ode_system(H, particleState);
            %
            %     % 3. 定义一个用于平滑输出的时间向量，而不是简单的起止点
            %     output_dt = 1e-4; % 输出步长 (s), 可根据需要调整以获得更平滑的曲线
            %     tspan = t_start:output_dt:obj.params.total_time;
            %
            %     % 4. 设置事件函数: 当温度达到目标时停止 (通过焓来判断)
            %     options = odeset('Events', @(t, H) obj.target_temp_event(t, H, target_temp));
            %     
            %     % 5. 调用ODE求解器求解 H(t)
            %     [t_history, H_history] = ode45(ode_func, tspan, H0, options);
            %     
            %     % 6. 将焓的历史记录转换为完整的状态历史记录
            %     num_steps = length(t_history);
            %     if num_steps > 0
            %         state_history(1, num_steps) = ParticleState(); % 预分配
            %         for i = 1:num_steps
            %             state_history(i) = obj.physicalModel.get_state_from_enthalpy(H_history(i), particleState);
            %         end
            %     else
            %         state_history = ParticleState.empty;
            %     end
        % end
        
        function [t_history, state_history] = solve(obj, particleState, t_start, target_temp)
            % 求解统一的加热/熔化阶段, 直到达到目标温度
            % 在温度达到反应阈值后考虑质量变化和氧化层破裂
            
            % 将初始状态保存为模板，供事件函数使用
            obj.initialStateTemplate = particleState;
            obj.oxidation_started = false; % 重置氧化反应标记
            
            % 初始化状态向量 [H, m_mg, m_mgo, m_c, r_c, r_p]
            H0 = obj.physicalModel.get_enthalpy_from_state(particleState);
            y0 = [H0; particleState.m_mg; particleState.m_mgo; particleState.m_c; particleState.r_c; particleState.r_p];
            
            % 定义求解ODE系统的函数
            ode_func = @(t, y) obj.extended_ode_system(t, y, particleState);

            % 定义一个用于平滑输出的时间向量
            output_dt = 1e-4; % 输出步长 (s), 可根据需要调整以获得更平滑的曲线
            tspan = t_start:output_dt:obj.params.total_time;

            % 设置事件函数: 当温度达到目标时停止
            options = odeset('Events', @(t, y) obj.target_temp_event(t, y, target_temp), ...
                           'RelTol', 1e-5, 'AbsTol', 1e-8, ...
                           'OutputFcn', @(t, y, flag) obj.output_function(t, y, flag));
            
            % 调用ODE求解器求解扩展系统
            [t_history, y_history] = ode45(ode_func, tspan, y0, options);
            
            % 将解的历史记录转换为完整的状态历史记录
            num_steps = length(t_history);
            if num_steps > 0
                state_history(1, num_steps) = ParticleState(); % 预分配
                for i = 1:num_steps
                    state_history(i) = obj.convert_solution_to_state(y_history(i,:), particleState);
                end
            else
                state_history = ParticleState.empty;
            end

            % 显示氧化层破裂历史（如果有）
            if ~isempty(obj.oxide_break_model.break_history.time)
                fprintf('统一加热阶段中发生了 %d 次氧化层破裂\n', ...
                    length(obj.oxide_break_model.break_history.time));
                obj.oxide_break_model.visualize_break_history();
            end

            fprintf('统一加热阶段结束\n');
        end

        function status = output_function(obj, t, y, flag)
            % 用于ODE求解器的输出函数，用于记录中间状态
            status = 0; % 继续求解
            
            % 在初始化和求解完成时忽略处理
            if strcmp(flag, 'init') || strcmp(flag, 'done')
                return;
            end
            
            % 每隔一定的步数输出信息
            persistent last_output_time;
            if isempty(last_output_time)
                last_output_time = -inf;
            end
            
            % 只在时间间隔足够大时输出
            current_time = t(end);
            if current_time - last_output_time < 1e-4
                return;
            end
            
            last_output_time = current_time;
            
            % 提取最新状态
            latest_y = y(:, end);
            H_current = latest_y(1);
            m_mg = latest_y(2);
            m_mgo = latest_y(3);
            m_c = latest_y(4);
            r_c = latest_y(5);
            r_p = latest_y(6);
            
            % 构建临时状态
            tempState = obj.build_temp_state(H_current, m_mg, m_mgo, m_c, r_c, r_p, obj.initialStateTemplate);
            
            % 检查是否考虑氧化层破裂
            consider_break = true;
            if isfield(obj.params, 'consider_oxide_break')
                consider_break = obj.params.consider_oxide_break;
            end
            
            if consider_break
                % 每隔一段时间，检查氧化层破裂状态
                dt = 1e-9;  % 时间步长估计
                [is_broken, break_factor, stress_info] = obj.oxide_break_model.calculate_oxide_break(tempState, current_time, dt);
                
                % 更新反应面积因子
                obj.oxide_break_model.update_reaction_area_factor(tempState, break_factor);
            end
            
            % 输出当前温度和氧化状态
            if mod(round(current_time/dt), 50) == 0 % 每隔50个时间步输出一次
                fprintf('t=%.5f s: T=%.1f K, r_p=%.2e m, r_c=%.2e m, 氧化层厚度=%.2e m\n', ...
                    current_time, tempState.T_p, r_p, r_c, tempState.oxide_thickness);
                
                if tempState.T_p >= obj.oxidation_threshold_temp
                    fprintf('  应力状态: 热应力=%.2e Pa, 反应应力=%.2e Pa, 总应力=%.2e Pa, 强度比=%.2f\n', ...
                        stress_info.thermal_stress, stress_info.reaction_stress, ...
                        stress_info.accumulated_stress, stress_info.stress_ratio);
                end
            end

        
            
        end

        function dydt = extended_ode_system(obj, t, y, referenceState)
            % 扩展的ODE系统：同时求解能量和质量变化
            % y = [H, m_mg, m_mgo, m_c, r_c, r_p]

            % 提取当前状态
            H_current = y(1);
            m_mg_current = y(2);
            m_mgo_current = y(3);
            m_c_current = y(4);
            r_c_current = y(5);
            r_p_current = y(6);

            % 从当前值构建临时状态
            tempState = obj.build_temp_state(H_current, m_mg_current, m_mgo_current, m_c_current, r_c_current, r_p_current, referenceState);

            % 检查是否考虑氧化层破裂
            consider_break = isfield(obj.params, 'consider_oxide_break') && obj.params.consider_oxide_break;
            
            if consider_break && tempState.T_p < obj.params.materials.Mg.boiling_point
                dt = 1e-5; % 假定的时间步长
                %fprintf(' 考虑氧化层破裂 计算\n');
                [is_broken, break_factor, ~] = obj.oxide_break_model.calculate_oxide_break(tempState, t, dt);
                obj.oxide_break_model.update_reaction_area_factor(tempState, break_factor);
            else
                %fprintf('  不考虑氧化层破裂\n');
                is_broken = false;
                break_factor = 0;
            end

            % 检查是否达到氧化反应温度阈值
            if tempState.T_p >= obj.oxidation_threshold_temp
                if ~obj.oxidation_started
                    fprintf('  t=%.4f s: 温度达到 %.1f K，开始考虑氧化反应\n', t, tempState.T_p);
                    obj.oxidation_started = true;
                end

                % 计算氧化反应速率 (基于CO2扩散)
                oxidation_rates = obj.calculate_oxidation_rates(tempState);

                % 考虑氧化反应的热流和质量变化
                dHdt = obj.physicalModel.calculate_heat_flux_with_oxidation(tempState, oxidation_rates);
                dmg_dt = oxidation_rates.dmg_dt;
                dmgo_dt = oxidation_rates.dmgo_dt;
                dc_dt = oxidation_rates.dc_dt;
                drc_dt = oxidation_rates.drc_dt;
                drp_dt = oxidation_rates.drp_dt;
                %fprintf('  t=%.4f s: 温度达到 %.1f K,考虑氧化反应带来的质量变化: dmg_dt=%.6e kg/s, dmgo_dt=%.6e kg/s, dc_dt=%.6e kg/s, drc_dt=%.6e m/s, drp_dt=%.6e m/s\n', t, tempState.T_p, dmg_dt, dmgo_dt, dc_dt, drc_dt, drp_dt);
                % 计算热量并存储在临时状态中
                % 1. 对流换热
                h_conv = obj.params.k_gas ;
                q_conv = h_conv / tempState.r_p * (obj.params.ambient_temperature - tempState.T_p);
                tempState.heat_convection = q_conv * (4 * pi * tempState.r_p^2);
                
                % 2. 辐射换热
                q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - tempState.T_p^4);
                tempState.heat_radiation = q_rad * (4 * pi * tempState.r_p^2);
                
                % 3. 反应热
                tempState.heat_reaction_surface = oxidation_rates.reaction_heat ;
                tempState.heat_reaction = tempState.heat_reaction_surface + 0;
                
                % 4. 总热量
                tempState.heat_total = tempState.heat_convection + tempState.heat_radiation + tempState.heat_reaction;
                
                % 考虑破裂后反应速率的增加
                if is_broken
                    %printf('考虑氧化层破裂对表面反应的影响 加入速率增加系数');
                    increase_factor = 1 + break_factor*5 ;  % 破裂后反应速率增加
                    dmg_dt = dmg_dt * increase_factor;
                    dmgo_dt = dmgo_dt * increase_factor;
                    dc_dt = dc_dt * increase_factor;
                    drc_dt = drc_dt * increase_factor;
                    drp_dt = drp_dt * increase_factor;
                    
                    % 额外的反应热
                    extra_heat = oxidation_rates.reaction_heat * (increase_factor - 1);
                    dHdt = dHdt + extra_heat;
                    tempState.heat_reaction = tempState.heat_reaction + extra_heat;
                    tempState.heat_total = tempState.heat_total + extra_heat;
                end
            else
                % 仅考虑传热，无氧化反应
                dHdt = obj.physicalModel.calculate_heat_flux(tempState);
                dmg_dt = 0;
                dmgo_dt = 0;
                dc_dt = 0;
                drc_dt = 0;
                drp_dt = 0;
                
                % 计算热量并存储在临时状态中
                % 1. 对流换热
                h_conv = obj.params.k_gas ;
                q_conv = h_conv / tempState.r_p * (obj.params.ambient_temperature - tempState.T_p);
                tempState.heat_convection = q_conv * (4 * pi * tempState.r_p^2);
                % 2. 辐射换热
                q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - tempState.T_p^4);
                tempState.heat_radiation = q_rad * (4 * pi * tempState.r_p^2);
                
                % 3. 反应热 
                tempState.heat_reaction_surface = 0;

                tempState.heat_reaction = 0;
                
                % 4. 总热量
                tempState.heat_total = tempState.heat_convection + tempState.heat_radiation;
            end

            % 组装导数向量
            dydt = [dHdt; dmg_dt; dmgo_dt; dc_dt; drc_dt; drp_dt];
        end

        function rates = calculate_oxidation_rates(obj, particleState)
            % 计算表面氧化反应的速率
            % 使用与气相燃烧阶段相似的积分公式计算
            % 单区域双边界问题：边界为颗粒核心和氧化层外表面
            
            % 基本参数
            T_p = particleState.T_p; % 颗粒温度
            r_p = particleState.r_p; % 颗粒外半径
            r_c = particleState.r_c; % 颗粒核心半径
            materials = obj.params.materials;
            
            % 1. 设置佩克莱数（CO2扩散的特征量）
            % 考虑孔隙率和弯曲度因子修正扩散系数
            rho_D_gas = obj.params.rho_D_gas; % 环境中的质量扩散系数参数
            porosity = obj.params.material_properties.oxide_porosity;
            tortuosity = obj.params.material_properties.oxide_tortuosity; 
            % porosity = 0.2; % 氧化层孔隙率
            % tortuosity = 2.5; % 弯曲度因子
            
            % 修正后的扩散系数
            rho_D_oxide = rho_D_gas * porosity / tortuosity;
            
            % 2. 估计CO2的质量流率（基于核心表面的反应）
            % 计算反应动力学限制
            k_reaction = obj.calculate_reaction_rate_constant(T_p);
            C_CO2_surf = obj.calculate_CO2_concentration(T_p);
            reaction_area = 4 * pi * r_c^2 * particleState.reaction_area_factor;
            
            % 反应限制的CO2消耗速率 [mol/s]
            dn_CO2_reaction_dt = k_reaction * C_CO2_surf * reaction_area;
            
            % 3. 计算扩散限制 - 使用积分公式
            % 与气相燃烧求解中相似，但简化为单区域问题
            % 介质中扩散通量的解析解 (单区域，固定浓度边界)
            % N = 4*pi*D*r_p*r_c/(r_p-r_c) * (C_surf - C_core)
            geo_factor = 4 * pi * r_p * r_c / (r_p - r_c + eps);
            dn_CO2_diffusion_dt = rho_D_oxide * geo_factor * C_CO2_surf; % 核心表面浓度为0（消耗）
            % 4. 取扩散和反应的最小值作为控制步骤
            dn_CO2_dt = min(dn_CO2_diffusion_dt, dn_CO2_reaction_dt);
            %fprintf('  dn_CO2_diffusion_dt %.6e : dn_CO2_reaction_dt %.6e \n', dn_CO2_diffusion_dt, dn_CO2_reaction_dt);
            % 5. 计算反应物和产物的质量变化
            % 反应: Mg + CO2 → MgO + C
            dn_Mg_dt = -dn_CO2_dt; % 等摩尔反应
            dn_MgO_dt = dn_CO2_dt;
            dn_C_dt = dn_CO2_dt;
            
            % 转换为质量变化率
            rates.dmg_dt = dn_Mg_dt * materials.Mg.molar_mass;
            rates.dmgo_dt = dn_MgO_dt * materials.MgO.molar_mass;
            rates.dc_dt = dn_C_dt * materials.C.molar_mass;
            
            % 6. 计算几何变化
            if T_p < 923
                rho_mg = materials.Mg.density_low;
            else
                rho_mg = materials.Mg.density_high;
            end
            
            rho_mgo = materials.MgO.density;
            rho_c = materials.C.density;
            
            % 核心半径变化
            dV_mg_dt = rates.dmg_dt / rho_mg;
            rates.drc_dt = dV_mg_dt / (4 * pi * r_c^2 + eps);
            
            % 外半径变化
            dV_mgo_dt = rates.dmgo_dt / rho_mgo;
            dV_c_dt = rates.dc_dt / rho_c;
            dV_products_dt = dV_mgo_dt + dV_c_dt;
            rates.drp_dt = dV_products_dt / (4 * pi * r_p^2 + eps);
            
            % 7. 计算反应热
            reaction_heat_value = obj.params.reaction_heat_face;
            rates.reaction_heat = rates.dmg_dt * reaction_heat_value;
            
            % 8. 诊断输出
            %if abs(rates.dmg_dt) > 1e-12
            %    fprintf('  氧化反应: r_c=%.2e m, r_p=%.2e m, dmg/dt=%.2e kg/s\n', ...
            %        r_c, r_p, rates.dmg_dt);
            %    if dn_CO2_diffusion_dt < dn_CO2_reaction_dt
            %        fprintf('  扩散控制: dn_diff=%.2e mol/s, dn_reac=%.2e mol/s\n', ...
            %            dn_CO2_diffusion_dt, dn_CO2_reaction_dt);
            %    else
            %        fprintf('  反应控制: dn_reac=%.2e mol/s, dn_diff=%.2e mol/s\n', ...
            %            dn_CO2_reaction_dt, dn_CO2_diffusion_dt);
            %    end
            %end
        end
        
        function k_reaction = calculate_reaction_rate_constant(obj, T_p)
            % 计算反应速率常数（阿伦尼乌斯方程）
            A_pre = 1.376e4; % 指前因子
            E_a = 1.324e5; % 活化能 [J/mol]
            
            %if isfield(obj.params, 'reaction_pre_exponential')
                A_pre = obj.params.reaction_pre_exponential;
            %else
            %    fprintf('警告: 加热阶段表面异相反应未设置指前因子\n');
            %end
            
            %if isfield(obj.params, 'reaction_activation_energy')
                E_a = obj.params.reaction_activation_energy;
            %else
            %    fprintf('警告: 加热阶段表面异相反应未设置活化能\n');
            %end
            
            R = obj.params.R_u; % 通用气体常数
            k_reaction = A_pre * exp(-E_a / (R * T_p));
        end
        
        function C_CO2 = calculate_CO2_concentration(obj, T_p)
            % 计算CO2浓度
            % 使用理想气体定律和考虑温度效应
            
            % 环境条件
            R_u = obj.params.R_u;
            T_amb = obj.params.ambient_temperature;
            p_amb = obj.params.ambient_pressure;
            
            % 环境中CO2的摩尔分数
            x_CO2 = obj.params.ambient_gas_composition.CO2;
            
            % CO2分压
            p_CO2 = p_amb * x_CO2;
            
            % 摩尔浓度 [mol/m³]
            C_CO2_mol = p_CO2 / (R_u * T_amb);
            
            % 质量浓度 [kg/m³]
            C_CO2 = C_CO2_mol * obj.params.materials.CO2.molar_mass;
            
            % 温度效应修正（氧化层表面浓度受温度影响）
            %T_ref = obj.params.initial_temperature;
            %depletion_factor = 0.2; % 温度影响因子
            
            %if isfield(obj.params, 'CO2_depletion_factor')
            %    depletion_factor = obj.params.CO2_depletion_factor;
            %end
            
            % 浓度随温度衰减
            %C_CO2 = C_CO2 * exp(-depletion_factor * (T_p - T_ref) / T_ref);
            %C_CO2 = max(C_CO2, 0.01 * C_CO2); % 保留最小浓度避免为零
        end
        
        function tempState = build_temp_state(obj, H, m_mg, m_mgo, m_c, r_c, r_p, referenceState)
            % 从当前的焓和质量构建临时状态对象
            tempState = referenceState.copy();
            
            % 更新质量和几何
            tempState.m_mg = m_mg;
            tempState.m_mgo = m_mgo;
            tempState.m_c = m_c;
            tempState.r_c = r_c;
            tempState.r_p = r_p;
            tempState.oxide_thickness = r_p - r_c;
            
            % 从焓反算温度（需要PhysicalModel支持）
            tempState = obj.physicalModel.get_state_from_enthalpy_and_mass(H, tempState);
            
            % 初始化热量字段
            if ~isfield(tempState, 'heat_convection')
                tempState.heat_convection = 0;
                tempState.heat_radiation = 0;
                tempState.heat_reaction = 0;
                tempState.heat_total = 0;
                tempState.heat_reaction_surface = 0 ;
            end
        end
        
        function state = convert_solution_to_state(obj, y_row, referenceState)
            % 将ODE解转换为ParticleState对象
            H = y_row(1);
            m_mg = y_row(2);
            m_mgo = y_row(3);
            m_c = y_row(4);
            r_c = y_row(5);
            r_p = y_row(6);
            
            state = obj.build_temp_state(H, m_mg, m_mgo, m_c, r_c, r_p, referenceState);
        end
        
        function [value, isterminal, direction] = target_temp_event(obj, ~, y, target_temp)
            % 修改后的事件函数，使用扩展状态向量
            H_current = y(1);
            m_mg = y(2);
            m_mgo = y(3);
            m_c = y(4);
            r_c = y(5);
            r_p = y(6);
            
            tempState = obj.build_temp_state(H_current, m_mg, m_mgo, m_c, r_c, r_p, obj.initialStateTemplate);
            
            value = tempState.T_p - target_temp;
            isterminal = 1;
            direction = 1;
        end

        % 保留原始的ODE系统作为备用
        function dHdt = ode_system(obj, H, particleState)
            % 统一加热阶段的ODE系统: dH/dt = Q_total
            
            % a) 从当前的总焓H反算出瞬时状态, 主要是为了得到温度
            tempState = obj.physicalModel.get_state_from_enthalpy(H, particleState);
            
            % b) 根据这个瞬时状态(主要是温度)计算净热流
            dHdt = obj.physicalModel.calculate_heat_flux(tempState);
        end
        
        % 保留原始的事件函数作为备用
        function [value, isterminal, direction] = original_target_temp_event(obj, ~, H, target_temp)
            % 事件函数: 当根据焓算出的温度达到目标温度时触发
            
            % a) 从当前的总焓H反算出瞬时状态 (使用保存的模板)
            tempState = obj.physicalModel.get_state_from_enthalpy(H, obj.initialStateTemplate);
            
            % b) 比较当前温度和目标温度
            value = tempState.T_p - target_temp;
            isterminal = 1;
            direction = 1;
        end
    end
end 