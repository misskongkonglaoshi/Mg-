classdef VaporizationODE < handle
    % VaporizationODE 气相燃烧阶段ODE求解器类
    % 从StageManager中提取的气相燃烧阶段ODE相关代码，进行独立封装
    
    properties
        params               % 参数对象
        physicalModel        % 物理计算模型
        vaporization_rate_info_cache % 缓存BVP结果
        
        % 可视化相关属性
        visualize_realtime = false % 是否开启实时可视化
        realtime_fig         % 实时可视化图形句柄
        realtime_axes        % 实时可视化子图句柄
        realtime_history     % 实时可视化历史数据
        last_visual_update_time = 0 % 上次可视化更新的时间
    end
    
    methods
        function obj = VaporizationODE(params, physicalModel)
            % 构造函数：初始化气相燃烧阶段ODE求解器
            obj.params = params;
            obj.physicalModel = physicalModel;
            obj.vaporization_rate_info_cache = {}; % 初始化缓存
            
            % 设置可视化参数
            obj.visualize_realtime = params.visualize_realtime;
        end
        
        function [t_stage, y_stage] = solve(obj, pState, current_time)
            % 求解气相燃烧阶段ODE系统
            % 输入:
            %   pState - 初始颗粒状态
            %   current_time - 当前时间
            % 输出:
            %   t_stage - 时间历史
            %   y_stage - 状态向量历史
            
            % 将颗粒温度锁定在沸点
            T_boil = obj.params.materials.Mg.ignition_temp;
            pState.T_p = T_boil;
            
            % 创建VaporizationStage对象用于计算反应速率
            vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
            
            % 设置ODE的初始状态向量 [m_mg, m_mgo, m_c, r_c, r_p, T_p]
            y0 = [pState.m_mg, pState.m_mgo, pState.m_c, pState.r_c, pState.r_p, pState.T_p];
            
            % 根据参数选择使用固定步长
            fprintf('使用固定步长求解器, 步长: %.2e s\n', obj.params.vaporization_fixed_timestep);
            t_span = [current_time, obj.params.total_time];
            [t_stage, y_stage] = obj.solve_with_fixed_timestep(...
                @(t,y) obj.vaporization_ode(t, y, vaporization_solver, T_boil), ...
                t_span, y0, obj.params.vaporization_fixed_timestep);
        end
        
        function [t, y] = solve_with_fixed_timestep(obj, fun, t_span, y0, step_size)
            % 使用固定步长显式欧拉法求解ODE
            % 输入:
            %   fun: ODE右侧函数句柄 @(t,y)
            %   t_span: 求解时间区间 [t_start, t_end]
            %   y0: 初始状态向量
            %   step_size: 固定时间步长
            
            t_start = t_span(1);
            t_end = t_span(2);
            
            % 生成时间步数组
            t = (t_start:step_size:t_end)';
            if t(end) < t_end
                t = [t; t_end];
            end
            
            % 初始化结果数组
            y = zeros(length(t), length(y0));
            y(1,:) = y0;
            
            % 提取T_boil用于输出函数调用
            T_boil = obj.params.materials.Mg.ignition_temp;
            
            % 初始化调用
            obj.vaporization_output_fcn(t(1), y(1,:)', 'init');
            
            % 初始化代数方程组解的缓存变量
            algebraic_solution_cache = [];
            vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
            
            % 时间推进
            for i = 1:length(t)-1
                dt = t(i+1) - t(i);
                
                % 调用事件函数检查是否应该停止
                [value, isterminal, ~] = obj.vaporization_events(t(i), y(i,:)');
                if isterminal && value <= 0
                    fprintf('在t=%.4f s检测到终止事件，提前结束计算\n', t(i));
                    % 截断时间和状态数组
                    t = t(1:i);
                    y = y(1:i,:);
                    break;
                end
                
                % 调用输出函数进行状态可视化和记录
                obj.vaporization_output_fcn(t(i), y(i,:)', '');
                
                % 创建一个临时的ParticleState，用于缓存信息提取
                tempState = obj.build_temp_state_from_vector(y(i,:)', T_boil);
                
                % 使用缓存的代数解直接调用求解器
                rate_info = vaporization_solver.solve_reaction_rates(tempState, algebraic_solution_cache);
                
                % 更新缓存供下一时间步使用
                if isfield(rate_info, 'algebraic_solution') && ~isempty(rate_info.algebraic_solution)
                    algebraic_solution_cache = rate_info.algebraic_solution;
                end

                if ~rate_info.success
                    fprintf('在t=%.4f s代数方程求解失败，提前结束计算\n', t(i));
                    % 截断时间和状态数组到当前位置
                    t = t(1:i);
                    y = y(1:i,:);
                    break; % 跳出时间推进循环
                end
                
                % 使用rate_info计算右侧函数
                dydt = obj.compute_dydt_from_rate_info(t(i), y(i,:)', rate_info, T_boil);
                
                % 更新下一个状态
                y(i+1,:) = y(i,:) + dt * dydt';
                
                % 保存到缓存，供可视化使用
                if isempty(obj.vaporization_rate_info_cache)
                    obj.vaporization_rate_info_cache = {};
                end
                obj.vaporization_rate_info_cache{end+1} = rate_info;
                
                % 打印进度
                if mod(i, 100) == 0
                    if ~isempty(algebraic_solution_cache)
                        cache_status = '已使用缓存';
                    else
                        cache_status = '无缓存';
                    end
                    fprintf('固定步长求解进度: %.2f%% (t=%.4fs, 缓存状态: %s)\n', ...
                        100*(t(i)-t_start)/(t_end-t_start), t(i), cache_status);
                end
            end
            
            % 最后一步输出函数调用
            obj.vaporization_output_fcn(t(end), y(end,:)', 'done');
        end

        function dydt = vaporization_ode(obj, t, y, solver, T_p_const)
            % 气相燃烧阶段ODE系统右侧函数
            % y(1)=m_mg, y(2)=m_mgo, y(3)=m_c, y(4)=r_c, y(5)=r_p, y(6)=T_p
            
            % 为防止数值问题，增加保护性判断
            if y(1) < 1e-15 || y(4) < 1e-9 || y(5) < 1e-9
                dydt = zeros(6,1);
                return;
            end

            % 1. 从当前状态向量y构建一个临时的ParticleState对象
            tempState = obj.build_temp_state_from_vector(y, T_p_const);
                
            % 2. 获取质量变化率
            try
                % 调用求解器获取反应速率
                rate_info = solver.solve_reaction_rates(tempState, []);
                
                % 改进: 检查求解是否成功
                if ~rate_info.success
                    error('扩散燃烧部分求解不收敛：在t=%.4f s处无法获得有效解', t);
                end
                
                dmdt_mg  = rate_info.dmdt_mg;
                dmdt_mgo = rate_info.dmdt_mgo;
                dmdt_c   = rate_info.dmdt_c;
                
                % 3. 计算表面温度变化率
                % 获取热量信息
                Q_total = rate_info.heat_ox;  % 总热量流入（W）
                
                % 计算氧化层有效热容
                C_oxide = obj.calculate_oxide_heat_capacity(y);
                
                % 表面温度变化率 (dT/dt = Q/C)
                if C_oxide > 1e-10
                    dTp_dt = Q_total / C_oxide;
                else
                    dTp_dt = 0;  % 避免除零
                end
                
                % 将热量信息添加到rate_info
                if ~isfield(rate_info, 'heat_convection')
                    % 默认热量信息
                    rate_info.heat_convection = 0;
                    rate_info.heat_radiation = 0;
                    rate_info.heat_reaction = 0;
                    rate_info.heat_reaction_surface = 0;
                    rate_info.heat_total = Q_total;
                end
                
            catch ME
                fprintf('t=%.4f s 时反应速率求解失败: %s\n', t, ME.message);
                % 重新抛出错误，这将中止ODE求解器
                rethrow(ME);
            end

            % 3. 计算半径变化率
            % 使用相同的密度参数确保一致性
            rho_mg = obj.params.materials.Mg.density_high;
            rho_mgo = obj.params.materials.MgO.density;
            rho_c = obj.params.materials.C.density;

            % 核心半径变化率 (基于Mg消耗)
            dVdt_mg = dmdt_mg / rho_mg;
            drc_dt = dVdt_mg / (4 * pi * y(4)^2 + eps);

            % 外部颗粒半径变化率 (基于产物沉积)
            dVdt_mgo = dmdt_mgo / rho_mgo;
            dVdt_c = dmdt_c / rho_c;
            dVdt_product_deposition = dVdt_mgo + dVdt_c;
            drp_dt = dVdt_product_deposition / (4 * pi * y(5)^2 + eps);

            % 4. 组合成完整的导数向量
            dydt = [dmdt_mg; dmdt_mgo; dmdt_c; drc_dt; drp_dt; dTp_dt];
        end
        
        function dydt = compute_dydt_from_rate_info(obj, t, y, rate_info, T_p_const)
            % 基于已有rate_info计算ODE右侧函数值
            
            % 初始化返回值
            dydt = zeros(size(y));
            
            % 设置质量变化率
            dydt(1) = rate_info.dmdt_mg;     % Mg质量变化率
            dydt(2) = rate_info.dmdt_mgo;    % MgO质量变化率
            dydt(3) = rate_info.dmdt_c;      % C质量变化率
           
  
            % 核心半径变化率
            rho_mg = obj.params.materials.Mg.density_low;
            dVdt_mg = dydt(1) / rho_mg;
            dydt(4) = dVdt_mg / (4 * pi * y(4)^2 + eps); % r_c变化率
            
            % 外部半径变化率
            rho_mgo = obj.params.materials.MgO.density;
            rho_c = obj.params.materials.C.density;
            dVdt_mgo = dydt(2) / rho_mgo;
            dVdt_c = dydt(3) / rho_c;
            dVdt_product_deposition = dVdt_mgo + dVdt_c;
            dydt(5) = dVdt_product_deposition / (4 * pi * y(5)^2 + eps); % r_p变化率
            
            % 温度变化率计算
            if isfield(rate_info, 'heat_ox')
                C_oxide = obj.calculate_oxide_heat_capacity(y);
                if C_oxide > 1e-10
                    dydt(6) = rate_info.heat_ox / C_oxide; % 温度变化率
                else
                    dydt(6) = 0;
                end
            else
                dydt(6) = 0;
            end
        end
        
        function C_oxide = calculate_oxide_heat_capacity(obj, y)
            % 计算氧化层热容量
            % 提取状态向量值
            m_mgo = y(2);
            m_c = y(3);
            T_p = y(6);
            
            % 计算氧化层质量分数
            m_total = m_mgo + m_c;
            if m_total > 0
                mass_frac_mgo = m_mgo / m_total;
                mass_frac_c = m_c / m_total;
            else
                mass_frac_mgo = 1;
                mass_frac_c = 0;
            end
            
            % 获取各组分比热
            cp_mgo_molar = obj.physicalModel.calc_cp_solid_mgo(T_p);  % J/(mol·K)
            cp_c_molar = obj.physicalModel.calc_cp_solid_c(T_p);      % J/(mol·K)
            
            % 转换为质量比热
            cp_mgo_mass = cp_mgo_molar / obj.params.materials.MgO.molar_mass;  % J/(kg·K)
            cp_c_mass = cp_c_molar / obj.params.materials.C.molar_mass;        % J/(kg·K)
            
            % 计算混合比热
            cp_mix = mass_frac_mgo * cp_mgo_mass + mass_frac_c * cp_c_mass;
            
            % 计算热容量
            porosity = obj.params.material_properties.oxide_porosity;
            m_oxide_eff = m_total; 
            C_oxide = m_oxide_eff * cp_mix;  % [J/K]
        end
        
        function [value, isterminal, direction] = vaporization_events(obj, t, y)
            % ODE事件函数: 当Mg质量降至初始值的阈值以下时停止
            
            % 获取当前Mg质量
            current_mg_mass = y(1);
            rho_mg = obj.params.materials.Mg.density_low;
            % 获取初始Mg质量
            initial_mg_mass = (obj.params.initial_diameter/2)^3*4/3 * pi * rho_mg;
            
            % 计算质量比例
            mass_ratio = current_mg_mass / initial_mg_mass;
            
            % 当质量比例低于阈值时触发事件
            threshold = obj.params.combustionfinish_ratio;
            value = mass_ratio - threshold;
            
            isterminal = 1;    % 触发后停止积分
            direction = -1;    % 质量比例降低穿过阈值时触发
        end
        
        function status = vaporization_output_fcn(obj, t, y, flag)
            % 用于ODE求解器的输出函数，用于记录中间状态
            status = 0; % 默认继续积分
            
            switch flag
                case 'init'
                    % 初始化实时可视化
                    if obj.visualize_realtime
                        obj.init_realtime_visualization();
                    end
                    
                case ''
                    % 不再重复计算，只记录当前状态
                    
                    % 添加对温度值的监控，确保温度在物理合理范围内
                    if length(y) >= 6
                        T_p_current = y(6);  % 当前表面温度
                        
                        % 检查温度是否在合理范围内
                        if isnan(T_p_current) || T_p_current < 0 || T_p_current > 5000
                            fprintf('警告: 在t=%.4f ms处检测到异常表面温度: %.2f K\n', t*1000, T_p_current);
                        end
                        
                        % 每50次迭代输出一次温度信息
                        persistent output_counter;
                        if isempty(output_counter)
                            output_counter = 0;
                        end
                        
                        output_counter = output_counter + 1;
                        if mod(output_counter, 50) == 0
                            fprintf('t=%.4f ms: 核心温度=%.2f K, 表面温度=%.2f K, 温差=%.2f K\n', t*1000, obj.params.materials.Mg.ignition_temp, T_p_current, T_p_current-obj.params.materials.Mg.ignition_temp);
                        end
                    end
                    
                    % 实时可视化颗粒状态
                    if obj.visualize_realtime
                        % 限制可视化更新频率，避免过度绘图导致性能问题
                        if t(end) - obj.last_visual_update_time > 1e-4
                            obj.visualize_realtime_state(t(end), y(:,end));
                            obj.last_visual_update_time = t(end);
                        end
                    end
                    
                case 'done'
                    fprintf('积分完成，共缓存%d个速率信息\n', length(obj.vaporization_rate_info_cache));
                    % 在求解完成时保存可视化结果
                    if obj.visualize_realtime && ishandle(obj.realtime_fig)
                        % 获取颗粒初始直径(μm)，而不是当前直径
                        d0_um = obj.params.initial_diameter * 2 * 1e6;
                        filename = sprintf('%.1fμm_气相燃烧结果.png', d0_um);
                        saveas(obj.realtime_fig, filename);
                        fprintf('实时可视化结果已保存为: %s\n', filename);
                    end
            end
        end
        
        function tempState = build_temp_state_from_vector(obj, y_vec, T_p_const)
            % 从状态向量构建临时ParticleState对象
            tempState = ParticleState(obj.params);
            tempState.m_mg  = y_vec(1);
            tempState.m_mgo = y_vec(2);
            tempState.m_c   = y_vec(3);
            tempState.r_c   = y_vec(4);
            tempState.r_p   = y_vec(5);
            
            % 使用状态变量中的温度作为氧化层表面温度
            if length(y_vec) >= 6
                tempState.T_p = y_vec(6);
            else
                tempState.T_p = T_p_const;  % 默认使用沸点温度
            end
            
            tempState.oxide_thickness = tempState.r_p - tempState.r_c;
            tempState.T_c = T_p_const;  % 金属核心温度锁定在沸点
            tempState.melted_fraction = 1.0;
        end
        
        function [state_history, flame_info_history] = convert_state_vector_to_state_history(obj, t_stage, y_stage, T_p_const)
            % 从ODE结果和缓存构建状态历史
            num_steps = size(y_stage, 1);
            if num_steps == 0
                state_history = ParticleState.empty;
                flame_info_history = [];
                return;
            end
            
            % 预分配
            state_history(1, num_steps) = ParticleState();
            flame_info_history(1, num_steps) = struct('r_f', NaN, 'T_f', NaN);
            
            % 检查缓存是否为空
            rate_info_cache = obj.vaporization_rate_info_cache;
            if isempty(rate_info_cache) || length(rate_info_cache) == 0
                fprintf('警告: 缓存为空，创建默认缓存项\n');
                % 创建一个默认的rate_info结构
                default_rate_info = struct('success', false, 'r_flame', NaN, 'T_flame', NaN);
                rate_info_cache = cell(1, num_steps);
                for i = 1:num_steps
                    rate_info_cache{i} = default_rate_info;
                end
            % 检查缓存大小与步数是否匹配
            elseif num_steps > length(rate_info_cache)
                % 扩展缓存以匹配步数
                fprintf('补充缺失的缓存项(%d → %d)\n', length(rate_info_cache), num_steps);
                if length(rate_info_cache) > 0
                    last_valid = rate_info_cache{end};
                else
                    % 如果缓存为空但长度不为0
                    last_valid = struct('success', false, 'r_flame', NaN, 'T_flame', NaN);
                end
                for i = length(rate_info_cache)+1:num_steps
                    rate_info_cache{i} = last_valid;
                end
            elseif num_steps < length(rate_info_cache)
                % 截断缓存
                rate_info_cache = rate_info_cache(1:num_steps);
            end
            
            % 填充状态历史和火焰信息
            for i = 1:num_steps
                % 从状态向量填充ParticleState对象
                currentState = y_stage(i, :);
                tempState = obj.build_temp_state_from_vector(currentState, T_p_const);
                
                % 确保温度值在合理范围内
                if tempState.T_p < 0 || isnan(tempState.T_p)
                    tempState.T_p = T_p_const; % 使用沸点温度作为后备值
                end
                
                state_history(i) = tempState;
                
                % 从缓存中读取火焰信息
                rate_info = rate_info_cache{i};
                
                % 确保关键字段存在
                if isfield(rate_info, 'r_flame')
                    flame_info_history(i).r_f = rate_info.r_flame;
                elseif isfield(rate_info, 'r_f')
                    flame_info_history(i).r_f = rate_info.r_f;
                else
                    flame_info_history(i).r_f = NaN;
                end
                
                if isfield(rate_info, 'T_flame')
                    flame_info_history(i).T_f = rate_info.T_flame;
                elseif isfield(rate_info, 'T_f')
                    flame_info_history(i).T_f = rate_info.T_f;
                else
                    flame_info_history(i).T_f = NaN;
                end
                
                % 复制热量数据
                if isfield(rate_info, 'heat_convection')
                    flame_info_history(i).heat_convection = rate_info.heat_convection;
                    flame_info_history(i).heat_radiation = rate_info.heat_radiation;
                    flame_info_history(i).heat_reaction = rate_info.heat_reaction;
                    flame_info_history(i).heat_reaction_surface = rate_info.heat_reaction_surface;
                    flame_info_history(i).heat_total = rate_info.heat_total;
                end
            end
        end
        
        % 可视化相关方法
        function init_realtime_visualization(obj)
            % 创建或获取实时可视化图形
            if isempty(obj.realtime_fig) || ~ishandle(obj.realtime_fig)
                obj.realtime_fig = figure('Name', '颗粒燃烧实时监控', 'Position', [100, 100, 1000, 800]);
                
                % 创建子图
                obj.realtime_axes.temp = subplot(2, 2, 1);
                title('温度变化');
                xlabel('时间 (ms)');
                ylabel('温度 (K)');
                hold(obj.realtime_axes.temp, 'on');
                grid(obj.realtime_axes.temp, 'on');
                
                obj.realtime_axes.radius = subplot(2, 2, 2);
                title('半径变化');
                xlabel('时间 (ms)');
                ylabel('半径 (μm)');
                hold(obj.realtime_axes.radius, 'on');
                grid(obj.realtime_axes.radius, 'on');
                
                obj.realtime_axes.mass = subplot(2, 2, 3);
                title('质量变化');
                xlabel('时间 (ms)');
                ylabel('质量 (μg)');
                hold(obj.realtime_axes.mass, 'on');
                grid(obj.realtime_axes.mass, 'on');
                
                % 将热量图改为质量变化率图
                obj.realtime_axes.mass_rate = subplot(2, 2, 4);
                title('质量变化率');
                xlabel('时间 (ms)');
                % 初始设置左侧y轴
                yyaxis(obj.realtime_axes.mass_rate, 'left');
                ylabel('Mg变化率 (μg/ms)');
                % 初始设置右侧y轴
                yyaxis(obj.realtime_axes.mass_rate, 'right');
                ylabel('MgO与C变化率 (μg/ms)');
                hold(obj.realtime_axes.mass_rate, 'on');
                grid(obj.realtime_axes.mass_rate, 'on');
                
                % 初始化历史数据
                obj.realtime_history = struct(...
                    'time', [], ...
                    'T_p', [], ...
                    'T_c', [], ...
                    'r_p', [], ...
                    'r_c', [], ...
                    'm_mg', [], ...
                    'm_mgo', [], ...
                    'm_c', [], ...
                    'dmdt_mg', [], ...    % 添加Mg质量变化率
                    'dmdt_mgo', [], ...   % 添加MgO质量变化率
                    'dmdt_c', [], ...     % 添加C质量变化率
                    'heat_convection', [], ...
                    'heat_radiation', [], ...
                    'heat_reaction', [], ...
                    'heat_total', []);
            end
        end

        function visualize_realtime_state(obj, t, y)
            % 从状态向量y创建临时颗粒状态
            T_boil = obj.params.materials.Mg.ignition_temp;
            tempState = obj.build_temp_state_from_vector(y, T_boil);
            
            % 如果有缓存的速率信息，获取最新的热量数据和质量变化率数据
            rate_info = struct('heat_convection', 0, 'heat_radiation', 0, 'heat_reaction', 0, 'heat_total', 0, ...
                              'dmdt_mg', 0, 'dmdt_mgo', 0, 'dmdt_c', 0);
            if ~isempty(obj.vaporization_rate_info_cache)
                latest_rate_info = obj.vaporization_rate_info_cache{end};
                if isfield(latest_rate_info, 'heat_convection')
                    rate_info.heat_convection = latest_rate_info.heat_convection;
                    rate_info.heat_radiation = latest_rate_info.heat_radiation;
                    rate_info.heat_reaction = latest_rate_info.heat_reaction;
                    rate_info.heat_total = latest_rate_info.heat_total;
                end
                if isfield(latest_rate_info, 'dmdt_mg')
                    rate_info.dmdt_mg = latest_rate_info.dmdt_mg;
                    rate_info.dmdt_mgo = latest_rate_info.dmdt_mgo;
                    rate_info.dmdt_c = latest_rate_info.dmdt_c;
                end
            end
            
            % 更新历史数据
            obj.realtime_history.time(end+1) = t*1000; % 转换为ms
            obj.realtime_history.T_p(end+1) = tempState.T_p;
            obj.realtime_history.T_c(end+1) = tempState.T_c;
            obj.realtime_history.r_p(end+1) = tempState.r_p*1e6; % 转换为μm
            obj.realtime_history.r_c(end+1) = tempState.r_c*1e6; % 转换为μm
            obj.realtime_history.m_mg(end+1) = tempState.m_mg*1e6; % 转换为μg
            obj.realtime_history.m_mgo(end+1) = tempState.m_mgo*1e6; % 转换为μg
            obj.realtime_history.m_c(end+1) = tempState.m_c*1e6; % 转换为μg
            obj.realtime_history.dmdt_mg(end+1) = rate_info.dmdt_mg*1e9; % 转换为μg/ms
            obj.realtime_history.dmdt_mgo(end+1) = rate_info.dmdt_mgo*1e9; % 转换为μg/ms
            obj.realtime_history.dmdt_c(end+1) = rate_info.dmdt_c*1e9; % 转换为μg/ms
            obj.realtime_history.heat_convection(end+1) = rate_info.heat_convection;
            obj.realtime_history.heat_radiation(end+1) = rate_info.heat_radiation;
            obj.realtime_history.heat_reaction(end+1) = rate_info.heat_reaction;
            obj.realtime_history.heat_total(end+1) = rate_info.heat_total;
            
            % 更新温度图
            cla(obj.realtime_axes.temp);
            plot(obj.realtime_axes.temp, obj.realtime_history.time, obj.realtime_history.T_p, 'r-', 'LineWidth', 2);
            hold(obj.realtime_axes.temp, 'on');
            plot(obj.realtime_axes.temp, obj.realtime_history.time, obj.realtime_history.T_c, 'b--', 'LineWidth', 1.5);
            legend(obj.realtime_axes.temp, '表面温度', '核心温度', 'Location', 'best');
            title(obj.realtime_axes.temp, sprintf('温度变化 (t=%.2f ms)', t*1000));
            
            % 更新半径图
            cla(obj.realtime_axes.radius);
            plot(obj.realtime_axes.radius, obj.realtime_history.time, obj.realtime_history.r_p, 'r-', 'LineWidth', 2);
            hold(obj.realtime_axes.radius, 'on');
            plot(obj.realtime_axes.radius, obj.realtime_history.time, obj.realtime_history.r_c, 'b--', 'LineWidth', 1.5);
            plot(obj.realtime_axes.radius, obj.realtime_history.time, obj.realtime_history.r_p - obj.realtime_history.r_c, 'g:', 'LineWidth', 1.5);
            legend(obj.realtime_axes.radius, '颗粒半径', '核心半径', '氧化层厚度', 'Location', 'best');
            title(obj.realtime_axes.radius, sprintf('半径变化 (t=%.2f ms)', t*1000));
            
            % 更新质量图
            cla(obj.realtime_axes.mass);
            plot(obj.realtime_axes.mass, obj.realtime_history.time, obj.realtime_history.m_mg, 'r-', 'LineWidth', 2);
            hold(obj.realtime_axes.mass, 'on');
            plot(obj.realtime_axes.mass, obj.realtime_history.time, obj.realtime_history.m_mgo, 'b--', 'LineWidth', 1.5);
            plot(obj.realtime_axes.mass, obj.realtime_history.time, obj.realtime_history.m_c, 'g:', 'LineWidth', 1.5);
            total_mass = obj.realtime_history.m_mg + obj.realtime_history.m_mgo + obj.realtime_history.m_c;
            plot(obj.realtime_axes.mass, obj.realtime_history.time, total_mass, 'k-.', 'LineWidth', 1);
            legend(obj.realtime_axes.mass, 'Mg质量', 'MgO质量', 'C质量', '总质量', 'Location', 'best');
            title(obj.realtime_axes.mass, sprintf('质量变化 (t=%.2f ms)', t*1000));
            
            % 更新质量变化率图（使用双y轴）
            cla(obj.realtime_axes.mass_rate);
            
            % 左侧Y轴显示Mg变化率
            yyaxis(obj.realtime_axes.mass_rate, 'left');
            plot(obj.realtime_axes.mass_rate, obj.realtime_history.time, obj.realtime_history.dmdt_mg, 'r-', 'LineWidth', 2);
            ylabel('Mg变化率 (μg/ms)', 'Color', 'r');
            set(obj.realtime_axes.mass_rate, 'YColor', 'r');
            
            % 右侧Y轴显示MgO和C变化率
            yyaxis(obj.realtime_axes.mass_rate, 'right');
            plot(obj.realtime_axes.mass_rate, obj.realtime_history.time, obj.realtime_history.dmdt_mgo, 'b-', 'LineWidth', 1.5);
            hold(obj.realtime_axes.mass_rate, 'on');
            plot(obj.realtime_axes.mass_rate, obj.realtime_history.time, obj.realtime_history.dmdt_c, 'g-', 'LineWidth', 1.5);
            ylabel('MgO与C变化率 (μg/ms)', 'Color', 'b');
            set(obj.realtime_axes.mass_rate, 'YColor', 'b');
            
            legend(obj.realtime_axes.mass_rate, 'Mg变化率', 'MgO变化率', 'C变化率', 'Location', 'best');
            title(obj.realtime_axes.mass_rate, sprintf('质量变化率 (t=%.2f ms)', t*1000));
            grid(obj.realtime_axes.mass_rate, 'on');
            
            % 更新图形显示
            drawnow;
        end
    end
end