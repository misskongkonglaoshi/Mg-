classdef StageManager < handle
    % StageManager 阶段管理器
    % 管理不同阶段的切换和求解 (动态调度器)
    % 
    % 气相燃烧阶段使用VaporizationODE类进行ODE求解，该类负责气相燃烧的所有计算过程
    
    properties
        model         % 颗粒模型 (ParticleModel)
        params        % 参数对象 (Parameters)
        physicalModel % 物理计算模型 (PhysicalModel)
        thermo_reader % 热力学数据读取器 (ThermoReader)
        
        % 结果存储
        results       % 包含所有历史记录的结构体
        vaporization_rate_info_cache % 气化阶段ODE求解中缓存的BVP结果，从VaporizationODE类获取
        heat_history  % 热量历史记录
        
        % 添加实时可视化相关属性
        realtime_fig         % 实时可视化图形句柄
        realtime_axes        % 实时可视化子图句柄
        realtime_history     % 实时可视化历史数据
        visualize_realtime = false % 是否开启实时可视化
        last_visual_update_time = 0 % 上次可视化更新的时间
    end
    
    methods
        function obj = StageManager(params)
            % 构造函数：初始化阶段管理器
            if nargin > 0
                obj.params = params;
                
                % 初始化依赖的模型
                obj.thermo_reader = ThermoReader('thermo.dat', obj.params);
                obj.physicalModel = PhysicalModel(obj.params, obj.thermo_reader);
                obj.model = ParticleModel(obj.params);
                
                % 初始化结果存储
                obj.results.time = [];
                obj.results.particleStateHistory = [];
                obj.results.stageHistory = {};
                obj.results.flame_radius = [];      % 新增: 初始化火焰半径历史
                obj.results.flame_temperature = []; % 新增: 初始化火焰温度历史
                obj.results.oxide_break_history = struct('time', [], 'temperature', [], 'stress', [], 'break_index', []);
                
                % 初始化热量历史记录
                obj.heat_history = struct('time', [], 'convection', [], 'radiation', [], 'reaction', [], 'reaction_surface', [], 'total', []);
            end
        end
        
        function results = solve(obj)
            % 求解整个过程
            
            % 清空历史记录
            obj.results.time = [];
            obj.results.particleStateHistory = [];
            obj.results.stageHistory = {};
            obj.results.flame_radius = [];      % 新增: 清空火焰半径历史
            obj.results.flame_temperature = []; % 新增: 清空火焰温度历史
            
            % 设置实时可视化参数
            obj.visualize_realtime = obj.params.visualize_realtime;
            
            current_time = 0;
            
            % 在初始化部分
            obj.params.current_stage = '';
            obj.params.consider_oxide_break = false;
            finish = false ;
            % 循环求解各阶段
            while current_time < obj.params.total_time
                % 1. 确定当前阶段
                current_stage_name = obj.model.determineStage();
                
                % 2. 根据阶段动态创建求解器并求解
                fprintf('当前时间: %.4fs, 进入 %s 阶段...\n', current_time, current_stage_name);
                
                t_stage = [];
                state_stage_history = ParticleState.empty;
                flame_info_history = []; % 新增: 为火焰信息历史记录提供默认空值

                switch current_stage_name
                    case 'heating_and_melting'
                        % 使用统一的、基于焓的求解器
                        solver = UnifiedHeatingStage(obj.params, obj.physicalModel);
                        target_temp = obj.params.materials.Mg.ignition_temp;
                        fprintf('--- 进入统一加热/熔化阶段 ---\n');
                        fprintf('  (目标: 点火温度 %.1f K)\n', target_temp);
                        
                        % 加热阶段才考虑氧化层破裂
                        obj.params.current_stage = 'heating_and_melting';
                        obj.params.consider_oxide_break = true;
                        
                        % 前置条件检查
                        if obj.model.particleState.T_p >= target_temp
                            fprintf('  (跳过: 当前温度 %.1f K >= 目标温度 %.1f K)\n', obj.model.particleState.T_p, target_temp);
                            t_stage = [current_time, current_time];
                            state_stage_history = [obj.model.particleState, obj.model.particleState];
                        else
                            fprintf('  加热: 当前温度 %.1f K ， 目标温度 %.1f K  继续加热\n', obj.model.particleState.T_p, target_temp);
                            [t_stage, state_stage_history] = solver.solve(obj.model.particleState, current_time, target_temp);
                        end

                    case 'vaporization'
                        % 1. 将颗粒温度锁定在沸点
                        T_boil = obj.params.materials.Mg.ignition_temp;
                        obj.model.particleState.T_p = T_boil;
                        fprintf('--- 进入气相燃烧阶段, 颗粒温度锁定在 %.1f K ---\n', T_boil);

                        % 气化阶段不考虑氧化层破裂
                        obj.params.current_stage = 'vaporization';
                        obj.params.consider_oxide_break = false;
                        
                        % 2. 初始化求解器
                        vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
                        
                        % 3. 创建VaporizationODE对象并传递参数
                        ode_solver = VaporizationODE(obj.params, obj.physicalModel);
                        
                        % 4. 调用VaporizationODE求解气相燃烧阶段
                        [t_stage, y_stage] = ode_solver.solve(obj.model.particleState, current_time);
                        
                        % 5. 将状态向量历史转换为状态对象历史, 并记录火焰信息
                        [state_stage_history, flame_info_history] = ode_solver.convert_state_vector_to_state_history(t_stage, y_stage, T_boil);
                        
                        % 6. 从VaporizationODE获取缓存以供结果处理
                        obj.vaporization_rate_info_cache = ode_solver.vaporization_rate_info_cache;

                end

                if isempty(t_stage)
                    fprintf('阶段 %s 未产生有效的时间步, 仿真可能已停滞。\n', current_stage_name);
                    break;
                end
                
                current_time = t_stage(end);
                fprintf('阶段 %s 完成, 结束时间: %.6f ms\n', current_stage_name, current_time*1000);
                switch current_stage_name
                    case 'vaporization'
                        finish =true;
                end
                
                % 3. 记录结果 (传入火焰信息)
                obj.recordResults(t_stage, state_stage_history, current_stage_name, flame_info_history);

                % 4. 关键: 用当前阶段的最终状态更新主模型
                obj.model.particleState = state_stage_history(end);
                
                % 5. 检查是否需要结束
                if current_time >= obj.params.total_time||finish
                    break;
                end
            end
            
            results = obj.prepareResults();
            Visualization = obj.params.visualization ;
            if Visualization
                visualize_combustion_process(results);
            end
            record = obj.params.recordResults;
            if record 
                write_results_for_origin(results);
            end

        end
        
        % 辅助函数: 三元操作符模拟
        function result = iif(condition, true_value, false_value)
            if condition
                result = true_value;
            else
                result = false_value;
            end
        end

        function recordResults(obj, t, state_history, stage_name, flame_info_history)
            % --- 开始修改: 增加对火焰信息和热量的记录 ---
            if nargin < 5
                flame_info_history = []; % 向后兼容
            end

            if isempty(t)
                return;
            end
            
            % 排除重复的第一个时间点
            if ~isempty(obj.results.time) && t(1) == obj.results.time(end)
                t = t(2:end);
                state_history = state_history(2:end);
                if ~isempty(flame_info_history)
                    flame_info_history = flame_info_history(2:end);
                end
            end

            if isempty(t)
                return;
            end
            
            obj.results.time = [obj.results.time; t(:)];
            obj.results.particleStateHistory = [obj.results.particleStateHistory, state_history];
            obj.results.stageHistory = [obj.results.stageHistory; repmat({stage_name}, length(t), 1)];
            
            % 记录火焰信息
            if ~isempty(flame_info_history)
                % V12: 确保flame_info_history是一个结构体数组
                if isstruct(flame_info_history)
                    obj.results.flame_radius = [obj.results.flame_radius; [flame_info_history.r_f]'];
                    obj.results.flame_temperature = [obj.results.flame_temperature; [flame_info_history.T_f]'];
                    
                    % 记录热量数据 (气相燃烧阶段)
                    if strcmp(stage_name, 'vaporization')
                        n_steps = length(flame_info_history);
                        heat_convection = zeros(n_steps, 1);
                        heat_radiation = zeros(n_steps, 1);
                        heat_reaction = zeros(n_steps, 1);
                        heat_reaction_surface = zeros(n_steps, 1);
                        heat_total = zeros(n_steps, 1);
                        
                        for i = 1:n_steps
                            % 直接使用flame_info_history中的热量数据
                            if isfield(flame_info_history(i), 'heat_convection')
                                heat_convection(i) = flame_info_history(i).heat_convection;
                                heat_radiation(i) = flame_info_history(i).heat_radiation;
                                heat_reaction(i) = flame_info_history(i).heat_reaction;
                                heat_reaction_surface(i) = flame_info_history(i).heat_reaction_surface;
                                heat_total(i) = flame_info_history(i).heat_total;
                            else
                                % 如果没有热量数据，计算基本热量
                                state = state_history(i);
                                T_p = state.T_p;
                                r_p = state.r_p;
                                
                                % 对流热量
                                q_conv = obj.params.h_conv * (obj.params.ambient_temperature - T_p);
                                heat_convection(i) = q_conv * (4 * pi * r_p^2);
                                
                                % 辐射热量
                                q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T_p^4);
                                heat_radiation(i) = q_rad * (4 * pi * r_p^2);
                                
                                % 无法计算反应热，设为0
                                heat_reaction(i) = 0;
                                heat_total(i) = heat_convection(i) + heat_radiation(i);
                            end
                        end
                        
                        % 添加到热量历史
                        obj.heat_history.time = [obj.heat_history.time; t(:)];
                        obj.heat_history.convection = [obj.heat_history.convection; heat_convection];
                        obj.heat_history.radiation = [obj.heat_history.radiation; heat_radiation];
                        obj.heat_history.reaction = [obj.heat_history.reaction; heat_reaction];
                        obj.heat_history.reaction_surface = [obj.heat_history.reaction_surface; heat_reaction_surface];
                        obj.heat_history.total = [obj.heat_history.total; heat_total];
                    end
                end
            else
                % 如果没有火焰信息，用NaN填充
                obj.results.flame_radius = [obj.results.flame_radius; nan(length(t), 1)];
                obj.results.flame_temperature = [obj.results.flame_temperature; nan(length(t), 1)];
            end
            
            % 记录加热阶段的热量数据
            if strcmp(stage_name, 'heating_and_melting')
                n_steps = length(t);
                heat_convection = zeros(n_steps, 1);
                heat_radiation = zeros(n_steps, 1);
                heat_reaction_surface = zeros(n_steps, 1);
                heat_reaction = zeros(n_steps, 1);
                heat_total = zeros(n_steps, 1);
                
                for i = 1:n_steps
                    % 计算对流换热
                    state = state_history(i);
                    T_p = state.T_p;
                    r_p = state.r_p;
                    
                    % 对流热量
                    h_conv = obj.params.k_gas ;
                    q_conv = h_conv / r_p * (obj.params.ambient_temperature - T_p);
                    heat_convection(i) = q_conv * (4 * pi * r_p^2);
                    
                    % 辐射热量
                    q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T_p^4);
                    heat_radiation(i) = q_rad * (4 * pi * r_p^2);
                    
                    % 反应热 (加热阶段通常没有反应热，除非温度足够高)
                    heat_reaction(i) = 0;
                    heat_reaction_surface(i) = 0 ;
                    if T_p >= obj.params.T_reaction_begin && i > 1
                        % 估算反应热（基于Mg消耗量）
                        if i > 1
                            dmg_dt = (state.m_mg - state_history(i-1).m_mg) / (t(i) - t(i-1));
                            if dmg_dt < 0 % 确保是消耗
                             %   n_mol = abs(dmg_dt) / obj.params.materials.Mg.molar_mass;
                                heat_reaction(i) = abs(dmg_dt * obj.params.reaction_heat_face);
                                heat_reaction_surface(i) = abs(dmg_dt * obj.params.reaction_heat_face);
                            end
                        end
                    end
                    
                    heat_total(i) = heat_convection(i) + heat_radiation(i) + heat_reaction(i);
                end
                
                % 添加到热量历史
                obj.heat_history.time = [obj.heat_history.time; t(:)];
                obj.heat_history.convection = [obj.heat_history.convection; heat_convection];
                obj.heat_history.radiation = [obj.heat_history.radiation; heat_radiation];
                obj.heat_history.reaction = [obj.heat_history.reaction; heat_reaction];
                obj.heat_history.reaction_surface = [obj.heat_history.reaction_surface; heat_reaction_surface];
                obj.heat_history.total = [obj.heat_history.total; heat_total];
            end
            % --- 结束修改 ---
            
            % 调试输出
            record = false ;
            if record
                if strcmp(stage_name, 'heating_and_melting')
                    fprintf('加热阶段记录结果：\n');
                    for i = 1:min(3, length(state_history)) % 只打印前3个点
                        fprintf('  时间=%.4fs, T=%.1fK, oxide=%.2e m\n', ...
                            t(i), state_history(i).T_p, state_history(i).oxide_thickness);
                    end
                end
            end
        end

        % 新增函数: 初始化实时可视化
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
        
        function results = prepareResults(obj)
            % 辅助函数: 将历史记录转换为易于绘图和分析的格式
            num_entries = length(obj.results.time);
            if num_entries == 0
                results = struct();
                return;
            end

            % 初始化数组
            results.time = obj.results.time;
            results.stage = obj.results.stageHistory; % 保留原始的计算阶段
            results.physical_stage = cell(num_entries, 1); % 新增物理阶段描述
            results.temperature = zeros(num_entries, 1);
            results.melted_fraction = zeros(num_entries, 1);
            results.mass_mg = zeros(num_entries, 1);
            results.mass_mgo = zeros(num_entries, 1);
            results.mass_c = zeros(num_entries, 1); % 新增: 碳质量
            results.radius = zeros(num_entries, 1);
            results.oxide_thickness = zeros(num_entries, 1);
            results.flame_radius = []; % 新增: 火焰半径
            results.flame_temperature = []; % 新增: 火焰温度

            % 热量数据数组
            results.heat_convection = zeros(num_entries, 1);
            results.heat_radiation = zeros(num_entries, 1);
            results.heat_reaction = zeros(num_entries, 1);
            results.heat_reaction_surface = zeros(num_entries, 1);
            results.heat_total = zeros(num_entries, 1);
            
            % 变化率数据数组
            results.dmdt_mg = zeros(num_entries, 1);
            results.dmdt_mgo = zeros(num_entries, 1);
            results.dmdt_c = zeros(num_entries, 1);
            results.drdt_c = zeros(num_entries, 1);
            results.drdt_p = zeros(num_entries, 1);

            % 填充数据
            for i = 1:num_entries
                state = obj.results.particleStateHistory(i);
                results.temperature(i) = state.T_p;
                results.melted_fraction(i) = state.melted_fraction;
                results.mass_mg(i) = state.m_mg;
                results.mass_mgo(i) = state.m_mgo;
                results.mass_c(i) = state.m_c; % 新增
                results.radius(i) = state.r_p;
                results.oxide_thickness(i) = state.oxide_thickness;
                
                % 根据温度和熔化分数追认物理阶段
                if strcmp(results.stage{i}, 'heating_and_melting')
                    T = state.T_p;
                    X = state.melted_fraction;
                    T_melt = obj.params.materials.Mg.melting_point;

                    if T < T_melt - 1e-6
                        results.physical_stage{i} = 'Preheating';
                    elseif abs(T - T_melt) < 1e-6 && X < 1.0
                        results.physical_stage{i} = 'Melting';
                    else
                        results.physical_stage{i} = 'Liquid Heating';
                    end
                else
                    % 对于其他阶段(如vaporization), 直接使用计算阶段名
                    results.physical_stage{i} = results.stage{i};
                end
            end
            
            % 直接从记录中附加火焰数据
            if isfield(obj.results, 'flame_radius')
                results.flame_radius = obj.results.flame_radius;
                results.flame_temperature = obj.results.flame_temperature;
            end
            
            % 使用记录的热量历史数据
            if ~isempty(obj.heat_history.time)
                % 将热量历史数据插值到结果时间点
                results.heat_convection = interp1(obj.heat_history.time, obj.heat_history.convection, results.time, 'linear', 'extrap');
                results.heat_radiation = interp1(obj.heat_history.time, obj.heat_history.radiation, results.time, 'linear', 'extrap');
                results.heat_reaction = interp1(obj.heat_history.time, obj.heat_history.reaction, results.time, 'linear', 'extrap');
                results.heat_reaction_surface = interp1(obj.heat_history.time, obj.heat_history.reaction_surface, results.time, 'linear', 'extrap');
                results.heat_total = interp1(obj.heat_history.time, obj.heat_history.total, results.time, 'linear', 'extrap');
            end
            
            % 计算变化率数据
            % 对于气相燃烧阶段，从vaporization_rate_info_cache中获取变化率
            if ~isempty(obj.vaporization_rate_info_cache) && strcmp(results.stage{end}, 'vaporization')
                % 仅处理气相燃烧阶段的最后一部分数据
                vap_indices = find(strcmp(results.stage, 'vaporization'));
                if ~isempty(vap_indices)
                    % 确保缓存长度足够
                    cache_len = length(obj.vaporization_rate_info_cache);
                    for i = 1:min(length(vap_indices), cache_len)
                        idx = vap_indices(i);
                        if idx <= num_entries
                            cache_idx = min(i, cache_len);
                            rate_info = obj.vaporization_rate_info_cache{cache_idx};
                            
                            % 提取变化率数据（如果存在）
                            if isfield(rate_info, 'dmdt_mg')
                                results.dmdt_mg(idx) = rate_info.dmdt_mg * 1e6;  % 转为μg/s
                                results.dmdt_mgo(idx) = rate_info.dmdt_mgo * 1e6;
                                results.dmdt_c(idx) = rate_info.dmdt_c * 1e6;
                            end
                        end
                    end
                end
            end
            
            % 计算半径变化率
            % 使用中心差分计算变化率
            for i = 2:num_entries-1
                dt = results.time(i+1) - results.time(i-1);
                if dt > 0
                    % 计算Mg核半径变化率
                    dr_c = obj.results.particleStateHistory(i+1).r_c - obj.results.particleStateHistory(i-1).r_c;
                    results.drdt_c(i) = (dr_c / dt) * 1e6;  % 转为μm/s
                    
                    % 计算颗粒半径变化率
                    dr_p = obj.results.particleStateHistory(i+1).r_p - obj.results.particleStateHistory(i-1).r_p;
                    results.drdt_p(i) = (dr_p / dt) * 1e6;  % 转为μm/s
                end
            end
            
            % 处理首尾两点（前向差分和后向差分）
            if num_entries > 1
                % 首点
                dt = results.time(2) - results.time(1);
                if dt > 0
                    dr_c = obj.results.particleStateHistory(2).r_c - obj.results.particleStateHistory(1).r_c;
                    results.drdt_c(1) = (dr_c / dt) * 1e6;
                    
                    dr_p = obj.results.particleStateHistory(2).r_p - obj.results.particleStateHistory(1).r_p;
                    results.drdt_p(1) = (dr_p / dt) * 1e6;
                end
                
                % 尾点
                dt = results.time(end) - results.time(end-1);
                if dt > 0
                    dr_c = obj.results.particleStateHistory(end).r_c - obj.results.particleStateHistory(end-1).r_c;
                    results.drdt_c(end) = (dr_c / dt) * 1e6;
                    
                    dr_p = obj.results.particleStateHistory(end).r_p - obj.results.particleStateHistory(end-1).r_p;
                    results.drdt_p(end) = (dr_p / dt) * 1e6;
                end
            end

            % 新增: 计算各阶段的时间统计
            stage_times = calculate_stage_durations(results);
            results.stage_times = stage_times;
        end
    end
end 

function visualize_combustion_process(results)
    % 创建四子图布局
    figure('Position', [100, 100, 1200, 900], 'Name', '颗粒燃烧过程综合分析');
    
    % 1. 温度变化
    subplot(2, 2, 1);
    plot(results.time*1000, results.temperature, 'r-', 'LineWidth', 2);
    hold on;
    if isfield(results, 'flame_temperature') && ~isempty(results.flame_temperature)
        valid_idx = ~isnan(results.flame_temperature);
        if any(valid_idx)
            plot(results.time(valid_idx)*1000, results.flame_temperature(valid_idx), 'r--', 'LineWidth', 1.5);
        end
    end
    xlabel('时间 (s)');
    ylabel('温度 (K)');
    title('温度变化');
    grid on;
    legend('颗粒温度', '火焰温度');
    
    % 2. 尺寸变化
    subplot(2, 2, 2);
    yyaxis left
    plot(results.time*1000, results.radius*1e6, 'b-', 'LineWidth', 2);
    ylabel('颗粒半径 (μm)');
    
    yyaxis right
    if isfield(results, 'oxide_thickness')
        plot(results.time*1000, results.oxide_thickness*1e6, 'g-', 'LineWidth', 1.5);
    end
    ylabel('氧化层厚度 (μm)');
    
    if isfield(results, 'flame_radius') && ~isempty(results.flame_radius)
        yyaxis left
        hold on;
        valid_idx = ~isnan(results.flame_radius);
        if any(valid_idx)
            plot(results.time(valid_idx)*1000, results.flame_radius(valid_idx)*1e6, 'r--', 'LineWidth', 1.5);
        end
    end
    
    xlabel('时间 (ms)');
    title('尺寸变化');
    grid on;
    legend('颗粒半径', '氧化层', '火焰半径');
    
    % 3. 质量变化
    subplot(2, 2, 3);
    plot(results.time*1000, results.mass_mg*1e6, 'r-', 'LineWidth', 2);
    hold on;
    plot(results.time*1000, results.mass_mgo*1e6, 'b-', 'LineWidth', 2);
    if isfield(results, 'mass_c')
        plot(results.time*1000, results.mass_c*1e6, 'k-', 'LineWidth', 2);
        total_mass = results.mass_mg + results.mass_mgo + results.mass_c;
    else
        total_mass = results.mass_mg + results.mass_mgo;
    end
    plot(results.time*1000, total_mass*1e6, 'g--', 'LineWidth', 1.5);
    xlabel('时间 (ms)');
    ylabel('质量 (μg)');
    title('组分质量变化');
    grid on;
    if isfield(results, 'mass_c')
        legend('Mg金属', 'MgO', '碳', '总质量');
    else
        legend('Mg金属', 'MgO', '总质量');
    end
    
    % 4. 热量来源分析 - 使用双Y轴
    subplot(2, 2, 4);

    % 左Y轴：对流和辐射热
    yyaxis left
    plot(results.time*1000, results.heat_convection, 'b-', 'LineWidth', 2, 'DisplayName', '对流换热');
    hold on;
    plot(results.time*1000, results.heat_radiation, 'r-', 'LineWidth', 2, 'DisplayName', '辐射换热');
    ylabel('对流和辐射热量 (W)');
    set(gca, 'YColor', [0 0 0.7]); % 深蓝色

    % 右Y轴：反应热和总热量
    yyaxis right
    plot(results.time*1000, results.heat_reaction, 'g-', 'LineWidth', 2, 'DisplayName', '反应热');
    plot(results.time*1000, results.heat_total, 'k--', 'LineWidth', 2, 'DisplayName', '总热量');
    ylabel('反应热和总热量 (W)');
    set(gca, 'YColor', [0 0.7 0]); % 深绿色

    % 共享设置
    xlabel('时间 (ms)');
    title('热量来源分析');
    grid on;
    legend('Location', 'best');

    % 添加物理阶段标记
    mark_physical_stages(results);
    
    % 整体标题
    sgtitle('颗粒燃烧过程综合分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 新增: 创建阶段时间统计图
    if isfield(results, 'stage_times')
        figure('Name', '各物理阶段时间统计', 'Position', [100, 500, 800, 400]);
        
        % 提取各阶段的时间
        stage_names = {'Preheating', 'Melting', 'Liquid_Heating', 'vaporization', 'Unclassified'};
        stage_labels = {'预热阶段', '熔融阶段', '液相加热', '气相燃烧', '未分类'};
        stage_durations = zeros(1, length(stage_names));
        for i = 1:length(stage_names)
            if isfield(results.stage_times, stage_names{i})
                stage_durations(i) = results.stage_times.(stage_names{i}) * 1000; % 转换为ms
            end
        end
        
        % 绘制柱状图
        bar(stage_durations, 0.6);
        set(gca, 'XTickLabel', stage_labels);
        ylabel('持续时间 (ms)');
        title('颗粒燃烧各物理阶段时间统计');
        grid on;
        
        % 显示具体数值
        for i = 1:length(stage_durations)
            if stage_durations(i) > 0
                text(i, stage_durations(i) + max(stage_durations)*0.03, ...
                     sprintf('%.2f ms', stage_durations(i)), ...
                     'HorizontalAlignment', 'center');
            end
        end
        
        % 添加文本框显示总时间
        total_time_ms = results.stage_times.Total * 1000;
        classified_time_ms = (results.stage_times.Total - results.stage_times.Unclassified) * 1000;
        annotation('textbox', [0.15, 0.8, 0.3, 0.15], ...
                  'String', {sprintf('总时间: %.2f ms', total_time_ms), ...
                             sprintf('已分类时间: %.2f ms (%.1f%%)', classified_time_ms, ...
                                    classified_time_ms/total_time_ms*100)}, ...
                  'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
                  'EdgeColor', 'black', 'FontWeight', 'bold');
    end
    
    % 在visualize_combustion_process中添加氧化层破裂可视化
    if isfield(results, 'oxide_break_history') && ~isempty(results.oxide_break_history.time)
        % 确保只展示加热阶段的破裂数据
        valid_indices = [];
        for i = 1:length(results.oxide_break_history.time)
            % 找出每个破裂时间点对应的阶段
            break_time = results.oxide_break_history.time(i);
            % 找到时间最接近的索引
            [~, time_idx] = min(abs(results.time - break_time));
            if time_idx <= length(results.stage) && strcmp(results.stage{time_idx}, 'heating_and_melting')
                valid_indices = [valid_indices, i];
            end
        end
        
        if ~isempty(valid_indices)
            figure('Name', '氧化层破裂历史', 'Position', [100, 100, 800, 600]);
            plot(results.oxide_break_history.time(valid_indices), results.oxide_break_history.temperature(valid_indices), 'r-', 'LineWidth', 2);
            hold on;
            plot(results.oxide_break_history.time(valid_indices), results.oxide_break_history.stress(valid_indices), 'b-', 'LineWidth', 2);
            xlabel('时间 (s)');
            ylabel('温度 (K) 和 应力 (Pa)');
            title('氧化层破裂历史 (仅加热阶段)');
            grid on;
            legend('温度', '应力');
        end
    end
end

function mark_physical_stages(results)
    % 在每个子图上标记物理阶段边界
    if ~isfield(results, 'physical_stage') || isempty(results.physical_stage)
        return;
    end
    
    % 获取阶段转换点
    stages = results.physical_stage;
    transitions = [];
    transition_names = {};
    
    for i = 1:length(stages)-1
        if ~strcmp(stages{i}, stages{i+1})
            transitions = [transitions; results.time(i)*1000];
            transition_names{end+1} = sprintf('%s→%s', stages{i}, stages{i+1});
        end
    end
    
    % 在每个子图上添加阶段分隔线
    subplots = get(gcf, 'Children');
    for i = 1:length(subplots)
        if strcmp(get(subplots(i), 'Type'), 'axes')
            axes(subplots(i));
            ylim_current = get(gca, 'YLim');
            
            for j = 1:length(transitions)
                line([transitions(j) transitions(j)], ylim_current, ...
                     'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                
                % 添加标签，位置在y轴中部
                text(transitions(j), mean(ylim_current), transition_names{j}, ...
                     'Rotation', 90, 'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'middle', 'FontSize', 8);
            end
        end
    end
end 

function stage_times = calculate_stage_durations(results)
    % 计算各物理阶段的持续时间
    % 初始化阶段时间结构体
    stage_times = struct(...
        'Preheating', 0, ...
        'Melting', 0, ...
        'Liquid_Heating', 0, ...
        'vaporization', 0, ...
        'Total', 0, ...
        'Unclassified', 0);  % 添加未分类时间字段
    
    if isempty(results.time) || length(results.time) < 2
        return;
    end
    
    % 获取物理阶段和时间点
    physical_stages = results.physical_stage;
    times = results.time;
    
    % 当前阶段和起始时间
    current_stage = physical_stages{1};
    stage_start = times(1);
    total_classified_time = 0;
    
    % 遍历所有时间点
    for i = 2:length(times)
        if ~strcmp(physical_stages{i}, current_stage)
            % 阶段变化，计算持续时间
            % 为避免间隙，认为阶段边界在当前时间点和前一点之间的中点
            mid_time = (times(i) + times(i-1)) / 2;
            duration = mid_time - stage_start;
            
            % 将阶段名中的空格替换为下划线，以便作为结构体字段名
            field_name = strrep(current_stage, ' ', '_');
            
            % 累加当前阶段的持续时间
            if isfield(stage_times, field_name)
                stage_times.(field_name) = stage_times.(field_name) + duration;
            else
                stage_times.(field_name) = duration;
            end
            total_classified_time = total_classified_time + duration;
            
            % 更新当前阶段和起始时间
            current_stage = physical_stages{i};
            stage_start = mid_time;
        end
        
        % 处理最后一个时间点
        if i == length(times)
            final_duration = times(end) - stage_start;
            field_name = strrep(current_stage, ' ', '_');
            if isfield(stage_times, field_name)
                stage_times.(field_name) = stage_times.(field_name) + final_duration;
            else
                stage_times.(field_name) = final_duration;
            end
            total_classified_time = total_classified_time + final_duration;
        end
    end
    
    % 计算总时间和未分类时间
    stage_times.Total = times(end) - times(1);
    stage_times.Unclassified = stage_times.Total - total_classified_time;
    
    % 确保未分类时间不为负（可能由于舍入误差）
    if stage_times.Unclassified < 0 && stage_times.Unclassified > -1e-10
        stage_times.Unclassified = 0;
    end
end 

function write_results_for_origin(results)
    % write_results_for_origin  将燃烧模拟结果写成 Origin 可直接读取的 CSV
    %   results  : 由主程序生成的结构体
    %   输出文件名 = result/初始颗粒直径(μm).csv （自动覆盖同名文件）
    
        %% 1. 自动生成文件名
        d0_um   = results.radius(1) * 2 * 1e6;      % 初始直径，μm
        fileName = sprintf('%.1fμm.csv', d0_um);
        
        % 确保results文件夹存在
        resultDir = 'result';
        if ~exist(resultDir, 'dir')
            mkdir(resultDir);
        end
        
        % 添加文件夹路径
        fullPath = fullfile(resultDir, fileName);
    
        %% 2. 准备标题（cell）和数据（double）
        headers = {'时间(ms)', '颗粒温度(K)'};
        data    = [results.time(:)*1e3, results.temperature(:)];  % 单位转换
    
        % 可选字段：存在就写，不存在填 NaN
        optFields = {
            'flame_temperature',        '火焰温度(K)';
            'radius',                   '颗粒半径(um)';
            'oxide_thickness',          '氧化层厚度(um)';
            'flame_radius',             '火焰半径(um)';
            'mass_mg',                  'Mg质量(ug)';
            'mass_mgo',                 'MgO质量(ug)';
            'mass_c',                   '碳质量(ug)';
            'heat_convection',          '对流换热(W)';
            'heat_radiation',           '辐射换热(W)';
            'heat_reaction',            '反应热(W)';
            'heat_reaction_surface',    '表面反应热(W)';
            'heat_total',               '总热量(W)';
            % 新增变化率数据字段
            'dmdt_mg',                  'Mg质量变化率(ug/s)';
            'dmdt_mgo',                 'MgO质量变化率(ug/s)';
            'dmdt_c',                   'C质量变化率(ug/s)';
            'drdt_c',                   'Mg核半径变化率(um/s)';
            'drdt_p',                   '颗粒半径变化率(um/s)';
        };
    
        nRow = numel(results.time);
        for k = 1:size(optFields,1)
            key   = optFields{k,1};
            label = optFields{k,2};
    
            if isfield(results, key) && ~isempty(results.(key))
                v = results.(key)(:);
                % 单位转换
                if contains(label, 'um'); v = v * 1e6;  end
                if contains(label, 'ug'); v = v * 1e6;  end
            else
                v = NaN(nRow,1);
            end
    
            headers{end+1} = label;          % 标题 cell
            data = [data, v];                % 数值矩阵
        end
    
        %% 3. 写 CSV：首行标题（cell），其余行数据（double）
        T = [headers; num2cell(data)];      % 组合成 (1+nRow) × nCol cell
        writecell(T, fullPath, 'Delimiter', 'comma', 'FileType', 'text', 'WriteMode', 'overwrite');
        fprintf('结果已写入 %s（共 %d 列），可直接用 Origin 打开。\n', fullPath, size(data,2));
    end