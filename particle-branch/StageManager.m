classdef StageManager < handle
    % StageManager 阶段管理器
    % 管理不同阶段的切换和求解 (动态调度器)
    
    properties
        model         % 颗粒模型 (ParticleModel)
        params        % 参数对象 (Parameters)
        physicalModel % 物理计算模型 (PhysicalModel)
        thermo_reader % 热力学数据读取器 (ThermoReader)
        
        % 结果存储
        results       % 包含所有历史记录的结构体
        vaporization_rate_info_cache % V12: 用于在气化阶段ODE求解中缓存BVP结果
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
                        
                        % 2. 初始化求解器和结果缓存
                        vaporization_solver = VaporizationStage(obj.params, obj.physicalModel);
                        obj.vaporization_rate_info_cache = {}; % V12: 初始化缓存
                        
                        % 3. 设置ODE的初始状态向量 [m_mg, m_mgo, m_c, r_c, r_p, T_p]
                        pState = obj.model.particleState;
                        y0 = [pState.m_mg, pState.m_mgo, pState.m_c, pState.r_c, pState.r_p, pState.T_p];

                        % 4. 根据参数选择使用固定步长或ODE45求解
                        if obj.params.use_fixed_timestep
                            % 使用固定步长求解
                            fprintf('使用固定步长求解器, 步长: %.2e s\n', obj.params.vaporization_fixed_timestep);
                            t_span = [current_time, obj.params.total_time];
                            [t_stage, y_stage] = obj.solve_with_fixed_timestep(...
                                @(t,y) obj.vaporization_ode(t, y, vaporization_solver, T_boil), ...
                                t_span, y0, obj.params.vaporization_fixed_timestep);
                        else
                            % 使用原有的ODE45求解
                            fprintf('使用ODE45自适应步长求解器\n');
                            t_span = [current_time, obj.params.total_time];
                            output_fcn_handle = @(t,y,flag) obj.vaporization_output_fcn(t, y, flag, vaporization_solver, T_boil);
                            ode_options = odeset('RelTol', 1e-6, 'Events', @(t,y) obj.vaporization_events(t,y), 'OutputFcn', output_fcn_handle,'MaxStep', 0.001, 'InitialStep', 0.0005);
                            
                            % 使用ODE45求解
                            [t_stage, y_stage] = ode45(@(t,y) obj.vaporization_ode(t, y, vaporization_solver, T_boil), t_span, y0, ode_options);
                        end
                        
                        % 6. 将状态向量历史转换为状态对象历史, 并记录火焰信息
                        [state_stage_history, flame_info_history] = obj.convert_state_vector_to_state_history(t_stage, y_stage, T_boil, obj.vaporization_rate_info_cache);
                        
                        % 7. 清理缓存
                        obj.vaporization_rate_info_cache = {};

                        % --- 结束修改 ---
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
        
        % 新增: 使用固定步长的显式欧拉法求解ODE
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
            obj.vaporization_output_fcn(t(1), y(1,:)', 'init', [], T_boil);
            
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
                obj.vaporization_output_fcn(t(i), y(i,:)', '', [], T_boil);
                
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
                dydt = compute_dydt_from_rate_info(obj, t(i), y(i,:)', rate_info, T_boil);
                
                % 更新下一个状态
                y(i+1,:) = y(i,:) + dt * dydt';
                
                % 保存到全局缓存，供可视化使用
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
            obj.vaporization_output_fcn(t(end), y(end,:)', 'done', [], T_boil);
        end

        function dydt = compute_dydt_from_rate_info(obj, t, y, rate_info, T_p_const)
            % 基于已有rate_info计算ODE右侧函数值
            % 这样就不需要重新调用求解器了
            
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
                C_oxide = calculate_oxide_heat_capacity(obj, y);
                if C_oxide > 1e-10
                    dydt(6) = rate_info.heat_ox / C_oxide; % 温度变化率
                    DT = obj.params.vaporization_fixed_timestep ;
                    DT_reaction = rate_info.heat_reaction_total / C_oxide;
                    DT_reaction_face= rate_info.heat_reaction_surface / C_oxide;
                    DT_reaction_flame= rate_info.heat_reaction_flame / C_oxide;
                    fprintf('  mg质量变化所需热量 : %.3e \n', rate_info.dmdt_mg * obj.params.materials.Mg.L_evap_Mg /...
                        obj.params.materials.Mg.molar_mass * DT);
                    fprintf('  heat_ox: %.3e \n', rate_info.heat_ox);
                    fprintf('  C_oxide: %.3e \n', C_oxide);
                    fprintf('  dydt(6): %.3e \n', dydt(6));
                    fprintf('  计算采用: %.3e \n', dydt(6)*DT);
                    fprintf('  总反应反应部分带来温升: %.3e \n', DT_reaction*DT);
                    fprintf('  表面反应部分带来温升: %.3e \n', DT_reaction_face*DT);
                    fprintf('  火焰反应部分带来温升: %.3e \n', DT_reaction_flame*DT);
                    fprintf('  heat_reaction_flame: %.3e \n', rate_info.heat_reaction_flame*DT);

                else
                    dydt(6) = 0;
                end
            else
                dydt(6) = 0;
            end
        end
        
        % 辅助函数: 计算氧化层热容量
        function C_oxide = calculate_oxide_heat_capacity(obj, y)
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
        
        % 辅助函数: 三元操作符模拟
        function result = iif(condition, true_value, false_value)
            if condition
                result = true_value;
            else
                result = false_value;
            end
        end

        function dydt = vaporization_ode(obj, t, y, solver, T_p_const)
            % --- V11 重构: 求解包含质量和半径的5变量状态向量 ---
            % y(1)=m_mg, y(2)=m_mgo, y(3)=m_c, y(4)=r_c, y(5)=r_p, y(6)=T_p
            
            % 检查是否已经有缓存的计算结果
            persistent last_t last_rate_info last_algebraic_solution
            fprintf('ODE 求解器在 t=%.4f ms 调用了 vaporization_ode\n', t*1000);
    ...
           % % 如果当前时间点已计算过且差异很小，直接使用缓存结果
           % if ~isempty(last_t) && abs(t - last_t) < 1e-10 && ~isempty(last_rate_info)
           %     % 使用缓存的结果
           %     fprintf('使用缓存结果\n');
           %     dmdt_mg = last_rate_info.dmdt_mg;
            %    dmdt_mgo = last_rate_info.dmdt_mgo;
            %    dmdt_c = last_rate_info.dmdt_c;
            %else
                % 为防止数值问题，增加保护性判断
                %if y(1) < 1e-15 || y(4) < 1e-9 || y(5) < 1e-9
                 %   dydt = zeros(5,1);
                %    return;
                %end

                % 1. 从当前状态向量y构建一个临时的ParticleState对象
                % 注意：这里不再调用update_geometry_from_mass
                tempState = obj.model.particleState.copy();
                tempState.m_mg = y(1);
                tempState.m_mgo = y(2);
                tempState.m_c = y(3);
                tempState.r_c = y(4);
                tempState.r_p = y(5);
                tempState.T_p = y(6);   % 使用当前温度而不是固定温度
                tempState.oxide_thickness = tempState.r_p - tempState.r_c;
                tempState.T_c = T_p_const;  % 金属核心温度锁定在沸点
                
                % 2. 获取质量变化率
                try
                    % 调用求解器并传入上一步的代数解作为缓存
                    rate_info = solver.solve_reaction_rates(tempState, last_algebraic_solution);
                    
                    % 改进: 检查求解是否成功
                    if ~rate_info.success
                        error('扩散燃烧部分求解不收敛：在t=%.4f s处无法获得有效解', t);
                    end
                    
                    % 保存代数方程解用于下一时间步
                    if isfield(rate_info, 'algebraic_solution')
                        last_algebraic_solution = rate_info.algebraic_solution;
                    end
                    
                    dmdt_mg  = rate_info.dmdt_mg;
                    dmdt_mgo = rate_info.dmdt_mgo;
                    dmdt_c   = rate_info.dmdt_c;
                    
                    % 3. 计算表面温度变化率
                    % 获取热量信息
                    Q_total = rate_info.heat_ox;  % 总热量流入（W）
                    
                    % 计算颗粒热容 (只考虑氧化层)
                    % 核心部分温度固定，不需考虑热容
                    porosity = obj.params.material_properties.oxide_porosity;
                    
                    % 获取各组分密度
                    rho_mgo = obj.params.materials.MgO.density;
                    rho_c = obj.params.materials.C.density;
                    
                    % 计算氧化层体积与碳沉积体积
                    V_total = (4/3)*pi*(y(5)^3 - y(4)^3);  % 氧化层总体积
                    
                    % 计算碳质量分数 (基于总氧化层质量)
                    m_mgo = y(2);  % 氧化镁总质量
                    m_c = y(3);    % 碳总质量
                    m_total = m_mgo + m_c;  % 氧化层总质量
                    
                    % 计算各组分体积占比
                    if m_total > 0
                        mass_frac_c = m_c / m_total;
                        mass_frac_mgo = m_mgo / m_total;
                    else
                        mass_frac_c = 0;
                        mass_frac_mgo = 1;
                    end
                    
                    % 使用PhysicalModel中的方法计算摩尔比热，转换为质量比热
                    cp_mgo_molar = obj.physicalModel.calc_cp_solid_mgo(y(6));  % J/(mol·K)
                    cp_c_molar = obj.physicalModel.calc_cp_solid_c(y(6));      % J/(mol·K)
                    
                    % 转换为质量比热
                    cp_mgo_mass = cp_mgo_molar / obj.params.materials.MgO.molar_mass;  % J/(kg·K)
                    cp_c_mass = cp_c_molar / obj.params.materials.C.molar_mass;        % J/(kg·K)
                    
                    % 计算氧化层有效质量 (考虑孔隙度)
                    m_oxide_eff = (1-porosity) * ((m_mgo/rho_mgo) + (m_c/rho_c)) * (rho_mgo*mass_frac_mgo + rho_c*mass_frac_c);
                    m_oxide_eff = m_total ;
                    %m_oxide_eff = 1 ;
                    % 计算混合比热
                    cp_mix = mass_frac_mgo * cp_mgo_mass + mass_frac_c * cp_c_mass;
                    
                    % 计算氧化层有效热容Q_total  j /s
                    C_oxide = m_oxide_eff * cp_mix;  % [J/K]
                    %fprintf(' Q_total: %.3e \n', Q_total);
                    %fprintf(' C_oxide: %.3e \n', C_oxide);
                    %fprintf(' cp_mix: %.3e \n', cp_mix);
                    %fprintf(' m_oxide_eff: %.3e \n', m_oxide_eff);
                    % 表面温度变化率 (dT/dt = Q/C)
                    if C_oxide > 1e-10
                        dTp_dt = Q_total / C_oxide;
                    else
                        dTp_dt = 0;  % 避免除零
                    end
                    fprintf('t=%.3f ms: cp_mgo=%.2f J/(kg·K), cp_c=%.2f J/(kg·K), cp_mix=%.2f J/(kg·K)\n', t*1000, cp_mgo_mass, cp_c_mass, cp_mix);
                    fprintf('  m_mgo=%.2e kg, m_c=%.2e kg, C_oxide=%.2e J/K, dT/dt=%.5e K/s\n', m_mgo, m_c, C_oxide, dTp_dt);
                    % 输出调试信息
                    if isfield(obj.params, 'debug') && obj.params.debug && mod(round(t*1000), 100) == 0
                        fprintf('t=%.3f ms: cp_mgo=%.2f J/(kg·K), cp_c=%.2f J/(kg·K), cp_mix=%.2f J/(kg·K)\n', t*1000, cp_mgo_mass, cp_c_mass, cp_mix);
                        fprintf('  m_mgo=%.2e kg, m_c=%.2e kg, C_oxide=%.2e J/K, dT/dt=%.2f K/s\n', m_mgo, m_c, C_oxide, dTp_dt);
                    end
                    
                    % 更新缓存
                    last_t = t;
                    last_rate_info = rate_info;
                    
                    % 同时更新全局缓存，供outputfcn使用
                    if isempty(obj.vaporization_rate_info_cache)
                        obj.vaporization_rate_info_cache = {};
                    end
                    obj.vaporization_rate_info_cache{end+1} = rate_info;
                catch ME
                    fprintf('t=%.4f s 时反应速率求解失败: %s\n', t, ME.message);
                    % 重新抛出错误，这将中止ODE求解器
                    rethrow(ME);
                end
            %end

            % 3. 新增核心逻辑: 计算半径变化率
           % rho_mg = obj.params.materials.Mg.density;

            %if tempState.T_p < 923
            %    rho_mg = obj.params.materials.Mg.density_low;
            %else
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
            %fprintf('  t=%.4f ms, dmdt_mg=%.2e, dmdt_mgo=%.2e, drp/dt=%.2e\n', t*1000, dmdt_mg, dmdt_mgo, dmdt_c);
        end
        
        function [value, isterminal, direction] = vaporization_events(obj, t, y)
            % ODE事件函数: 当Mg质量降至初始值的1%以下时停止
            
            % 获取当前Mg质量 (y(1)是状态向量中的Mg质量)
            current_mg_mass = y(1);
            fprintf('  t=%.4f ms, current_mg_mass=%.2e\n', t*1000, current_mg_mass);
            % 获取初始Mg质量 (从结果的第一个时间点获取)
            %if isempty(obj.results.time)
            %    % 如果还没有结果，使用当前模型的初始值
            %    initial_mg_mass = obj.model.particleState.m_mg;
            %else
            %    % 使用结果中记录的初始值
             %   initial_mg_mass = obj.results.mass_mg(1);
            %end
            initial_mg_mass = obj.model.particleState.m_mg;
            fprintf('  t=%.4f ms, initial_mg_mass=%.2e\n', t*1000, initial_mg_mass);
            % 计算质量比例
            mass_ratio = current_mg_mass / initial_mg_mass;
            
            % 当质量比例低于1%时触发事件
            threshold = obj.params.combustionfinish_ratio ;
            % threshold = 0.01; % 1%阈值
            value = mass_ratio - threshold;
            
            isterminal = 1;    % 触发后停止积分
            direction = -1;    % 质量比例降低穿过阈值时触发
            %fprintf('t=%.6f ms\n, Mg比例=%.6f, threshold=%.3f%%\n', t *1000, mass_ratio*100, threshold*100);
            % 输出调试信息
            %if mass_ratio < 0.05  % 当接近阈值时输出信息
            %    fprintf('  t=%.6f ms: Mg质量比例 = %.4f%% (%.2e/%.2e kg)\n', ...
            %        t*1000, mass_ratio*100, current_mg_mass, initial_mg_mass);
            %
                
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

        function [state_history, flame_info_history] = convert_state_vector_to_state_history(obj, t_stage, y_stage, T_p_const, rate_info_cache)
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
            if isempty(rate_info_cache) || length(rate_info_cache) == 0
                fprintf('警告: 缓存为空，创建默认缓存项\n');
                % 创建一个默认的rate_info结构
                default_rate_info = obj.create_default_rate_info(obj.build_temp_state_from_vector(y_stage(1,:), T_p_const));
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
                    % 如果缓存为空但长度不为0（可能是{[]}这种情况）
                    last_valid = obj.create_default_rate_info(obj.build_temp_state_from_vector(y_stage(1,:), T_p_const));
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
                % 1. 从状态向量填充ParticleState对象，正确处理温度
                currentState = y_stage(i, :);
                
                % 使用新的build_temp_state_from_vector函数，将表面温度与核心温度分开处理
                % 状态向量中的第6个元素是表面温度T_p，而核心温度T_c保持为沸点温度
                tempState = obj.build_temp_state_from_vector(currentState, T_p_const);
                
                % 确保温度值在合理范围内
                if tempState.T_p < 0 || isnan(tempState.T_p)
                    tempState.T_p = T_p_const; % 如果温度无效，使用沸点温度作为后备值
                end
                
                state_history(i) = tempState;
                
                % 2. 从缓存中读取火焰信息
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
                
                % 3. 复制热量数据
                if isfield(rate_info, 'heat_convection')
                    flame_info_history(i).heat_convection = rate_info.heat_convection;
                    flame_info_history(i).heat_radiation = rate_info.heat_radiation;
                    flame_info_history(i).heat_reaction = rate_info.heat_reaction;
                    flame_info_history(i).heat_reaction_surface = rate_info.heat_reaction_surface;
                    flame_info_history(i).heat_total = rate_info.heat_total;
                end
            end
        end


        % --- V12 新增: ODE求解相关的辅助函数 ---
        function status = vaporization_output_fcn(obj, t, y, flag, solver, T_p_const)
            % 用于ODE求解器的输出函数，用于记录中间状态
            status = 0; % 默认继续积分
            
            switch flag
                case 'init'
                    % 初始化时不需要额外操作，因为vaporization_ode已经计算并缓存了结果
                    % 初始化实时可视化
                    if obj.visualize_realtime
                        obj.init_realtime_visualization();
                    end
                    
                case ''
                    % 不再重复计算，只记录当前状态
                    % 缓存已在vaporization_ode中更新
                    
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
                            fprintf('t=%.4f ms: 核心温度=%.2f K, 表面温度=%.2f K, 温差=%.2f K\n', t*1000, T_p_const, T_p_current, T_p_current-T_p_const);
                        end
                    end
                    
                    % 实时可视化颗粒状态
                    if obj.visualize_realtime
                        % 限制可视化更新频率，避免过度绘图导致性能问题
                        if t(end) - obj.last_visual_update_time > 1e-4
                            obj.visualize_realtime_state(t(end), y(:,end), T_p_const, solver);
                            obj.last_visual_update_time = t(end);
                        end
                    end
                    
                case 'done'
                    fprintf('积分完成，共缓存%d个速率信息\n', length(obj.vaporization_rate_info_cache));
                    % 在求解完成时保存可视化结果
                    if obj.visualize_realtime && ishandle(obj.realtime_fig)
                        saveas(obj.realtime_fig, 'particle_state_history.png');
                        fprintf('实时可视化结果已保存为: particle_state_history.png\n');
                    end
            end
        end

        function rate_info = create_default_rate_info(obj, state)
            % 创建一个合理的默认rate_info结构
            rate_info = struct();
            rate_info.success = false;
            rate_info.dmdt_mg = -1e-10; % 很小的负值表示缓慢消耗
            rate_info.dmdt_mgo = 1e-11;
            rate_info.dmdt_c = 1e-12;
            rate_info.r_flame = state.r_p * 3; % 默认火焰半径是颗粒的3倍
            rate_info.T_flame = state.T_p + 1000; % 默认火焰温度比颗粒高1000K
            
            % 确保字段名一致性
            rate_info.r_f = rate_info.r_flame;
            rate_info.T_f = rate_info.T_flame;
            
            % 添加其他必要字段
            rate_info.m_dot = abs(rate_info.dmdt_mg);
            rate_info.m_dot_CO = 0;
            rate_info.r_f_to_r_p_ratio = rate_info.r_flame / state.r_p;
        end

        function tempState = build_temp_state_from_vector(obj, y_vec, T_p_const)
             % 辅助函数: 从状态向量构建一个临时的ParticleState对象
             tempState = ParticleState(obj.params);
             tempState.m_mg  = y_vec(1);
             tempState.m_mgo = y_vec(2);
             tempState.m_c   = y_vec(3);
             tempState.r_c   = y_vec(4);
             tempState.r_p   = y_vec(5);
            
                  % 使用状态变量中的温度作为氧化层表面温度，不再使用固定值
     if length(y_vec) >= 6
         tempState.T_p = y_vec(6);
     else
         tempState.T_p = T_p_const;  % 向后兼容，默认使用沸点温度
     end
    
     tempState.oxide_thickness = tempState.r_p - tempState.r_c;
     tempState.T_c = T_p_const;  % 始终保持金属核心温度锁定在沸点
             tempState.melted_fraction = 1.0;
        end
        % --- 结束 V12 新增 ---


        %%%%%  本质即从运行计算的结果中读取想要的参数结果 转化为易于出图的结构形式。
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

            % 新增: 计算各阶段的时间统计
            stage_times = calculate_stage_durations(results);
            results.stage_times = stage_times;
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
        
        % 新增函数: 可视化实时状态
        function visualize_realtime_state(obj, t, y, T_p_const, solver)
            % 从状态向量y创建临时颗粒状态
            tempState = obj.build_temp_state_from_vector(y, T_p_const);
            
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