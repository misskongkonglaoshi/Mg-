classdef DiffusionControl < handle
    % DiffusionControl 扩散控制模型
    % 处理颗粒在扩散控制阶段的反应过程
    
    properties
        % 反应物和产物的扩散系数
        D_CO2    % CO2的扩散系数
        D_MgO    % MgO的扩散系数
        D_CO2_base % CO2的基础扩散系数（不受氧化层影响时）
        
        % 颗粒和环境参数
        r_p      % 颗粒半径
        T        % 温度
        P        % 压强
        C_CO2_b  % 环境中CO2的浓度
        
        % 反应层参数
        delta    % 产物层厚度
        porosity % 产物层孔隙率
        tau      % 产物层迂曲度
        
        % 颗粒状态
        particle_state
    end
    
    methods
        function obj = DiffusionControl(particle_state, params)
            % 构造函数：初始化扩散控制模型
            % params: 包含所需参数的结构体
            % particle_state: 颗粒状态对象
            
            % 保存颗粒状态
            obj.particle_state = particle_state;
            
            % 初始化扩散系数
            if isfield(params, 'D_CO2')
                obj.D_CO2 = params.D_CO2;
                obj.D_CO2_base = params.D_CO2;
            else
                obj.D_CO2 = 1e-5;  % 默认值 (m²/s)
                obj.D_CO2_base = 1e-5;
            end
            
            if isfield(params, 'D_MgO')
                obj.D_MgO = params.D_MgO;
            else
                obj.D_MgO = 1e-12;  % 默认值 (m²/s)
            end
            
            % 初始化颗粒和环境参数
            if isfield(params, 'initial_diameter')
                obj.r_p = params.initial_diameter / 2;
            else
                obj.r_p = 25e-6;  % 默认值 (m)
            end
            
            if isfield(params, 'ambient_temperature')
                obj.T = params.ambient_temperature;
            else
                obj.T = 1500;  % 默认值 (K)
            end
            
            if isfield(params, 'ambient_pressure')
                obj.P = params.ambient_pressure;
            else
                obj.P = 101325;  % 默认值 (Pa)
            end
            
            % 估算CO2浓度（理想气体近似）
            if isfield(params, 'ambient_gas_composition') && isfield(params.ambient_gas_composition, 'CO2')
                R = 8.314;  % 气体常数 (J/(mol·K))
                obj.C_CO2_b = params.ambient_gas_composition.CO2 * obj.P / (R * obj.T);
            else
                obj.C_CO2_b = 40;  % 默认值 (mol/m³)
            end
            
            % 初始化反应层参数
            if isfield(params, 'initial_oxide_thickness')
                obj.delta = params.initial_oxide_thickness;
            else
                obj.delta = 1e-7;  % 默认值 (m)
            end
            
            obj.porosity = 0.5;  % 默认孔隙率
            obj.tau = 2.0;      % 默认迂曲度
        end
        
        function updateParams(obj, oxide_thickness)
            % 更新与氧化层相关的参数
            % 输入:
            %   oxide_thickness: 氧化层厚度 (m)
            
            % 更新产物层厚度
            obj.delta = oxide_thickness;
            
            % 更新颗粒半径（如果颗粒状态存在）
            if ~isempty(obj.particle_state)
                obj.r_p = obj.particle_state.diameter / 2;
            end
            
            % 更新扩散系数 (简化模型：随氧化层厚度增加而减小)
            if oxide_thickness > 0
                reduction_factor = exp(-oxide_thickness / 1e-6);  % 使用1µm作为特征长度
                obj.D_CO2 = obj.D_CO2_base * reduction_factor;
                
                % 孔隙率随着氧化层增厚而减小
                obj.porosity = max(0.1, 0.5 * reduction_factor);
                
                % 迂曲度随着氧化层增厚而增加
                obj.tau = min(5.0, 2.0 + 3.0 * (1 - reduction_factor));
            else
                obj.D_CO2 = obj.D_CO2_base;
                obj.porosity = 0.5;
                obj.tau = 2.0;
            end
        end
        
        function dX_dt = calculateReactionRate(obj)
            % 计算扩散控制阶段的反应速率
            % 返回：反应转化率的时间导数 dX/dt
            
            % 计算有效扩散系数
            D_eff_CO2 = obj.D_CO2 * obj.porosity / obj.tau;
            
            % 计算扩散通量
            if obj.delta > 0
                J_CO2 = D_eff_CO2 * obj.C_CO2_b / obj.delta;
            else
                J_CO2 = D_eff_CO2 * obj.C_CO2_b / 1e-7;  % 避免除以零
            end
            
            % 计算反应速率（基于扩散通量）
            dX_dt = J_CO2 * 4 * pi * obj.r_p^2 / (obj.r_p^3 * 4/3 * pi);
            
            % 考虑温度和压强的影响
            dX_dt = dX_dt * sqrt(obj.T/298) * (1e5/obj.P);
            
            % 考虑氧化层对反应面积的影响
            if obj.delta > 0
                % 简化模型：氧化层增厚导致有效反应面积减小
                reduction_factor = exp(-obj.delta / 1e-6);  % 使用1µm作为特征长度
                dX_dt = dX_dt * reduction_factor;
            end
        end
        
        function updateProductLayer(obj, dt)
            % 更新产物层参数
            % dt: 时间步长
            
            % 计算产物层增长
            dX = obj.calculateReactionRate() * dt;
            obj.delta = obj.delta + dX * obj.r_p;
            
            % 更新孔隙率（简化模型）
            obj.porosity = max(0.1, obj.porosity - 0.01 * dX);
            
            % 更新迂曲度（简化模型）
            obj.tau = min(5.0, obj.tau + 0.1 * dX);
            
            % 更新扩散系数
            obj.updateParams(obj.delta);
        end
        
        function rate = calculate_diffusion_rate(obj, c_surface, c_bulk)
            % 计算扩散控制下的Mg消耗速率
            % 输入:
            %   c_surface: 表面CO2浓度 (mol/m³)
            %   c_bulk: 环境CO2浓度 (mol/m³)
            % 输出:
            %   rate: Mg消耗速率 (kg/s)
            
            % 如果未提供浓度，使用默认值
            if nargin < 3
                c_bulk = obj.C_CO2_b;
                c_surface = 0;  % 假设表面浓度为零
            end
            
            % 计算浓度梯度
            delta_c = c_bulk - c_surface;
            
            % 计算有效扩散系数
            D_eff = obj.D_CO2 * obj.porosity / obj.tau;
            
            % 计算扩散通量
            if obj.delta > 0
                J = D_eff * delta_c / obj.delta;  % (mol/(m²·s))
            else
                J = D_eff * delta_c / 1e-7;  % 避免除以零
            end
            
            % 计算总扩散速率
            if ~isempty(obj.particle_state)
                surface_area = 4 * pi * (obj.particle_state.diameter/2)^2;
            else
                surface_area = 4 * pi * obj.r_p^2;
            end
            
            % 考虑氧化层对反应面积的影响
            if obj.delta > 0
                % 简化模型：氧化层增厚导致有效反应面积减小
                reduction_factor = exp(-obj.delta / 1e-6);  % 使用1µm作为特征长度
                surface_area = surface_area * reduction_factor;
            end
            
            % 计算总扩散速率 (mol/s)
            mol_rate = J * surface_area;
            
            % 转换为Mg消耗速率 (kg/s)，考虑化学计量比
            % 2Mg + CO2 -> 2MgO + C
            molar_mass_mg = 24.305e-3;  % (kg/mol)
            rate = 2 * mol_rate * molar_mass_mg;
        end
    end
end 