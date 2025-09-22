 classdef Parameters < handle
    % Parameters 仿真参数管理类
    
    properties
        params        % 集中存储所有参数的结构体
        
        % 物理常数
        R_u = 8.314;         % 通用气体常数 (J/(mol·K))
        sigma = 5.67e-8;     % Stefan-Boltzmann常数 (W/(m²·K⁴))
        emissivity = 0.8;           % 发射率
        % 材料物理参数
        materials = struct(...
            'Mg', struct(...
                'density_low', 1740, ...        % 密度 (kg/m^3)   可丰富
                'density_high', 1580, ...  
                'heat_capacity', 1020, ...  % 比热容 (J/(kg·K))
                'latent_heat', 1.97e5, ...  % 熔化潜热 (J/kg)
                'L_evap_Mg', 147.1e3, ...    % 蒸发潜热 (J/mol)
                'molar_mass', 24.305e-3, ...% 摩尔质量 (kg/mol)
                'melting_point', 923, ...   % 熔点 (K)
                'Hf298', 147100, ...        % 气相标准生成焓 (J/mol)
                'ignition_temp', 1363 ...   % 点火温度 (K)
            ), ...
            'MgO', struct(...
                'density', 3650, ...        % 密度 (kg/m^3)
                'heat_capacity', 937, ...   % 比热容 (J/(kg·K))
                'molar_mass', 40.304e-3, ...% 摩尔质量 (kg/mol)
                'thermal_conductivity', 25.1, ...  % 热导率 (W/(m·K))
                'Hf298', -19400, ...        % 气相标准生成焓 (J/mol)
                'formation_enthalpy', -601700, ... % 标准生成焓 (J/mol)
                'diffusivity', 1e-12 ...    % 扩散系数 (m^2/s)
            ), ...
            'C', struct(...
                'density', 1800, ...
                'heat_capacity', 709, ...
                'molar_mass', 12.011e-3, ...
                'Hf298', 0 ...              % 标准生成焓 (J/mol), 稳定单质为0
            ), ...
            'CO', struct(...
                'molar_mass', 28.010e-3, ...
                'Hf298', -110500 ...        % 标准生成焓 (J/mol)
            ), ...
            'CO2', struct(...
                'molar_mass', 44.010e-3, ...
                'Hf298', -393500, ...
                'k_coeffs', [2.5e-3, 5.0e-5] ... % V33: CO2热导率线性多项式 k = a + b*T (a, b)
            ) ...
        );
        cp_mg_gas = 20.79 ; % j/mol/k
        % V33: 固相Mg, MgO, C的比热系数 (J/(kg·K)，线性多项式 Cp = a + b*T)
        solid_cp_coeffs = struct(...
            'Mg', struct('a', 950, 'b', 0.2), ...
            'MgO', struct('a', 900, 'b', 0.1), ...
            'C', struct('a', 600, 'b', 0.5) ...
        );

        % 蒸气压计算参数 
        antoine_coeffs_mg = struct(...
            'A', 7.378, ...
            'B', 1605, ...
            'C', 211 ...
        );

        % V32: 恒定的反应热 (J/kg),针对mg
        reaction_heats = struct(...
            'surface_reac_H', -2.04e7, ... % Mg(g)+CO(g)=MgO(s)+C(s) 非均相
            'flame_reac_H',   -1.92e7  ... % Mg(g)+CO2(g)=MgO(g)+CO(g)  均相
        );
        

        % 氧化层物理属性和破裂模型参数
        material_properties = struct(...
            'alpha_mg', 20e-6, ...          % 镁热膨胀系数 (1/K)
            'alpha_oxide', 9.7e-6, ...        % 氧化镁热膨胀系数 (1/K)
            'E_mg', 45e9, ...               % 镁弹性模量 (Pa)
            'E_oxide', 250e9, ...           % 氧化镁弹性模量 (Pa)
            'nu_mg', 0.35, ...              % 镁泊松比
            'nu_oxide', 0.25, ...           % 氧化镁泊松比
            'oxide_limit_stress', 300e6, ... % 氧化层极限应力 (Pa)
            'oxide_porosity', 0.2, ...      % 氧化层孔隙率
            'oxide_tortuosity', 5, ...    % 氧化层弯曲度因子  2.5
            'k_solid_MgO', 30, ...          % MgO导热系数(W/m/K)
            'k_solid_C', 100, ...           % C导热系数(W/m/K)
            'emissivity_oxide', 0.9 ...     % 氧化层辐射率
        );

       
        % 颗粒物理参数
        initial_diameter = 75.09e-6;  % 初始颗粒直径 (m)   75.09
        initial_temperature = 300;  % 初始颗粒温度 (K)
        initial_oxide_thickness = 1e-7;  % 初始氧化层厚度 (m)，默认为100nm  10e-7
        
        % 环境参数
        ambient_temperature = 1573;  % 环境温度 (K)
        ambient_pressure = 101325;   % 环境压力 (Pa)

        % 求解器配置
        solver_options = struct(...
            'RelTol', 1e-6, ...
            'AbsTol', 1e-8, ...
            'MaxStep', 1e-5 ...
        );

        % 气体组成 (摩尔分数)
        ambient_gas_composition = struct(...
            'CO2', 1.0, ...  % 纯CO2环境
            'O2', 0.0, ...
            'N2', 0.0 ...
        );

        % 气体物性参数
        cp_gas  = 1255.2;               % 气体比热估计值 [J/(kg·K)] 1400   2500  wen 1255.2
        k_gas = 0.0280241;            % 气体导热系数估计值 [W/(m·K)] 0.102682     0.3278
        gas_density = 0.2;           % 环境气体密度 (kg/m³)  1.1 
        gas_viscosity = 3e-5;        % 环境气体黏度 (Pa·s)
        rho_D_gas = 2.24501e-5;       % 质量扩散系数参数 kg/(m·s)  2.8688e-5    8.5 e-5
        %D_gas = 1.63e-5 ;
        %D_ox  = D_gas * material_properties.oxide_porosity / material_properties.oxide_tortuosity;
        %rho_D_ox = D_ox * gas_density;
        k_eff = 25.1 ;               % Mgo 导热
        k_ox = 20.1 ;
        k_solid = 25.1 ;
        % 气体物理属性 (温度依赖)
        gas_properties = struct(...
            'O2_k', 0.03, ... % O2热导率常数 (W/(m·K))
            'N2_k', 0.027 ... % N2热导率常数 (W/(m·K))
        );
    
        % 仿真控制参数
        time_step = 1e-6;           % 时间步长 (s)
        total_time = 0.5;           % 总仿真时间 (s)
        t_combustion = 0.1;         % 气相燃烧求解时长 (s)
        output_interval = 1000;     % 输出间隔（步数）
        combustionfinish_ratio = 0.01; % 仿真结束判据（mg的质量）
        flam_thickness = 2.5e-5;    % 气膜厚度 (m)
        T_reaction_begin = 923;     % 反应开始温度 (K)
        
        % 气相燃烧产物分配比例参数
        alpha_CO = 0.5;             % CO向内流动的比例（0-1之间，默认0.5表示均匀分配）
        alpha_MgO = 0.5;            % MgO向内流动的比例（0-1之间，默认0.5表示均匀分配）
        
        % 氧化层破裂控制参数
        CO2_depletion_factor = 0.2;  % CO2浓度随温度衰减因子
        oxide_break_parameters = struct(...
            'relaxation_factor', 0.95, ... % 应力松弛因子
            'recovery_time', 1e-5, ...     % 破裂后恢复期 (s)
            'reaction_rate_boost', 5.0 ... % 破裂后反应速率增强因子
        );

        % 表面异相反应 参数
        reaction_pre_exponential = 1.376e4; % 指前因子
        reaction_activation_energy = 1.324e5; % 活化能 [J/mol]
        
        reaction_heat_face = -1.66e7;     % 表面异相反应热 J/kg
        H_dep_MgO = 139e4 ; % J/kg
        % 组分质量流率相对于总质量流率
        alfa_mg = 1.5;                   % Mg比例
        alfa_co = -0.25;                 % CO比例
        alfa_co2 = 0.0;                  % CO2比例
        alfa_mgo = -0.25;                % MgO比例
        



        k_mg= 156.0;
        k_mgo = 25.1;
        k_c = 23.8;

        r_inf_particle = 10 ;
        % 控制标志
        visualization = false;            % 是否可视化结果
        visualize_residual = false;       % 是否可视化残差
        recordResults = true;            % 是否记录结果
        debug = true;                   % 是否开启调试模式
        current_stage = '';              % 当前阶段标记
        consider_oxide_break = false;     % 是否考虑氧化层破裂
        plot_radial_distributions = true ;
        visualize_realtime = true;      % 是否开启实时可视化
        visualize_realtime_interval = 1e-4; % 实时可视化更新时间间隔(s)
        use_solution_cache = false;       % 是否启用代数方程求解缓存机制
        
        % 气相燃烧求解控制参数
        use_fixed_timestep = true;      % 是否使用固定步长求解
        vaporization_fixed_timestep = 1e-4;  % 固定时间步长(s)
    end
    
    methods
        function obj = Parameters()
            % 构造参数对象
            % 将类属性复制到params结构体中
            try
                % 创建参数结构体并从类属性复制
                fieldNames = properties(obj);
                params = struct();
                
                for i = 1:length(fieldNames)
                    fieldName = fieldNames{i};
                    if ~strcmp(fieldName, 'params') % 跳过params自身
                        params.(fieldName) = obj.(fieldName);
                    end
                end
                
                % 存储参数
                obj.params = params;
                
            catch ME
                fprintf('参数初始化错误: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        function validate(obj)
            % 验证参数的有效性
            assert(obj.initial_diameter > 0, '初始直径必须大于0');
            assert(obj.initial_temperature > 0, '初始温度必须大于0');
            assert(obj.materials.Mg.melting_point > 0, '熔点必须大于0');
            assert(obj.materials.Mg.ignition_temp > obj.materials.Mg.melting_point, ...
                '点火温度必须高于熔点');
            assert(obj.ambient_temperature > 0, '环境温度必须大于0');
            assert(obj.ambient_pressure > 0, '环境压力必须大于0');
            assert(obj.emissivity > 0 && obj.emissivity <= 1, '发射率必须在0到1之间');
            assert(obj.time_step > 0, '时间步长必须大于0');
            assert(obj.total_time > 0, '总仿真时间必须大于0');
            assert(obj.output_interval > 0, '输出间隔必须大于0');
            assert(obj.initial_oxide_thickness >= 0, '初始氧化层厚度必须大于等于0');
        end
    end
end 