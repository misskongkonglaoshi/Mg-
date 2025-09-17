classdef PhysicalModel < handle
    % PhysicalModel 物理模型计算中心
    % 集中管理所有与物理、化学、热力学相关的计算。
    
    properties
        params          % 参数对象
        thermo_reader   % 热力学数据读取器
    end
    
    methods
        function obj = PhysicalModel(params, thermo_reader)
            % 构造函数
            %
            % 输入:
            %   params: 参数对象
            %   thermo_reader: ThermoReader 实例
            if nargin > 0
            obj.params = params;
            obj.thermo_reader = thermo_reader;
            end
        end
        
        function Re = calculate_reynolds_number(obj, particleState)
            % 计算颗粒雷诺数
            rho_gas = obj.params.gas_density;
            mu_gas = obj.params.gas_viscosity;
            v_rel = obj.params.relative_velocity;  % 颗粒相对于气体的速度

            Re = rho_gas * v_rel * 2 * particleState.r_p / mu_gas;
        end
        
        function C_CO2 = calculate_ambient_CO2_concentration(obj)
            % 计算环境中的CO2浓度
            % 使用理想气体定律: c = p/(R*T)
            R_u = obj.params.R_u;
            T_amb = obj.params.ambient_temperature;
            p_amb = obj.params.ambient_pressure;
            
            % CO2摩尔分数
            x_CO2 = obj.params.ambient_gas_composition.CO2;
            
            % CO2分压
            p_CO2 = p_amb * x_CO2;
            
            % CO2摩尔浓度 [mol/m³]
            C_CO2_mol = p_CO2 / (R_u * T_amb);
            
            % 转换为质量浓度 [kg/m³]
            C_CO2 = C_CO2_mol * obj.params.materials.CO2.molar_mass;
        end

        function C_CO2_surf = calculate_surface_CO2_concentration(obj, T_p)
            % 计算颗粒表面的CO2浓度（考虑温度效应）
            C_CO2_inf = obj.calculate_ambient_CO2_concentration();
            T_ref = obj.params.ambient_temperature;
            
            % 随温度指数衰减的经验关系
            surface_depletion_factor = 0.5; % 默认表面消耗因子
            if isfield(obj.params, 'surface_depletion_factor')
                surface_depletion_factor = obj.params.surface_depletion_factor;
            end
            
            C_CO2_surf = C_CO2_inf * exp(-surface_depletion_factor * (T_p - T_ref) / T_ref);
            C_CO2_surf = max(C_CO2_surf, 0.01 * C_CO2_inf);  % 设置最小值避免为零
        end

        function k_reaction = calculate_reaction_rate_constant(obj, T_p)
            % 计算表面反应速率常数（阿伦尼乌斯方程）   (氧化阶段mg二氧化碳异相反应)
            % 默认反应参数
            A_pre = 1.376e4;  % 指前因子
            E_a = 132400;  % 活化能 [J/mol]
        
            R = obj.params.R_u;  % 通用气体常数
            k_reaction = A_pre * exp(-E_a / (R * T_p));
        end
        
        function rates = calculate_oxidation_rates(obj, particleState)
            % 计算表面氧化反应的速率
            % 基于CO2在氧化层中的扩散和反应动力学
            
            T_p = particleState.T_p;
            r_p = particleState.r_p;
            r_c = particleState.r_c;
            oxide_thickness = r_p - r_c;

            % 获取材料参数
            materials = obj.params.materials;
            rho_mg = materials.Mg.density_low;
            rho_mgo = materials.MgO.density;
            rho_c = materials.C.density;
            
            % 1. 计算CO2扩散到颗粒核心的通量
            % 使用努塞尔特数和相关传质系数
            %Re = obj.calculate_reynolds_number(particleState);
            %Pr = obj.params.Pr;  % 确保参数中有普朗特数
            %Nu = 2.0 + 0.6 * Re^0.5 * Pr^(1/3);  % 球形颗粒努塞尔特数
            
            % CO2在环境和氧化层中的扩散系数 (温度相关)
            D_CO2_ambient = obj.params.rho_D_gas / obj.params.gas_density;  % 从参数中获取
            
            % 氧化层中的扩散系数会受到多孔介质结构影响，通常比自由气体中小
            porosity = 0.2;  % 氧化层孔隙率估计值
            tortuosity = 2.5; % 弯曲度因子估计值
            D_CO2_oxide = D_CO2_ambient * porosity / tortuosity;
            
            % 环境CO2浓度和颗粒表面CO2浓度
            C_CO2_inf = obj.calculate_ambient_CO2_concentration();
            C_CO2_surf = obj.calculate_surface_CO2_concentration(T_p);
            
            % 使用壳层扩散模型计算通量 (考虑球形几何)
            % 对于球壳扩散，通量 = D_eff * (C_outer - C_inner) * 4*pi*r_c*r_p/(r_p-r_c)
            geo_factor = 4 * pi * r_c * r_p / (r_p - r_c + eps);
            N_CO2 = D_CO2_oxide * (C_CO2_surf - 0) * geo_factor;  % 内表面浓度为0（全部消耗）
            
            % 2. 考虑化学反应动力学限制
            % 表面反应: Mg + CO2 → MgO + C
            k_reaction = obj.calculate_reaction_rate_constant(T_p);
            reaction_area = 4 * pi * r_c^2 * particleState.reaction_area_factor;  % 有效反应面积
            N_CO2_reaction = k_reaction * C_CO2_surf * reaction_area;  % 反应限制的CO2消耗速率
            
            % 取扩散和反应的最小值（控制步骤）
            N_CO2_effective = min(N_CO2, N_CO2_reaction);
            
            % 3. 计算质量变化率
            % 化学计量比: 1 mol Mg + 1 mol CO2 → 1 mol MgO + 1 mol C
            dn_CO2_dt = N_CO2_effective;  % CO2消耗的摩尔速率 [mol/s]
            dn_Mg_dt = -dn_CO2_dt;        % Mg消耗 (负值)
            dn_MgO_dt = -dn_Mg_dt;        % MgO生成 (正值)
            dn_C_dt = -dn_Mg_dt;          % C生成 (正值)
            
            % 转换为质量变化率 [kg/s]
            rates.dmg_dt = dn_Mg_dt * materials.Mg.molar_mass;
            rates.dmgo_dt = dn_MgO_dt * materials.MgO.molar_mass;
            rates.dc_dt = dn_C_dt * materials.C.molar_mass;
            
            % 4. 计算几何变化率
            % 核心半径变化（基于Mg消耗）
            dV_mg_dt = rates.dmg_dt / rho_mg;
            rates.drc_dt = dV_mg_dt / (4 * pi * r_c^2 + eps);
            
            % 颗粒外半径变化（基于MgO和C沉积）
            dV_mgo_dt = rates.dmgo_dt / rho_mgo;
            dV_c_dt = rates.dc_dt / rho_c;
            dV_product_deposition = dV_mgo_dt + dV_c_dt;
            rates.drp_dt = dV_product_deposition / (4 * pi * r_p^2 + eps);
            
            % 5. 
            
            rates.reaction_heat = rates.dmg_dt  * obj.params.reaction_heat_face;  % 反应热 [W]
        end
        
        function cp_total = calculate_total_heat_capacity(obj, particleState)
            % 计算颗粒的总热容 (J/K)
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            % 输出:
            %   cp_total: 总热容 (J/K)
            
            % 计算各组分的摩尔数
            n_mg = particleState.m_mg / obj.params.materials.Mg.molar_mass;
            n_mgo = particleState.m_mgo / obj.params.materials.MgO.molar_mass;
            n_total = n_mg + n_mgo;
            
            if n_total > 1e-20
                % 计算摩尔分数
                x_mg = n_mg / n_total;
                x_mgo = n_mgo / n_total;
                
                % 计算当前温度下的摩尔比热容
                cp_mg_mol = obj.calc_cp_solid_mg(particleState.T_p);
                cp_mgo_mol = obj.calc_cp_solid_mgo(particleState.T_p);
                
                % 按摩尔分数加权计算总的摩尔比热容
                cp_mol = x_mg * cp_mg_mol + x_mgo * cp_mgo_mol;
                
                % 转换为总热容
                cp_total = cp_mol * n_total;
            else
                cp_total = 0;
            end
            
            % 添加碳的贡献
            cp_total = cp_total + particleState.m_c * obj.calc_cp_solid_c(particleState.T_p);
        end
        
        function cp = calc_cp_solid_mg(~, T)
            % 计算固态镁的摩尔比热容 (J/(mol·K))
            % 使用线性多项式 Cp = a + b*T
            t = T /1000;
            if T < 923
                a = 26.541;
                b = -1.533;
                c = 8.062 ;
                d = 0.572 ;
                e = -0.174 ;
                cp = a+ b * t + c * t^2 + d * t ^ 3 + e/t^2 ;
            else
                a = 34.309 ;
                b = -7.471e-10 ;
                c = 6.146e-10 ;
                d = -1.598e-10 ;
                e = -1.152e-11 ;
                cp = a+ b * t + c * t^2 + d * t ^ 3 + e/t^2 ;
            end
        %    a = 24.869;  % J/(mol·K)
        %    b = 0.00313; % J/(mol·K²)
        %    cp = a + b * T;
        end
        
        function cp = calc_cp_solid_mgo(~, T)
            % 计算氧化镁的摩尔比热容 (J/(mol·K))
            % 使用线性多项式 Cp = a + b*T
            t = T /1000;
            a = 47.26 ;
            b = 5.682 ;
            c = -0.873 ;
            d = 0.104 ;
            e = -1.054 ;
            cp = a+ b * t + c * t^2 + d * t ^ 3 + e/t^2 ;
        
        %    a = 42.59;   % J/(mol·K)
        %    b = 0.00735; % J/(mol·K²)
        %    cp = a + b * T;
        end
        
        function cp = calc_cp_solid_c(~, T)
            % 计算固态碳的质量比热容 (J/(mol·K))
            cp = 8.6186 ;
        end
        
        function cp = calculate_specific_heat(obj, species, T)
            % 使用thermo_reader计算特定温度下的摩尔比热容
            %
            % 输入:
            %   species: 物质名称
            %   T: 温度 (K)
            % 输出:
            %   cp: 摩尔比热容 (J/(mol·K))
            
            % 根据不同物种调用对应的比热计算函数
            switch species
                case 'Mg'
                    cp = obj.calc_cp_solid_mg(T);
                case 'MgO'
                    cp = obj.calc_cp_solid_mgo(T);
                case 'C'
                    cp = obj.calc_cp_solid_c(T) ;
                otherwise
                    % 对于气相物种，尝试使用thermo_reader
                    if ~isempty(obj.thermo_reader) && obj.thermo_reader.has_species(species)
                        % 使用NASA 7系数多项式计算
                        cp = obj.thermo_reader.calculate_Cp(species, T);
                    else
                        % 回退到参数中的常数值
                        cp = obj.params.materials.(species).heat_capacity * obj.params.materials.(species).molar_mass;
                    end
            end
            
            % 防御性修复: 确保返回的是一个标量
            if ~isscalar(cp)
                cp = sum(cp);
            end
        end
        
        function cp = calc_cp_gas(obj, species, T)
            % 计算气体的摩尔比热容 (J/(mol·K))
            % 使用NASA 7系数多项式
            
            if ~isempty(obj.thermo_reader) && obj.thermo_reader.has_species(species)
                cp = obj.thermo_reader.calculate_Cp(species, T);
            else
                % 如果没有数据，使用默认值
                switch species
                    case 'CO2'
                        % 简化的CO2比热多项式
                        cp = 22.26 + 5.981e-2*T - 3.501e-5*T^2 + 7.469e-9*T^3;
                        cp = 60.038 ;
                    case 'CO'
                        % 简化的CO比热多项式
                        cp = 25.56 + 6.096e-3*T + 4.054e-6*T^2 - 2.671e-9*T^3;
                    case 'O2'
                        % 简化的O2比热多项式
                        cp = 25.48 + 1.520e-2*T - 7.155e-6*T^2 + 1.312e-9*T^3;
                    case 'N2'
                        % 简化的N2比热多项式
                        cp = 28.58 - 3.330e-3*T + 1.035e-5*T^2 - 3.729e-9*T^3;
                    otherwise
                        % 默认值
                        cp = 29.1; % 近似值 (J/(mol·K))
                end
            end
        end
        
        function Q_total = calculate_heat_flux(obj, particleState)
            % 计算颗粒与环境之间的总热通量 (W)
            %
            % 输入:
            %   particleState: 当前的颗粒状态对象
            % 输出:
            %   Q_total: 总热通量 (W), >0表示吸热
      
            T_p = particleState.T_p;
            A_p = 4 * pi * particleState.r_p^2;

            % 对流换热
            %%%%   
            h_conv = obj.params.k_gas / particleState.r_p;
            q_conv = h_conv * (obj.params.ambient_temperature - T_p);
            
            % 辐射换热
            q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - T_p^4);
            
            % 总热通量 (W)
            Q_total = (q_conv + q_rad ) * A_p;
        end
        

        function h_conv = calculate_h_conv(obj,particleState)
            T = particleState.Tp ;

            k_mix = calculate_h_conv();
            
            h_conv = 0.5 ;
        end


        function props = get_gas_properties(obj, T_gas)
            % 获取在特定温度下的气体属性
            props = obj.params.gas_properties;
            
            % 计算温度依赖的热导率
            props.k_gas = obj.calc_k_gas_mixture(T_gas);
            
            % 计算温度依赖的比热容
            props.Cp_gas = obj.calc_cp_gas_mixture(T_gas);
        end
        
        function k_mix = calc_k_gas_mixture(obj, T)
            % 计算混合气体的热导率 (W/(m·K))
            gas_comp = obj.params.ambient_gas_composition;
            gas_fields = fieldnames(gas_comp);
            
            k_mix = 0;
            total_fraction = 0;
            
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = gas_comp.(gas_name);
                
                if mole_fraction > 0
                    k_species = obj.calc_k_gas_species(gas_name, T);
                    k_mix = k_mix + mole_fraction * k_species;
                    total_fraction = total_fraction + mole_fraction;
                end
            end
            
            if total_fraction > 1e-9
                k_mix = k_mix / total_fraction;
            else
                k_mix = 0.026; % 默认值 (W/(m·K))
            end
        end
        
        function k = calc_k_gas_species(~, species, T)
            % 计算单个气体组分的热导率 (W/(m·K))
            switch species
                case 'CO2'
                    % CO2热导率多项式系数
                    a = -0.01183;
                    b = 1.0174e-4;
                    c = -2.2242e-8;
                    k = a + b * T + c * T ^2;
                case 'O2'
                    % O2热导率多项式系数
                    a = 0.0074;
                    b = 7.0e-5;
                    k = a + b * T;
                case 'N2'
                    % N2热导率多项式系数
                    a = 0.0066;
                    b = 6.5e-5;
                    k = a + b * T;
                case 'CO'
                    % CO热导率多项式系数
                    a = 0.0059;
                    b = 6.3e-5;
                    k = a + b * T;
                case 'Mg'
                    % 气态Mg热导率（简化）
                    k = 0.01;
                otherwise
                    k = 0.026; % 默认值
            end
        end
        
        function cp_mix = calc_cp_gas_mixture(obj, T)
            % 计算混合气体的质量比热容 (J/(kg·K))
            gas_comp = obj.params.ambient_gas_composition;
            gas_fields = fieldnames(gas_comp);
            
            cp_mol_mix = 0;
            total_fraction = 0;
            
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = gas_comp.(gas_name);
                
                if mole_fraction > 0
                    cp_mol_species = obj.calc_cp_gas(gas_name, T);
                    cp_mol_mix = cp_mol_mix + mole_fraction * cp_mol_species;
                    total_fraction = total_fraction + mole_fraction;
                end
            end
            
            if total_fraction > 1e-9
                cp_mol_mix = cp_mol_mix / total_fraction;
                
                % 转换为质量比热容
                M_mix = obj.get_ambient_mixture_molar_mass();
                cp_mix = cp_mol_mix / M_mix;
            else
                cp_mix = 1000; % 默认值 (J/(kg·K))
            end
        end
        
         % 根据环境气体组分计算加权平均的摩尔质量 (kg/mol)
        function M_mix = get_ambient_mixture_molar_mass(obj)
           
            M_mix = 0;
            total_fraction = 0;
            
            gas_fields = fieldnames(obj.params.ambient_gas_composition);
            for i = 1:length(gas_fields)
                gas_name = gas_fields{i};
                mole_fraction = obj.params.ambient_gas_composition.(gas_name);
                
                if mole_fraction > 0
                    if isfield(obj.params.materials, gas_name)
                        molar_mass_gas = obj.params.materials.(gas_name).molar_mass;
                        M_mix = M_mix + mole_fraction * molar_mass_gas;
                        total_fraction = total_fraction + mole_fraction;
                    else
                        warning('PhysicalModel:gasNotFound', ...
                            '在 materials 结构体中未找到气体 "%s" 的数据。', gas_name);
                    end
                end
            end
            
            if total_fraction > 1e-9 % 避免除零
                M_mix = M_mix / total_fraction;
            else
                warning('PhysicalModel:noGasComposition', ...
                    '环境气体组分总和为零，无法计算平均摩尔质量。');
                M_mix = 28.97e-3; % 回退到空气的近似值
            end
        end
        
        function H = get_enthalpy_from_state(obj, particleState)
            % 根据颗粒状态（温度和熔化分数）计算其总焓 (J)
            % 参考点: 固态物质在0K时焓为0
            
            T = particleState.T_p;
            X = particleState.melted_fraction;
            T_melt = obj.params.materials.Mg.melting_point;
            
            % 我们需要一个在0K到T_melt之间的平均总热容
            % 为简单起见, 我们使用在初始温度下的热容作为整个固相的常数热容
            tempState_solid = particleState.copy();
            tempState_solid.T_p = obj.params.initial_temperature;
            Cp_total_solid = obj.calculate_total_heat_capacity(tempState_solid);

            % 熔点时固相的焓
            H_solid_at_melt = Cp_total_solid * T_melt;
            
            % 熔化潜热
            L_m_total = particleState.m_mg * obj.params.materials.Mg.latent_heat;

            if X < 1e-9 % 固相
                H = Cp_total_solid * T;
            elseif X >= 1e-9 && X < 1 % 熔化中
                H = H_solid_at_melt + X * L_m_total;
            else % 液相
                tempState_liquid = particleState.copy();
                tempState_liquid.T_p = T_melt; % 使用熔点温度估算液相热容
                Cp_total_liquid = obj.calculate_total_heat_capacity(tempState_liquid);
                
                H_liquid_at_melt = H_solid_at_melt + L_m_total;
                H = H_liquid_at_melt + Cp_total_liquid * (T - T_melt);
            end
        end

        function pStateOut = get_state_from_enthalpy(obj, H, pStateIn)
            % 根据总焓 H 反算颗粒的状态（温度和熔化分数）
            pStateOut = pStateIn.copy();

            % 使用与get_enthalpy_from_state一致的逻辑和近似
            tempState_solid = pStateIn.copy();
            tempState_solid.T_p = obj.params.initial_temperature;
            Cp_total_solid = obj.calculate_total_heat_capacity(tempState_solid);
            
            T_melt = obj.params.materials.Mg.melting_point;
            
            H_solid_at_melt = Cp_total_solid * T_melt;
            L_m_total = pStateIn.m_mg * obj.params.materials.Mg.latent_heat;
            
            if L_m_total < 1e-12 % 防止除零
                H_liquid_at_melt = H_solid_at_melt;
            else
                H_liquid_at_melt = H_solid_at_melt + L_m_total;
            end

            if H < H_solid_at_melt
                pStateOut.T_p = H / Cp_total_solid;
                pStateOut.melted_fraction = 0;
            elseif H >= H_solid_at_melt && H <= H_liquid_at_melt
                pStateOut.T_p = T_melt;
                if L_m_total > 1e-12
                    pStateOut.melted_fraction = (H - H_solid_at_melt) / L_m_total;
                else
                    pStateOut.melted_fraction = 1.0;
                end
            else
                tempState_liquid = pStateIn.copy();
                tempState_liquid.T_p = T_melt;
                Cp_total_liquid = obj.calculate_total_heat_capacity(tempState_liquid);
                
                pStateOut.T_p = T_melt + (H - H_liquid_at_melt) / Cp_total_liquid;
                pStateOut.melted_fraction = 1.0;
            end
        end

        function pStateOut = get_state_from_enthalpy_and_mass(obj, H, pStateIn)
            % 根据总焓H和已知质量反算颗粒的状态（温度和熔化分数）
            % 在质量可能变化的情况下使用此函数
            pStateOut = pStateIn.copy();

            % 使用与get_enthalpy_from_state一致的逻辑和近似
            tempState_solid = pStateIn.copy();
            tempState_solid.T_p = obj.params.initial_temperature;
            Cp_total_solid = obj.calculate_total_heat_capacity(tempState_solid);
            
            T_melt = obj.params.materials.Mg.melting_point;
            
            H_solid_at_melt = Cp_total_solid * T_melt;
            L_m_total = pStateIn.m_mg * obj.params.materials.Mg.latent_heat;
            
            if L_m_total < 1e-12 % 防止除零
                H_liquid_at_melt = H_solid_at_melt;
            else
                H_liquid_at_melt = H_solid_at_melt + L_m_total;
            end

            if H < H_solid_at_melt
                pStateOut.T_p = H / Cp_total_solid;
                pStateOut.melted_fraction = 0;
            elseif H >= H_solid_at_melt && H <= H_liquid_at_melt
                pStateOut.T_p = T_melt;
                if L_m_total > 1e-12
                    pStateOut.melted_fraction = (H - H_solid_at_melt) / L_m_total;
                else
                    pStateOut.melted_fraction = 1.0;
                end
            else
                tempState_liquid = pStateIn.copy();
                tempState_liquid.T_p = T_melt;
                Cp_total_liquid = obj.calculate_total_heat_capacity(tempState_liquid);
                
                pStateOut.T_p = T_melt + (H - H_liquid_at_melt) / Cp_total_liquid;
                pStateOut.melted_fraction = 1.0;
            end
        end


        function dHdt = calculate_heat_flux_with_oxidation(obj, particleState, oxidation_rates)
            % 计算包含氧化反应热的总热流
            
            % 原有的热流（对流、辐射等）
            dHdt_base = obj.calculate_heat_flux(particleState);
            
            % 氧化反应热贡献
            dHdt_oxidation = oxidation_rates.reaction_heat;
            
            % 总热流
            dHdt = dHdt_base + dHdt_oxidation;
            
            if oxidation_rates.dmg_dt < 0  % 确实在发生氧化
                %fprintf(' 考虑氧化反应后:  氧化反应热: %.2e W, 其余基础热流（仅来自于环境的导热和辐射）: %.2e W\n', dHdt_oxidation, dHdt_base);
            end
        end
    end
end 