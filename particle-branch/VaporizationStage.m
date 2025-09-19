classdef VaporizationStage < handle
    % VaporizationStage (气化阶段求解器) 

    properties
        params          % 参数对象
        physicalModel   % 物理模型对象
        thermo          % 热力学数据读取器
    end

    methods
        function obj = VaporizationStage(params, physicalModel)
            % 构造函数
            obj.params = params;
            obj.physicalModel = physicalModel;
            obj.thermo = physicalModel.thermo_reader;
        end

        function rate_info = solve_reaction_rates(obj, pState, previous_solution)
            % --- 主求解函数 ---
            % 将所有依赖项打包到一个结构体中
            bvp_deps.pState = pState;
            bvp_deps.params = obj.params;
            bvp_deps.physicalModel = obj.physicalModel;
            bvp_deps.thermo = obj.thermo;
            porosity = obj.params.material_properties.oxide_porosity;
            % 添加氧化层参数
            bvp_deps.r_c = pState.r_c;
            bvp_deps.oxide_thickness = pState.oxide_thickness;
            
            fprintf('--- 代数方程求解一次 ---\n');
            try
                % 设置产物分配比例参数
                alpha_CO = 0.5;  % CO向内流动的比例，默认50%
                alpha_MgO = 0.5; % MgO向内流动的比例，默认50%
                
                % 如果参数中定义了产物分配比例，则使用参数中的值
                if isfield(obj.params, 'alpha_CO')
                    alpha_CO = obj.params.alpha_CO;
                end
                
                if isfield(obj.params, 'alpha_MgO')
                    alpha_MgO = obj.params.alpha_MgO;
                end
            
                % 扩散求解 结果的质量流率等 - 传递previous_solution作为缓存
                if nargin < 3
                    % 未提供previous_solution，调用无缓存版本
                    rate_info = solve_vaporization_algebraic(bvp_deps, alpha_CO, alpha_MgO);
                else
                    % 提供了previous_solution，传递给求解器
                    rate_info = solve_vaporization_algebraic(bvp_deps, alpha_CO, alpha_MgO, previous_solution);
                end
                
                % 添加热量计算
                % 1. 对流换热
                h_conv = obj.params.k_gas ;
                q_conv = h_conv / pState.r_p * (obj.params.ambient_temperature - pState.T_p);
                q_conv = 0 ;
                rate_info.heat_convection = q_conv * (4 * pi * pState.r_p^2 * porosity);
                
                % 2. 辐射换热
                q_rad = obj.params.emissivity * obj.params.sigma * (obj.params.ambient_temperature^4 - pState.T_p^4);
                rate_info.heat_radiation = q_rad * (4 * pi * pState.r_p^2 * porosity);
                
                % 3. 反应热 (基于Mg消耗率)
                if isfield(rate_info, 'dmdt_mg') && rate_info.dmdt_mg < 0
                    %n_mol = abs(rate_info.dmdt_mg) / obj.params.materials.Mg.molar_mass;
                    reaction_heat_value_flame = abs(obj.params.reaction_heats.flame_reac_H);
                    reaction_heat_value_face = abs(obj.params.reaction_heats.surface_reac_H)
                    % 如果反应释放热量，则为正值（向颗粒提供热量）    
                    reaction_heat_value_flame = abs(rate_info.m_dot_Mg_region1)*reaction_heat_value_flame;
                    reaction_heat_value_face = abs(rate_info.dmdt_Mg_surface_reaction) *reaction_heat_value_face;
                    
                    rate_info.heat_reaction_surface = reaction_heat_value_face;
                    rate_info.heat_reaction_flame = reaction_heat_value_flame;
                    rate_info.heat_reaction_total = reaction_heat_value_flame+reaction_heat_value_face;
                    rate_info.heat_reaction = rate_info.heat_reaction_total;
                    rate_info.heat_to_particle = reaction_heat_value_face + reaction_heat_value_flame*0.3152 ;
                else
                    rate_info.heat_reaction_total = 0;
                    rate_info.heat_to_particle = 0 ;
                    rate_info.heat_to_particle = 0 ;
                    reaction_heat_value_face = 0 ;
                    reaction_heat_value_flame= 0 ;
                    rate_info.heat_reaction = 0;
                    rate_info.heat_reaction_surface = 0;
                    rate_info.heat_reaction_flame = 0;
                    rate_info.heat_reaction_total = 0;
                
                end
                fprintf('  heat_reaction_flame: %.3e \n', rate_info.heat_reaction_flame);
                fprintf('  heat_radiation: %.3e \n', rate_info.heat_radiation);
                %fprintf('  heat_reaction_flame: %.3e \n', rate_info.heat_reaction_flame);
                %fprintf('  heat_reaction_surface: %.3e \n', rate_info.heat_reaction_surface);
                %fprintf('  heat_reaction_total: %.3e \n', rate_info.heat_reaction_total);
                %fprintf('  heat_to_particle: %.3e \n', rate_info.heat_to_particle);

                % 4. 总热量
                rate_info.heat_total_particle = rate_info.heat_convection + rate_info.heat_radiation + rate_info.heat_to_particle;
                  %传递到颗粒表面的热量  考虑不考虑反应呢？
                rate_info.heat_ox = rate_info.heat_total_particle - abs(rate_info.dmdt_mg)/obj.params.materials.Mg.molar_mass...
                 * obj.params.materials.Mg.L_evap_Mg;
                %rate_info.heat_ox =rate_info.heat_convection + rate_info.heat_radiation + reaction_heat_value_flame * 0;
                rate_info.heat_total = rate_info.heat_total_particle ;
                Q_Mg_evap = abs(rate_info.dmdt_mg)/obj.params.materials.Mg.molar_mass* obj.params.materials.Mg.L_evap_Mg ;
                fprintf('  Q_Mg_evap: %.3e \n', Q_Mg_evap);

                fprintf('  reaction_heat_value_flame: %.3e \n', reaction_heat_value_flame);
                fprintf('  ratio: %.3e \n', Q_Mg_evap/reaction_heat_value_flame);
                
            catch ME
                fprintf('--- 代数方程求解器失败 ---\n');
                fprintf('错误信息: %s\n', ME.message);
                
                % 在求解失败时输出关键尺度参数
                r_p = pState.r_p;
                r_inf = 5 * r_p; % 远场边界设为5倍颗粒半径
                fprintf('  > 关键参数: r_p = %.3e m, r_inf = %.3e m (尺度差异: %.1f 倍)\n', ...
                        r_p, r_inf, r_inf/r_p);

                fprintf('详细错误堆栈:\n');
                for i = 1:length(ME.stack)
                    fprintf('  > 文件: %s, 行: %d, 函数: %s\n', ...
                        ME.stack(i).file, ME.stack(i).line, ME.stack(i).name);
                end
                        
                % 重新抛出错误以终止整个计算
                rethrow(ME);
            end
        end
    end
end