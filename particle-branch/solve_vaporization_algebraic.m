function rate_info = solve_vaporization_algebraic(bvp_deps, alpha_CO, alpha_MgO, previous_solution)
    % SOLVE_VAPORIZATION_ALGEBRAIC - 使用区域特定解析解求解气相燃烧问题
    % 参数:
    %   bvp_deps - 求解依赖参数
    %   alpha_CO - CO向内流动的比例，默认0.5
    %   alpha_MgO - MgO向内流动的比例，默认50%
    %   previous_solution - 可选参数，上一时间步的解，用于缓存机制
    %   猜测值针对颗粒表面到氧化层展开、
    % =========================================================================
    % 步骤1: 提取颗粒状态参数、物理计算调用函数
    % =========================================================================
        pState  = bvp_deps.pState;
        params = bvp_deps.params;
        physicalModel = bvp_deps.physicalModel;
    
        porosity = params.material_properties.oxide_porosity;
        tortuosity = params.material_properties.oxide_tortuosity;
        rho_D_gas = params.rho_D_gas;

        
        % 气体物性的无量纲参数     co2  60 j/mol/k
        % 气体比热估计值 [J/(kg·K)]
        % co2  2000K
        % 气体导热系数估计值 [W/(m·K)]
        cp_gas = params.cp_gas;
        k_gas = params.k_gas;
        k_solid = params.k_solid;
        %k_ox = k_gas* porosity + (1-porosity) * k_solid;
        k_ox = k_gas * porosity ;
        k_ox_gas = k_gas * porosity;
        le = cp_gas / k_gas * rho_D_gas ;
        D_ox_gas = rho_D_gas * porosity / tortuosity;
        ratio_k_cp_ox = D_ox_gas ;
        ratio_k_cp_gas = rho_D_gas;

        r_p = pState.r_p;            % 颗粒半径（氧化层外表面）
        r_c = pState.r_c;            % 核心半径（金属核心表面）
        T_p = pState.T_p;            % 氧化层表面温度
        %T_p = 1600 ;
        T_c = pState.T_c;            % 金属核心温度 (新增)
        T_amb = params.ambient_temperature; % 环境温度
        r_inf = params.r_inf_particle * r_p;      % 使用参数中设置的无量纲远场边界

        m_mg = pState.m_mg;
        % 固定物性参数
        Le = 1.0;                    % 刘易斯数(简化值)
        mw = params.materials;       % 物质的摩尔质量
    
    fprintf('气相扩散燃烧求解开始：\n');
    fprintf('  > 颗粒半径: %.3e um, 外边界: %.3e um, 金属质量: %.5e kg\n', r_p*1e6, r_inf*1e6, m_mg);
    fprintf('  > 核心半径: %.3e um, 氧化层厚度: %.3e um\n', r_c*1e6, (r_p - r_c)*1e6);
    fprintf('  > 核心温度: %.2f K, 氧化层表面温度: %.2f K, 环境温度: %.2f K\n', T_c, T_p, T_amb);
    
    % =========================================================================
    % 步骤2: 生成无量纲初始猜测或使用缓存的解
    % =========================================================================
        % 检查是否启用缓存机制
        use_cache = true;
        if isfield(params, 'use_solution_cache')
            use_cache = params.use_solution_cache;
        end
        
        % 定义特征热量尺度 - 用于所有热量项的无量纲化
         
        k_0 = k_gas;  % 参考导热系数
        Q_0 = 4 * pi * r_c * k_0 * T_c;  % 特征热量尺度，基于金属核心 [W]
        
        % 火焰位置初始猜测 (无量纲)
        r_f_nd_guess = 1.5 * r_p / r_c;  % 无量纲火焰半径 r_f/r_c = 3.0*r_p/r_c
       
        
        % 初始pem数猜测
        Q_rad = params.emissivity * params.sigma * (4 * pi * r_p^2) * (T_amb^4 - T_p^4);
        Q_rad_f = params.emissivity * params.sigma * (4 * pi * r_p^2) * (2200^4 - T_p^4);
        

        % q_conv = params.h_conv * (T_amb -T_p);
        %fprintf('  > Q_rad = %.6e\n', Q_rad );
        q_conv = k_gas / r_p * (T_amb -T_p);
        %fprintf('  > Q_rad_f = %.6e\n', Q_rad_f);
        heat_convection = q_conv * (4 * pi * r_p^2);
        %heat_convection = 0 ;        
        Q_total = Q_rad + Q_rad_f + heat_convection;
        L_v = params.materials.Mg.L_evap_Mg / mw.Mg.molar_mass;
        %m_dot_guess = max(Q_rad / L_v, 1e-10);
        m_dot_guess = Q_total / L_v  *1.5 ;
        Pe_m_guess_1 = m_dot_guess / (4 * pi * r_c * rho_D_gas);
        fprintf('m_dot_guess  %.6e \n',m_dot_guess);
        fprintf('r_c  %.6e \n',r_c);
        fprintf('D_ox_gas  %.6e \n',D_ox_gas);
        %Pe_m_guess = m_dot_guess / (4 * pi * r_c * rho_D_gas );
    
    fprintf('  > 初始颗粒半径猜测: r_p = %.2f μm, r_c = %.2f μm\n', r_p*1e6, r_c*1e6);
    fprintf('  > 初始火焰半径猜测: r_f/r_c = %.2f\n', r_f_nd_guess);
    fprintf('  > 初始佩克莱数猜测: %.3e (m_dot = %.3e kg/s)\n', Pe_m_guess_1, m_dot_guess);
    
    % 估计初始总质量流率，用于无量纲化
    m_dot_total_region1 = m_dot_guess ;
    %m_dot_total = Pe_m_guess * (4 * pi * r_c * rho_D_gas);m_dot_guess
    
    r_f_guess = r_f_nd_guess * r_c ;  % 物理火焰半径基于核心半径

    % 使用缓存的解或生成新的初始猜测
    if use_cache && exist('previous_solution', 'var') && ~isempty(previous_solution)
        X0 = previous_solution;
        Pe_m_region1 = previous_solution(2);
        m_dot_region1_total = Pe_m_region1 * 4 * pi * r_c * rho_D_gas;
        m_dot_total_region1 = m_dot_region1_total;
        fprintf('  > 使用上一时间步的解作为初始猜测\n');
        r_f_nd= X0(1) ;
        fprintf('  > 上一步无量纲火焰半径: %.3e um\n', r_f_nd);
        fprintf('  > 这一步计算得到火焰位置距离颗粒表面距离: %.3e um\n', r_f_nd * r_c - r_p);
        %fprintf('  > 上一步Mg质量流率: %.3e \n', X0(12)*m_dot_total_region1);
    else
        % 生成完整初始猜测向量 (无量纲化)   所求向量内容
        X0 = generate_initial_guess(m_dot_total_region1,r_f_guess, Pe_m_guess_1, r_p, r_c, r_inf, T_p,T_c, T_amb, Le, mw, physicalModel, params,cp_gas...
        ,k_gas,rho_D_gas,porosity,tortuosity,D_ox_gas,ratio_k_cp_ox,ratio_k_cp_gas,k_ox,k_ox_gas);
        fprintf('  > 生成全新初始猜测\n');
    end
    
    % 分析初始猜测合理性
    %  analyze_initial_guess(X0, r_p);
    
    % =========================================================================
    % 步骤3: 逐步求解策略
    % =========================================================================
    try
        % --- 阶段1: 求解高度简化模型 ---
            fprintf('  > 阶段1/3: 求解高度简化模型...\n');   
            % 定义简化模型的求解环境
            solve_env1.simplified_mode = true;
            solve_env1.simplified_reaction = true;
            solve_env1.relaxation_level = 0.1;
            solve_env1.Pe_m_factor = 0.2;
            solve_env1.Le = Le;
            solve_env1.bvp_deps = bvp_deps;
            solve_env1.last_X = X0;
            solve_env1.last_norm = 1e10;
            solve_env1.iter_count = 0;
            solve_env1.alpha_CO = alpha_CO;
            solve_env1.alpha_MgO = alpha_MgO;
            solve_env1.dimensionless_mode = true;
            solve_env1.r_p = r_p;
            solve_env1.T_p = T_p;
            solve_env1.T_c = T_c;
            solve_env1.r_c = r_c;
            solve_env1.m_dot_total = m_dot_total_region1;
            solve_env1.Q_0 = Q_0;
            solve_env1.stage_number = 1;
            solve_env1.cp_gas = cp_gas ;
            solve_env1.k_gas = k_gas ;
            solve_env1.rho_D_gas = rho_D_gas ;
            solve_env1.porosity = porosity ;
            solve_env1.tortuosity = tortuosity ;
            solve_env1.D_ox_gas = D_ox_gas ;
            solve_env1.ratio_k_cp_ox = ratio_k_cp_ox ;
            solve_env1.k_ox = k_ox ;
            solve_env1.k_ox_gas = k_ox_gas ;
            
            % 简化的fsolve选项
            % options1 = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 10000, ...
            %            'FunctionTolerance', 1e-4, 'StepTolerance', 1e-8, 'MaxFunctionEvaluations', 10000, ...
            %            'OutputFcn', @(x,optimValues,state) monitor_progress(x,optimValues,state,solve_env1));
            options1 = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIterations', 10000, ...
                        'FunctionTolerance', 2e-3, 'StepTolerance', 1e-4, 'MaxFunctionEvaluations', 10000, ...
                        'OutputFcn', @(x,optimValues,state) monitor_progress(x,optimValues,state,solve_env1));
                    
            % 求解简化模型
            [X1, fval1, exitflag1] = fsolve(@(X) equations_system(X, solve_env1), X0, options1);
            
            % 评估阶段1解的质量
            [stage1_ok, quality1] = assess_solution_quality(X1, fval1, solve_env1);
            
            % 如果阶段1不收敛，直接返回失败
            if exitflag1 <= 0 || ~stage1_ok
                fprintf('阶段1求解未收敛，终止计算\n');
                fprintf('exitflag1: %d\n',exitflag1);
                rate_info = create_failed_result_struct();
                return;
            end
            Pe_m_region1 = X1(2);
            m_dot_total_region1 = Pe_m_region1 * 4 * pi * r_c * rho_D_gas;
            m_dot_region1_Mg = m_dot_total_region1 *X1(12);
            Y_T_surf_nd = X1(9);
            T_p_nd = Y_T_surf_nd * T_c;
            Y_flame_nd = X1(19);
            T_flame_nd = Y_flame_nd * T_c;
            fprintf('  > T_p_nd= %.6e \n', T_p_nd);
            fprintf('  > T_flame_nd= %.6e \n', T_flame_nd);
            fprintf('  > m_dot_region1_Mg= %.6e \n', m_dot_region1_Mg);
        % Y_T_inf_nd = X1(19);
        % T_f = Y_T_inf_nd * T_c;
        % fprintf('  > 阶段1求解完成，火焰温度: %.3e K\n', T_f);
        % X_final = convert_from_dimensionless(X1, r_p, r_c, T_p, T_c, m_dot_total);
        % % 提取主要参数
        % r_f = X_final(1);
        % Pe_m = X_final(2);
        % r_p_nd = r_p / r_c;
        % r_f_nd = r_f / r_c;  % 无量纲火焰半径，相对于金属核心半径
        % % 提取系数用于分析和后处理
        % [coeffs_0, coeffs_1, coeffs_2] = extract_solution_coeffs(X_final);
        % Y_T_inf_nd = X_final(19);
        % T_f = Y_T_inf_nd * T_c;
        % rate_info = package_rate_info(m_dot_region0_total, r_f, T_f, params, r_p, r_c, coeffs_0, coeffs_1, coeffs_2);
        % plot_radial_distributions(r_p, r_f, r_inf, Pe_m, coeffs_0, coeffs_1, coeffs_2, T_p, T_c, params, physicalModel,pState);
        % --- 阶段2: 求解中度简化模型 ---
            fprintf('  > 阶段2/3: 求解中度简化模型...\n');
            
            solve_env2.simplified_mode = true;
            solve_env2.simplified_reaction = true;
            solve_env2.relaxation_level = 0.5;
            solve_env2.Pe_m_factor = 0.5;
            solve_env2.Le = Le;
            solve_env2.bvp_deps = bvp_deps;
            solve_env2.last_X = X1;
            solve_env2.last_norm = 1e10;
            solve_env2.iter_count = 0;
            solve_env2.alpha_CO = alpha_CO;
            solve_env2.alpha_MgO = alpha_MgO;
            solve_env2.dimensionless_mode = true;
            solve_env2.r_p = r_p;
            solve_env2.T_p = T_p;
            solve_env2.T_c = T_c;
            solve_env2.r_c = r_c;
            solve_env2.m_dot_total = m_dot_total_region1;
            solve_env2.Q_0 = Q_0;
            solve_env2.stage_number = 2;
            solve_env2.cp_gas = cp_gas ;
            solve_env2.k_gas = k_gas ;
            solve_env2.rho_D_gas = rho_D_gas ;
            solve_env2.porosity = porosity ;
            solve_env2.tortuosity = tortuosity ;
            solve_env2.D_ox_gas = D_ox_gas ;
            solve_env2.ratio_k_cp_ox = ratio_k_cp_ox ;
            solve_env2.k_ox = k_ox ;
            solve_env2.k_ox_gas = k_ox_gas ;
            % 简化的fsolve选项
            options2 = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt','Display', 'off', 'MaxIterations', 1000000, ...
                        'FunctionTolerance', 1e-3, 'StepTolerance', 1e-4, 'MaxFunctionEvaluations', 1000000, ...
                        'OutputFcn', @(x,optimValues,state) monitor_progress(x,optimValues,state,solve_env2));
            
            % 使用阶段1的解作为初值
            [X2, fval2, exitflag2] = fsolve(@(X) equations_system(X, solve_env2), X1, options2);
            
            % 评估阶段2解的质量
            [stage2_ok, quality2] = assess_solution_quality(X2, fval2, solve_env2);
            
            % 简化逻辑：如果阶段2不收敛，直接返回失败
            if exitflag2 <= 0 || ~stage2_ok
                fprintf('阶段2求解未收敛，终止计算\n');
                fprintf('exitflag2: %d\n',exitflag2);
                rate_info = create_failed_result_struct();
                return;
            end
            Pe_m_region1 = X2(2);   
            m_dot_total_region1 = Pe_m_region1 * 4 * pi * r_c * rho_D_gas;
            m_dot_region1_Mg = m_dot_total_region1 *X2(12);
            Y_T_surf_nd = X2(9);
            T_p_nd = Y_T_surf_nd * T_c;
            Y_flame_nd = X2(19);
            T_flame_nd = Y_flame_nd * T_c;
            fprintf('  > T_p_nd= %.6e \n', T_p_nd);
            fprintf('  > T_flame_nd= %.6e \n', T_flame_nd);
            fprintf('  > m_dot_region1_Mg= %.6e \n', m_dot_region1_Mg);
        % Y_T_inf_nd = X2(19);
        % T_f = Y_T_inf_nd * T_c;
        % fprintf('  > 阶段2求解完成，火焰温度: %.3e K\n', T_f);
        % 
        % X_final = convert_from_dimensionless(X2, r_p, r_c, T_p, T_c, m_dot_total);
        % % 提取主要参数
        % r_f = X_final(1);
        % Pe_m = X_final(2);
        % r_p_nd = r_p / r_c;
        % r_f_nd = r_f / r_c;  % 无量纲火焰半径，相对于金属核心半径
        % % 提取系数用于分析和后处理
        % [coeffs_0, coeffs_1, coeffs_2] = extract_solution_coeffs(X_final);
        % Y_T_inf_nd = X_final(19);
        % T_f = Y_T_inf_nd * T_c;
        % rate_info = package_rate_info(m_dot_region0_total, r_f, T_f, params, r_p, r_c, coeffs_0, coeffs_1, coeffs_2);
        % plot_radial_distributions(r_p, r_f, r_inf, Pe_m, coeffs_0, coeffs_1, coeffs_2, T_p, T_c, params, physicalModel,pState);

        % --- 阶段3: 求解完整模型 ---
        fprintf('  > 阶段3/3: 求解完整模型...\n');
        
            solve_env3.simplified_mode = false; 
            solve_env3.simplified_reaction = false;
            solve_env3.relaxation_level = 0.9;
            solve_env3.Pe_m_factor = 1.0;
            solve_env3.Le = Le;
            solve_env3.bvp_deps = bvp_deps;
            solve_env3.last_X = X2;
            solve_env3.last_norm = 1e10;
            solve_env3.iter_count = 0;
            solve_env3.alpha_CO = alpha_CO;
            solve_env3.alpha_MgO = alpha_MgO;
            solve_env3.dimensionless_mode = true;
            solve_env3.r_p = r_p;
            solve_env3.T_p = T_p;
            solve_env3.T_c = T_c;
            solve_env3.r_c = r_c;
            solve_env3.m_dot_total = m_dot_total_region1;
            solve_env3.Q_0 = Q_0;
            solve_env3.stage_number = 3;
            solve_env3.cp_gas = cp_gas ;
            solve_env3.k_gas = k_gas ;
            solve_env3.rho_D_gas = rho_D_gas ;
            solve_env3.porosity = porosity ;
            solve_env3.tortuosity = tortuosity ;
            solve_env3.D_ox_gas = D_ox_gas ;
            solve_env3.ratio_k_cp_ox = ratio_k_cp_ox ;
            solve_env3.k_ox = k_ox ;
            solve_env3.k_ox_gas = k_ox_gas ;
            % 简化的fsolve选项
            options3 = optimoptions('fsolve','Algorithm', 'levenberg-marquardt', 'Display', 'off', 'MaxIterations', 100000, ...
                        'FunctionTolerance', 1e-4, 'StepTolerance', 1e-6, 'MaxFunctionEvaluations', 100000, ...
                        'Algorithm', 'levenberg-marquardt', ...
                        'OutputFcn', @(x,optimValues,state) monitor_progress(x,optimValues,state,solve_env3));
            
            % 使用阶段2的解作为初值
            [X_final_nd, fval3, exitflag3] = fsolve(@(X) equations_system(X, solve_env3), X2, options3);
            
            % 评估最终解的质量
            [stage3_ok, quality3] = assess_solution_quality(X_final_nd, fval3, solve_env3);

            %
            %X_final = convert_from_dimensionless(X_final_nd, r_p, r_c, T_p, T_c, m_dot_total);
            %r_f = X_final(1);
            %Pe_m = X_final(2);
            %[coeffs_0, coeffs_1, coeffs_2] = extract_solution_coeffs(X_final);
            %plot_radial_distributions(r_p, r_f, r_inf, Pe_m, coeffs_0, coeffs_1, coeffs_2, T_p, T_c, params, physicalModel,pState);

            % 如果阶段3不收敛，使用前一阶段的结果
            if exitflag3 <= 0 || ~stage3_ok
                fprintf('阶段3求解未收敛，终止计算\n');
                fprintf('exitflag3: %d\n',exitflag3);
                rate_info = create_failed_result_struct();
                return;
            end
            
            fprintf('    * 最终解质量评分: %.1f/100\n', quality3);
                    
        Y_T_inf_nd = X_final_nd(19);
        T_f = Y_T_inf_nd * T_c;
        fprintf('  > 阶段3求解完成，火焰温度: %.3e K\n', T_f);
        % 从无量纲解转换回物理量
        X_final = convert_from_dimensionless(X_final_nd, r_p, r_c, T_p, T_c, m_dot_total_region1);
        % 提取主要参数
        r_f = X_final(1) * r_c;
        Pe_m_region1 = X_final(2);
        Y_T_core_nd = X_final(3);
        Y_T_surf_nd = X_final(9);
        m_dot_total_region1 = Pe_m_region1 * 4 * pi * r_c * rho_D_gas;
        r_p_nd = r_p / r_c;
        r_f_nd = r_f / r_c;  % 无量纲火焰半径，相对于金属核心半径
        % 提取系数用于分析和后处理
        [coeffs_0, coeffs_1, coeffs_2] = extract_solution_coeffs(X_final);
        m_dot_region1_Mg= m_dot_total_region1 * coeffs_1.Y_Mg_frac;
        fprintf('  > mg质量m_dot_region1_Mg: %.3e K\n', m_dot_region1_Mg);
        fprintf('  > Mg核表面温度: %.3e K\n', Y_T_core_nd * T_c);
        fprintf('  > 颗粒表面温度: %.3e K\n', Y_T_surf_nd * T_c);
        Y_T_flame_nd = X_final(19);
        T_f = Y_T_flame_nd * T_c;
            % 检查结果是否合理
            success = true;
            if Pe_m_region1 < 0 || r_f <= r_p || isnan(T_f) || T_f <= 0
                fprintf('警告: 求解结果不合理，Pe_m=%.3e, r_f/r_p=%.2f, T_f=%.2f\n', ...
                    Pe_m_region1, r_f/r_p, T_f);
                success = false;
            end
            
        % 封装结果，传递所有系数用于组分流率计算
        rate_info = package_rate_info(m_dot_total_region1, r_f, T_f, params, r_p, r_c, coeffs_0, coeffs_1, coeffs_2,T_c);
        rate_info.success = success;
        
        % 存储最终解供下一时间步使用
        rate_info.algebraic_solution = X_final_nd;
        
        fprintf('区域特定解析解求解完成。火焰位置: %.3e m (r_f/r_p = %.2f), 蒸发速率: %.3e kg/s\n', ...
            r_f, r_f/r_p, rate_info.dmdt_mg);
        fprintf('火焰温度: %.3e K\n', T_f);
        if m_mg / rate_info.dmdt_mg > 0.005
            fprintf('警告: 蒸发速率与质量不匹配，m_mg/dmdt_mg = %.3e\n', m_mg / rate_info.dmdt_mg);
        end
        % 在成功求解后添加可选的可视化
        if exitflag3 > 0 && success
            try
                visualization = params.plot_radial_distributions ; 
                if visualization
                    plot_radial_distributions(r_p, r_f, r_inf, Pe_m_region1, coeffs_0, coeffs_1, coeffs_2, T_p, T_c, params, physicalModel,pState,cp_gas...
                        ,k_gas,rho_D_gas,porosity,tortuosity,D_ox_gas,ratio_k_cp_ox,ratio_k_cp_gas,k_ox,k_ox_gas);
                    fprintf('已生成求解结果的径向分布图\n');
                end
            catch ME
                fprintf('结果可视化失败，但不影响计算结果: %s\n', ME.message);
            end
        end
        
    catch ME
        fprintf('求解过程中发生错误: %s\n', ME.message);
        fprintf('错误位置: %s, 行: %d\n', ME.stack(1).file, ME.stack(1).line);
        fprintf('--- 气相扩散燃烧代数方程求解器失败 ---\n');
        rate_info = create_failed_result_struct();
        return;  % 确保在catch块中也直接返回
    end
    

end

% =========================================================================
% 辅助函数
% =========================================================================

function X0 = generate_initial_guess(m_dot_total,r_f, Pe_m, r_p, r_c,r_inf, T_p, T_c, T_amb, Le, mw, physicalModel, params,cp_gas...
    ,k_gas,rho_D_gas,porosity,tortuosity,D_ox_gas,ratio_k_cp_ox,ratio_k_cp_gas,k_ox,k_ox_gas);
    % 为所有未知数生成无量纲化的初始猜测向量
    % 注意：无量纲化基准已改为金属核心半径r_c和金属核心温度T_c
    T_f_guess = 2200 ;
    k = 0.2;
    % 气体物性的无量纲参数     co2  60 j/mol/k
    % 气体比热估计值 [J/(kg·K)]
    % co2  2000K
     % 气体导热系数估计值 [W/(m·K)]
    %fprintf('  > le_solid: = %.6e\n', le_solid);
    % 假定T_c = T_p (初始猜测阶段)

    r_p_nd = r_p / r_c;        % 无量纲颗粒外表面半径
    r_f_nd = r_f / r_c;        % 无量纲火焰半径
    r_inf_nd = r_inf / r_c;    % 无量纲远场边界
    T_amb_nd = T_amb / T_c;    % 无量纲环境温度（基于核心温度）
    T_p_nd = T_p / T_c;        % 无量纲氧化层表面温度（基于核心温度）



    r_ini = params.initial_diameter/2;
    oxide_ratio = (r_p - r_c) / r_ini;  % 氧化层厚  度与初始半径比


    %%%%%基于区域1进行展开 
        % 核心表面Mg质量分数
        pressure = params.ambient_pressure / 101325 ;
        A = 4.32e5;
        B = 1.77e4;
        P_Mg_surface = A * exp(-B / T_c);
        V_frac_Mg = P_Mg_surface / pressure;
        V_frac_CO = 1 - V_frac_Mg;
        Y_Mg_theory = V_frac_Mg * mw.Mg.molar_mass / (V_frac_Mg * mw.Mg.molar_mass + V_frac_CO * mw.CO.molar_mass);
        %fprintf('T_c  %.6e \n',T_c);
        %fprintf('P_Mg_surface  %.6e \n',P_Mg_surface);
        %fprintf('Y_Mg_theory  %.6e \n',Y_Mg_theory);
        %fprintf('V_frac_Mg  %.6e \n',V_frac_Mg);
        % Y_Mg_theory = physicalModel.cal_Mg_frac_surface(pressure,T_c);
        Y_Mg_core = Y_Mg_theory;  
        Y_CO_core = 1 - Y_Mg_core; 

    %表面反应
        %fprintf('Y_Mg_region0_surf  %.6e \n',Y_Mg_region0_surf);
        K_CO_surface = params.reaction_pre_exponential * exp(-params.reaction_activation_energy / (params.R_u * T_c));
        %M_mix_surface = Y_Mg_region0_surf * mw.Mg.molar_mass + Y_CO_region0_surf * mw.CO.molar_mass;
        M_mix_surface = Y_Mg_core * mw.Mg.molar_mass + Y_CO_core * mw.CO.molar_mass;

        A_surface = 4 * pi * r_c^2 ;
        %m_dot_surface_Mg = A_surface * Y_CO_region0_surf * pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass  *  K_CO_surface;
        m_dot_surface_Mg = A_surface * Y_CO_core * pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass  *  K_CO_surface;
        m_dot_surface_CO = m_dot_surface_Mg * mw.CO.molar_mass / mw.Mg.molar_mass;
        m_dot_surface_MgO = m_dot_surface_Mg * mw.MgO.molar_mass / mw.Mg.molar_mass;

    % 组分质量流率相对于总质量流率
        alfa_mg = 1.5;    % Mg比例
        alfa_co=  -0.8;   % CO比例
        alfa_co2 = 0.0;   % CO2比例
        alfa_mgo = -0.3;  % MgO比例
        m_dot_region1_total = m_dot_total;
        m_dot_region1_CO = m_dot_region1_total * alfa_co;
        m_dot_region1_Mg = m_dot_region1_total * alfa_mg;
        m_dot_region1_CO2 = m_dot_region1_total * alfa_co2;
        m_dot_region1_MgO = m_dot_region1_total * alfa_mgo;

        m_dot_region0_Mg = m_dot_region1_Mg;
        m_dot_region0_CO = -(abs(m_dot_region1_CO) - m_dot_surface_CO);
        m_dot_region0_total = m_dot_region0_Mg + m_dot_region0_CO;
        Y_CO_frac0 = m_dot_region0_CO / m_dot_region0_total;
        Y_Mg_frac0 = m_dot_region0_Mg / m_dot_region0_total ;

        m_dot_region2_Mg = 0;
        m_dot_region2_CO2 = -m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO2.molar_mass;
        m_dot_region2_CO = m_dot_region1_Mg / mw.Mg.molar_mass*mw.CO.molar_mass - abs(m_dot_region1_CO);
        m_dot_region2_MgO = m_dot_region1_Mg / mw.Mg.molar_mass*mw.MgO.molar_mass - abs(m_dot_region1_MgO);
        m_dot_region2_total = m_dot_region2_Mg + m_dot_region2_CO2+m_dot_region2_CO+m_dot_region2_MgO;


        Pe_m_region0 = m_dot_region0_total / (4 * pi * r_c * D_ox_gas );
        Pe_cp_region0 = m_dot_region0_total * cp_gas / (4 * pi * r_c * k_ox_gas);
        %Pe_cp_region0 = m_dot_region0_total  / (4 * pi * r_c * ratio_k_cp_ox);
        Pe_m_region1 = Pe_m ;
        Pe_cp_region1 = m_dot_region1_total * cp_gas / (4 * pi * r_c * k_gas);
        %Pe_cp_region1 = m_dot_region1_total / (4 * pi * r_c * ratio_k_cp_gas);

        Pe_m_region2 = m_dot_region2_total / (4 * pi * r_c * rho_D_gas );
        Pe_cp_region2 = m_dot_region2_total * cp_gas / (4 * pi * r_c * k_gas);
        %Pe_cp_region2 = m_dot_region2_total / (4 * pi * r_c * ratio_k_cp_gas);
    %Pe_cp_region0 = m_dot_region0_total  / (4 * pi * r_c * ratio_k_cp_ox);


    exp_arg_region0_Y = Pe_m_region0 * (1 - 1/r_p_nd);
    Y_Mg_region0_surf = Y_Mg_core + (Y_Mg_core - Y_Mg_frac0) * (exp(exp_arg_region0_Y) - 1);
    Y_CO_region0_surf = Y_CO_core + (Y_CO_core - Y_CO_frac0) * (exp(exp_arg_region0_Y) - 1);

% B0 计算
    L_v = mw.Mg.L_evap_Mg / mw.Mg.molar_mass;
    Q_evap = m_dot_region0_Mg * L_v;
    H_reac_f = params.reaction_heats.flame_reac_H;  % [J/kg]
    Q_flame = m_dot_region0_Mg * abs(H_reac_f);
    T_oxshell = T_p;  % 初始猜测中，氧化层表面温度等于T_p   
    %%%%%%%%%%此处猜测 较大温差时间 计算问题

    %Q_rad_oxshell = 4 * pi * r_c^2 * params.emissivity * params.sigma * porosity * (T_oxshell^4 - T_c^4);
    Q_rad_oxshell = 4 * pi * r_c^2 * params.emissivity * params.sigma * (T_oxshell^4 - T_c^4) ;
    
    %   表面反应热 全用来加热 
    H_reac_CO = params.reaction_heats.surface_reac_H;
    Q_reaction_face  = abs(m_dot_surface_Mg) * (-H_reac_CO) ;

    % mg不直接吸收 气相反应热 
    H_reac_Mg_flame= params.reaction_heats.flame_reac_H;
    %Q_reaction_flame = abs(m_dot_region0_mg) * (-H_reac_Mg_flame) ;
    Q_reaction_flame = 0 ;
    Q_reaction = Q_reaction_face + Q_reaction_flame * k;

    %对流部分   不作用于mg核
    q_conv = params.k_gas / r_p * (T_amb -T_p);
    heat_convection = q_conv * (4 * pi * r_p^2);
    heat_convection = 0 ;


    %  mg核表面能量守恒 蒸发=表面反应 + 辐射 + 镁右侧导热（此处为 kox * dT）
    %  A表示导热部分   来自于 右侧
    A =  Q_evap - Q_reaction - Q_rad_oxshell - heat_convection ; 
    %%%%%  区域0能量方程积分  在此处为 koxgas * dt  （B0推导针对能量方程展开 、此处对 A变化）
    A_porosity = A  * porosity ;   % 区域一积分后的导热  用于后续计算； 
    Q_sens_region0 = m_dot_region0_total * cp_gas * T_c ;

    % B0 
    B0 =  A_porosity - Q_sens_region0 ;
    % 无量纲化 两步   先显热部分 再进行温度处理  此处pemcp 中为koxgas
    B0_nd = B0 / (m_dot_region0_total * cp_gas * T_c );
    %fprintf(' 镁核表面 m_dot_region0_total: %.3e \n', m_dot_region0_total);
    fprintf(' 镁核表面温度梯度: %.3e \n', A_porosity);
    %fprintf(' 镁核表面 B0_nd: %.3e \n', B0_nd);
    %fprintf(' MG表面 Q_cod_region0: %.3e \n', Q_cod_region0);
    % 区域0初始猜测（氧化层内）
    Y_T_core_nd = 1.0;                % 核心表面温度为沸点温度(无量纲，基于T_c)               % 初始猜测B0参数
    

    % 1. 全局参数
    X0(1) = r_f_nd;            
    X0(2) = Pe_m_region1;                
    
    % 2. 区域0: 核心表面到氧化层表面
    X0(3) = Y_T_core_nd;             % 核心表面温度(无量纲)
    X0(4) = B0_nd;                   % 区域0的无量纲B参数
    X0(5) = Y_Mg_core;               % Mg核心表面质量分数
    X0(6) = Y_Mg_frac0;              % Mg在区域0质量流率中的比例
    X0(7) = Y_CO_core;               % CO核心表面质量分数
    X0(8) = Y_CO_frac0;              % CO在区域0质量流率中的比例
    

    exp_term_p = exp(Pe_cp_region0 * (1-1/r_p_nd)) - 1;
    T_p_nd = Y_T_core_nd + (Y_T_core_nd + B0_nd) * exp_term_p;
    T_p_region0 = T_p_nd * T_c;

    exp_arg_region0_Y = Pe_m_region0 * (1 - 1/r_p_nd);
    Y_Mg_region0_surf = Y_Mg_core + (Y_Mg_core - Y_Mg_frac0) * (exp(exp_arg_region0_Y) - 1);
    Y_CO_region0_surf = Y_CO_core + (Y_CO_core - Y_CO_frac0) * (exp(exp_arg_region0_Y) - 1);
    %fprintf('Y_Mg_core  %.6e \n',Y_Mg_core);
    %fprintf('Y_Mg_frac0  %.6e \n',Y_Mg_frac0);
    %fprintf('(exp(exp_arg_region0_Y) - 1)  %.6e \n',(exp(exp_arg_region0_Y) - 1));
    %fprintf('Y_Mg_region0_surf  %.6e \n',Y_Mg_region0_surf);
    %K_CO_surface = params.reaction_pre_exponential * exp(-params.reaction_activation_energy / (params.R_u * T_c));
    %M_mix_surface = Y_Mg_region0_surf * mw.Mg.molar_mass + Y_CO_region0_surf * mw.CO.molar_mass;
    %A_surface = 4 * pi * r_p^2 * porosity;
    %m_dot_surface_Mg = A_surface * Y_CO_region0_surf * pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass ...
    %                    *  K_CO_surface;
    %m_dot_surface_CO = m_dot_surface_Mg * mw.CO.molar_mass / mw.Mg.molar_mass;
    %m_dot_surface_MgO = m_dot_surface_Mg * mw.MgO.molar_mass / mw.Mg.molar_mass;
    % 3. 区域1: 氧化层表面到火焰面
 
            %    fprintf(' 初始猜测 Pe_cp_region1: %.3e \n', Pe_cp_region1);
            % 计算无量纲能量参数

            % MgO沉积热贡献  Q_dep_MgO 绝对值  J/kg
            H_dep_MgO = params.H_dep_MgO;
            Q_dep_MgO = abs(m_dot_region1_MgO)  * (-H_dep_MgO);
            %Q_dep_MgO = 0 ;

            % 火焰温度猜测
            T_amb = params.ambient_temperature ; 

            H_reac_CO = params.reaction_heats.surface_reac_H;
            Q_reaction_face  = abs(m_dot_surface_CO) * (-H_reac_CO) ;
            Q_reaction_face = 0 ;

            % 辐射热贡献
            emissivity = params.emissivity;
            sigma = params.sigma ;
            Q_rad_f = emissivity * sigma * (4 * pi * r_p^2) * (T_f_guess^4 - T_p^4);
            %Q_rad_f = 0;
            Q_rad_amb = params.emissivity * params.sigma * (4 * pi * r_p^2) * (T_amb^4 - T_p^4);
            Q_rad_amb = 0 ;
            Q_rad = Q_rad_f+ Q_rad_amb ;
            
            %对流部分 
            q_conv = params.k_gas / r_p * (T_amb -T_p);
            heat_convection = q_conv * (4 * pi * r_p^2);
            %heat_convection = 0 ;

            %  气相反应直接被颗粒吸收部分
            H_reac_f = params.reaction_heats.flame_reac_H;
            Q_reac_f  =  abs(m_dot_region1_Mg) *(-H_reac_f) * k ;
            %Q_reac_f = 0 ;

            % 类似A 此处用B表示能量
            %%%%%%%%%%%来源于 颗粒表面能量      （火焰+环境）辐射 +沉积 = 该表面处向两侧的热传导 + 向颗粒的辐射
            %%%%%%此处热传导中 颗粒左侧表面部分  koxgas * DT   右侧为 对应区域的kgas

            %左侧热传导计算 即 koxgas * DT    根据0区域积分后能量方程  此处为A_porosity 不是 A 
            % 
            Q_conv_left = (A_porosity + m_dot_region0_total * cp_gas * ( T_p - T_c ) ) ;
            % B表示颗粒右侧导热部分   即kgas* dT
            B = (Q_rad  + Q_reaction_face + Q_dep_MgO + heat_convection + Q_reac_f- ...
                    Q_conv_left - Q_rad_oxshell) ;
            fprintf(' 颗粒表面 Q_conv_left : %.3e \n', Q_conv_left);
            fprintf(' 颗粒表面 B : %.3e \n', B );
            %Q_cod_region0 = Q_evap - Q_reaction - Q_rad_oxshell ;   Q_sens_region0 = m_dot_region0_total * cp_gas * T_c / porosity;
            %   Q_cod_region0 +m_dot_region0_total*cp_gas*(T_p-T_c)/porosity   左侧梯度   求导热时 乘孔隙率
            % 显热贡献
            Q_sens_region1 = m_dot_region1_total * cp_gas * T_p;
            %fprintf(' 颗粒表面 Q_conv_left: %.3e: Q_conv_right: %.3e\n', Q_conv_left,Q_conv_right);
            %Q_sens_region1 = m_dot_region1_total * cp_gas * T_f_guess;    
            %fprintf(' 颗粒表面 Q_conv_right: %.3e: Q_sens_region1:\n', Q_conv_right,Q_sens_region1);
            % 计算X1   物理意义：
            B1 = B - Q_sens_region1;
            
            % 转换为无量纲B1
            B1_nd = B1 / (m_dot_region1_total * cp_gas * T_c);
            %fprintf(' 颗粒表面 m_dot_region1_total : %.3e \n', m_dot_region1_total);
            %fprintf(' 颗粒表面 Q_sens_region1 : %.3e \n', Q_sens_region1 );
            fprintf(' 颗粒表面 B1_nd : %.3e \n', B1_nd );
            % 表面温度的无量纲值 (归一化为1.0)
            Y_T_surf_nd = T_p_nd;
            m_dot_B1_cp_ratio_nd = B1_nd;
        
        X0(9) = Y_T_surf_nd;            % 表面无量纲温度(为1.0)
        X0(10) = m_dot_B1_cp_ratio_nd;   % 无量纲B1
        
        
        %fprintf('Y_Mg_core  %.6e \n',Y_Mg_core);
        %fprintf('Y_Mg_frac0  %.6e \n',Y_Mg_frac0);
        %fprintf('exp(exp_arg_region0_Y) - 1  %.6e \n',exp(exp_arg_region0_Y) - 1);
        %fprintf('(Y_Mg_core - Y_Mg_frac0) * (exp(exp_arg_region0_Y) - 1)  %.6e \n',(Y_Mg_core - Y_Mg_frac0) * (exp(exp_arg_region0_Y) - 1));
        %fprintf('Y_Mg_region0_surf  %.6e \n',Y_Mg_region0_surf);
        %fprintf('Y_CO_region0_surf  %.6e \n',Y_CO_region0_surf);
        alfa_Mg_nd = m_dot_region1_Mg / m_dot_region1_total;
        alfa_CO_nd = m_dot_region1_CO / m_dot_region1_total;
        alfa_CO2_nd = m_dot_region1_CO2 / m_dot_region1_total;
        alfa_MgO_nd = m_dot_region1_MgO / m_dot_region1_total;

            Y_Mg_surf = Y_Mg_region0_surf;
            Y_Mg_frac1 = alfa_Mg_nd;
            % CO2组分参数 (无量纲)
            Y_CO2_surf = 0.0;
            Y_CO2_frac1 = alfa_CO2_nd;
            
            % CO组分参数 (无量纲)
            Y_CO_surf = Y_CO_region0_surf;
            Y_CO_frac1 = alfa_CO_nd;
                
            % MgO组分参数 (无量纲)
            Y_MgO_surf = 0.0;
            Y_MgO_frac1 = alfa_MgO_nd;
            
            %Y_total = Y_CO2_surf + Y_CO_surf + Y_MgO_surf+Y_Mg_surf;
            %Y_Mg_surf = Y_Mg_surf / Y_total;
            %Y_CO2_surf = Y_CO2_surf / Y_total;
            %Y_CO_surf = Y_CO_surf / Y_total;
            %Y_MgO_surf = Y_MgO_surf / Y_total;

        X0(11) = Y_Mg_surf;
        X0(12) = Y_Mg_frac1;
        
    
        X0(13) = Y_CO2_surf;
        X0(14) = Y_CO2_frac1;
    
        X0(15) = Y_CO_surf;
        X0(16) = Y_CO_frac1;

        X0(17) = Y_MgO_surf;
        X0(18) = Y_MgO_frac1;
    
    % 4. 区域2: 火焰面到远场
        % 区域2的质量流率计算
            % 区域2中的组分分布 
            Y_Mg_frac2 = m_dot_region2_Mg /m_dot_region2_total;       % 区域2无Mg
            Y_CO2_frac2 = m_dot_region2_CO2 / m_dot_region2_total;     % 区域2 CO2质量流率比例
            Y_CO_frac2 = m_dot_region2_CO / m_dot_region2_total;      % 区域2 CO质量流率比例
            Y_MgO_frac2 = m_dot_region2_MgO / m_dot_region2_total;     % 区域2 MgO质量流率比例
 
            % 计算火焰面的无量纲温度 (区域1解析解)
            exp_term_f = exp(Pe_cp_region1 * (1/r_p_nd-1/r_f_nd)) - 1;
            T_f_nd = Y_T_surf_nd + (Y_T_surf_nd + m_dot_B1_cp_ratio_nd) * exp_term_f;

            debug = params.debug;
            if debug 
            fprintf('初始猜测结果：\n')
            fprintf(' m_dot_region2_CO2: %.3e \n', m_dot_region2_CO2);
            fprintf(' m_dot_region2_CO: %.3e \n', m_dot_region2_CO);
            fprintf(' m_dot_region2_MgO: %.3e \n', m_dot_region2_MgO);
            fprintf(' m_dot_region2_Mg: %.3e \n', m_dot_region2_Mg);
            fprintf(' m_dot_region2_total: %.3e \n', m_dot_region2_total);
            fprintf(' 初始猜测根据颗粒表面温度: %.3e K\n', T_p);
            fprintf(' 初始猜测根据颗粒表面能量 计算得到火焰面温度: %.3e K\n', T_f_nd*T_c);
            end 
                T_f = T_f_nd * T_p;  % 物理火焰温度
                r_f = r_f_nd * r_p;  % 物理火焰半径
                
            % 计算区域2的物理量B2参数
                
                % 辐射热贡献
                Q_rad_f_out = emissivity * sigma * (4 * pi * r_f^2) * (T_f^4 );
                Q_rad_f_out = emissivity * sigma * (4 * pi * r_f^2) * (T_f^4 - T_amb^4);
                Q_rad_f_in = emissivity * sigma *  (4 * pi * r_f^2) * (T_f^4 - T_p^4);
                Q_rad_total = Q_rad_f_out + Q_rad_f_in;

                % 反应热贡献   
                H_reac_f = params.reaction_heats.flame_reac_H;
                Q_reac_f  =  abs(m_dot_region1_Mg) *(-H_reac_f)* (1-k);
                %Q_reac_f  =  abs(m_dot_region1_Mg) *(-H_reac_f) * 0.5;
                % 颗粒表面到火焰面 积分后 表示   注意与B1的区别   B1 = B - Q_sens_region1
                Q_cond_left = m_dot_region1_total * cp_gas * (T_f - T_p) + B ;
                fprintf(' 火焰面 Q_cond_left: %.3e \n', Q_cond_left);
                fprintf(' m_dot_region1_Mg: %.3e \n', m_dot_region1_Mg);
                %%%%% C 表示火焰面右侧导热 部分
                C = -(Q_reac_f- Q_rad_total - Q_cond_left) ;
    
                Q_sens_region2 = m_dot_region2_total * cp_gas * T_f;
              
                B2 = C - Q_sens_region2;
                B2_nd = B2 / (m_dot_region2_total * cp_gas* T_c );

                %B2_nd = B2/ T_p;
                fprintf(' 火焰面 C: %.3e \n', C);
        X0(19) = T_f_nd;            % 火焰面无量纲温度
        X0(20) = B2_nd;             % 无量纲B2
        
        % 计算火焰面各组分质量分数 (使用区域1解析解)
        exp_term_Y = exp( Pe_m_region1 * ( 1/r_p_nd - 1/r_f_nd)) - 1;
        Y_Mg_f = Y_Mg_surf + (Y_Mg_surf - Y_Mg_frac1) * exp_term_Y;
        Y_CO2_f = Y_CO2_surf + (Y_CO2_surf - Y_CO2_frac1) * exp_term_Y;
        Y_CO_f = Y_CO_surf + (Y_CO_surf - Y_CO_frac1) * exp_term_Y;
        Y_MgO_f = Y_MgO_surf + (Y_MgO_surf - Y_MgO_frac1) * exp_term_Y;
        Y_Mg_f = 0 ;
        Y_CO2_f = 0 ;
        %fprintf('  Y_CO_f = %.3e, Y_MgO_f = %.3e\n', Y_CO_f, Y_MgO_f);
        Y_CO_f  =  Y_CO_f / (Y_CO_f + Y_MgO_f);
        Y_MgO_f = 1 - Y_CO_f ;
        %fprintf('  Y_CO_f = %.3e, Y_MgO_f = %.3e\n', Y_CO_f, Y_MgO_f);
        exp_term_inf = exp(Pe_cp_region2 * (1/r_f_nd-1/r_inf_nd)) - 1;
        T_inf_nd = T_f_nd + (T_f_nd + B2_nd) * exp_term_inf;

    % 区域2的组分参数
        X0(21) = Y_Mg_f;       % 火焰面Mg质量分数
        X0(22) = Y_Mg_frac2;   % 区域2 Mg质量流率比例
        X0(23) = Y_CO2_f;      % 火焰面CO2质量分数
        X0(24) = Y_CO2_frac2;  % 区域2 CO2质量流率比例
        X0(25) = Y_CO_f;       % 火焰面CO质量分数
        X0(26) = Y_CO_frac2;   % 区域2 CO质量流率比例
        X0(27) = Y_MgO_f;      % 火焰面MgO质量分数
        X0(28) = Y_MgO_frac2;  % 区域2 MgO质量流率比例
    
    % 初始猜测结果的可视化
    if debug 
        fprintf('区域0参数(无量纲):\n');
        fprintf('  Y_T_core = %.3e, B0 = %.3e\n', X0(3), X0(4));
        fprintf('  Y_Mg_core = %.3e, Y_Mg_frac0 = %.3e\n', X0(5), X0(6));
        fprintf('  Y_CO_core = %.3e, Y_CO_frac0 = %.3e\n', X0(7), X0(8));
        
        fprintf('区域1参数(无量纲):\n');
        fprintf('  Y_T_surf_nd = %.3e, B1 = %.3e\n', X0(9)*T_c, X0(10));
        fprintf('  Y_Mg_surf = %.3e, Y_Mg_frac = %.3e\n', X0(11), X0(12));
        fprintf('  Y_CO2_surf = %.3e, Y_CO2_frac = %.3e\n', X0(13), X0(14));
        fprintf('  Y_CO_surf = %.3e, Y_CO_frac = %.3e\n', X0(15), X0(16));
        fprintf('  Y_MgO_surf = %.3e, Y_MgO_frac = %.3e\n', X0(17), X0(18));
        
        fprintf('区域2参数(无量纲):\n');
        fprintf('  T_f_nd = %.3e, B2_nd = %.3e\n', X0(19)*T_c, X0(20));
        fprintf('  Y_Mg_flame = %.3e, Y_Mg_frac = %.3e\n', X0(21), X0(22));
        fprintf('  Y_CO2_flame = %.3e, Y_CO2_frac = %.3e\n', X0(23), X0(24));
        fprintf('  Y_CO_flame = %.3e, Y_CO_frac = %.3e\n', X0(25), X0(26));
        fprintf('  Y_MgO_flame = %.3e, Y_MgO_frac = %.3e\n', X0(27), X0(28));
        fprintf(' T_F: %.2f K\n', T_f_nd*T_c);
        fprintf(' T_inf: %.2f K\n', T_inf_nd*T_c);
    end
    % 可视化初始猜测
       %visualize_initial_guess(X0, r_p, r_f, r_inf, T_p, T_amb, params,  physicalModel);
end

function F = equations_system(X, env)
    % 构建完整的代数方程组
    % X: 未知数向量
    % env: 包含所有依赖和求解模式的环境结构体
    % 
    % 新解析解形式 (无量纲，基准为金属核心半径r_c和金属核心温度T_c):
    % 区域0: Y(r) = (Y_core - Y_frac) * (exp(Pe_m * ((1) - (1/r_nd))) - 1) + Y_core
    % 区域1: Y(r) = (Y_surf - Y_frac) * (exp(Pe_m * ((1/r_p_nd) - (1/r_nd))) - 1) + Y_surf
    % 区域2: Y(r) = (Y_flame - Y_frac) * (exp(Pe_m * ((1/r_f_nd) - (1/r_nd))) - 1) + Y_flame
   
    k = 0.15;
    % 提取求解环境
    bvp_deps = env.bvp_deps;
    simplified_mode = env.simplified_mode;
    simplified_reaction = env.simplified_reaction;
    relaxation = env.relaxation_level;
    Pe_m_factor = env.Pe_m_factor;
    Le = env.Le;
    m_dot_total = env.m_dot_total;

    cp_gas = env.cp_gas;
    k_gas = env.k_gas;
    rho_D_gas = env.rho_D_gas;
    porosity = env.porosity;
    tortuosity = env.tortuosity;
    D_ox_gas = env.D_ox_gas;
    ratio_k_cp_ox = env.ratio_k_cp_ox;
    ratio_k_cp_gas = rho_D_gas;
    k_ox = env.k_ox;
    k_ox_gas = env.k_ox_gas;
    % 提取参数
    params = bvp_deps.params;
    pState = bvp_deps.pState;
    physicalModel = bvp_deps.physicalModel;
  
    T_amb = params.ambient_temperature;
    mw = params.materials ;
    
    %fprintf('  > le_solid: = %.6e\n', le_solid);
    % 基本物理参数
    r_p = pState.r_p;
    r_c = pState.r_c;             % 核心半径
    T_p = pState.T_p;             % 氧化层表面温度
    %T_p = 1600 ;
    T_c = pState.T_c;             % 金属核心温度 (新增)
    T_amb_nd = T_amb / T_c;       % 无量纲环境温度（相对于核心温度）
    T_p_nd = T_p / T_c;           % 无量纲氧化层表面温度（相对于核心温度）
    r_p_nd = r_p / r_c;           % 无量纲颗粒表面半径（基于核心半径）
    r_c_nd = 1 ;
    
    %cp_gas = physicalModel.get_gas_properties(T_p).Cp_gas;
    %k_gas = physicalModel.get_gas_properties(T_p).k_gas;
    
    %fprintf('  cp_gas: %.5e, k_gas: %.6e\n', cp_gas, k_gas);
        dimensionless_mode = false;
        if isfield(env, 'dimensionless_mode')
            dimensionless_mode = env.dimensionless_mode;
        end

        if dimensionless_mode
            r_p = env.r_p;
            T_p = env.T_p;
            m_dot_total = env.m_dot_total;
        end

    % 提取未知变量
    r_f_nd = X(1);                % 无量纲火焰半径 (r_f/r_c)  
    %%% 火焰半径处理 无量纲来源于上次求解值  相较于xc变化不大即无影响？
    r_f = r_f_nd * r_c;           % 物理火焰半径
    Pe_m_region1 = X(2);                  % 佩克莱数
    
    % 区域0（氧化层内）参数
    Y_T_core_nd = X(3);           % 核心表面温度(无量纲)
    B0_nd = X(4);                 % 区域0的无量纲B参数
    Y_Mg_core = X(5);             % Mg核心表面质量分数
    Y_Mg_frac0 = X(6);            % Mg在区域0总质量流率中的比例
    Y_CO_core = X(7);             % CO核心表面质量分数
    Y_CO_frac0 = X(8);            % CO在区域0总质量流率中的比例
                
    %Y_frac0_total = Y_Mg_frac0 + Y_CO_frac0;
    %Y_Mg_frac0 = Y_Mg_frac0 /Y_frac0_total;
    %Y_CO_frac0 = Y_CO_frac0 /Y_frac0_total;
    % 保存当前值用于自适应松弛
    persistent last_iter_count;
    if isempty(last_iter_count)
        last_iter_count = 0;
    end
                    
    % 根据迭代轮次动态调整松弛因子，避免过度震荡
    if simplified_mode && isfield(env, 'iter_count')
        if env.iter_count > last_iter_count + 10
            relaxation = relaxation * 0.9;  % 如果迭代次数增加很多，减小松弛因子
        end
        last_iter_count = env.iter_count;
    end
    
    % 区域1（氧化层外表面到火焰面）参数
    Y_T_surf_nd = X(9);           % 表面温度 无量纲
    B1_nd = X(10);                % 区域1的无量纲B参数
    
    Y_Mg_surf = X(11);            % Mg表面质量分数
    Y_Mg_frac1 = X(12);           % Mg在区域1总质量流率中的比例
    
    Y_CO2_surf = X(13);           % CO2表面质量分数
    Y_CO2_frac1 = X(14);          % CO2在区域1总质量流率中的比例
    
    Y_CO_surf = X(15);            % CO表面质量分数
    Y_CO_frac1 = X(16);           % CO在区域1总质量流率中的比例
    
    Y_MgO_surf = X(17);           % MgO表面质量分数
    Y_MgO_frac1 = X(18);          % MgO在区域1总质量流率中的比例

    % 区域2: 火焰面到远场
    Y_T_flame_nd = X(19);         % 火焰面温度(无量纲)
    B2_nd = X(20);                % 区域2的无量纲B参数
    
    Y_Mg_flame = X(21);           % 火焰面Mg质量分数
    Y_Mg_frac2 = X(22);           % 区域2 Mg在总质量流率中的比例
    
    Y_CO2_flame = X(23);          % 火焰面CO2质量分数
    Y_CO2_frac2 = X(24);          % 区域2 CO2在总质量流率中的比例
    
    Y_CO_flame = X(25);           % 火焰面CO质量分数
    Y_CO_frac2 = X(26);           % 区域2 CO在总质量流率中的比例
    
    Y_MgO_flame = X(27);          % 火焰面MgO质量分数
    Y_MgO_frac2 = X(28);          % 区域2 MgO在总质量流率中的比例



    % 表面反应部分
    K_CO_surface = params.reaction_pre_exponential * exp(-params.reaction_activation_energy / (params.R_u * T_c));
    %M_mix_surface = Y_Mg_region0_surf * mw.Mg.molar_mass + Y_CO_region0_surf * mw.CO.molar_mass;
    M_mix_surface = Y_Mg_core * mw.Mg.molar_mass + Y_CO_core * mw.CO.molar_mass;
    %A_surface = 4 * pi * r_p^2 * porosity;
    A_surface = 4 * pi * r_c^2 ;
    pressure = params.ambient_pressure  / 101325; 
    %m_dot_surface_Mg = A_surface * Y_CO_region0_surf *  pressure / (params.R_u * T_p) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass *  K_CO_surface;
    m_dot_surface_Mg = A_surface * Y_CO_core *  pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass *  K_CO_surface;
    m_dot_surface_CO = m_dot_surface_Mg * mw.CO.molar_mass / mw.Mg.molar_mass;
    m_dot_surface_MgO = m_dot_surface_Mg * mw.MgO.molar_mass / mw.Mg.molar_mass;

    % 计算区域0的质量流率 (氧化层内)
    r_ini = params.initial_diameter/2 - params.initial_oxide_thickness;
    m_dot_region1_total = Pe_m_region1 *(4 * pi * r_c * rho_D_gas) ;
    m_dot_region1_Mg = m_dot_region1_total * Y_Mg_frac1;
    m_dot_region1_CO = m_dot_region1_total * Y_CO_frac1;
    m_dot_region1_MgO = m_dot_region1_total * Y_MgO_frac1;
    m_dot_region1_CO2 = m_dot_region1_total * Y_CO2_frac1;

    m_dot_region0_Mg = m_dot_region1_Mg;
    m_dot_region0_CO = -(abs(m_dot_region1_CO) - m_dot_surface_CO);
    m_dot_region0_total = m_dot_region0_Mg + m_dot_region0_CO;

   
    Y_Mg_frac0_cal = m_dot_region0_Mg/m_dot_region0_total;            % Mg在区域0总质量流率中的比例
    Y_CO_frac0_cal = m_dot_region0_CO/m_dot_region0_total;     
            
    %fprintf('m_dot_region0_CO  %.6e \n',m_dot_region0_CO);
    %fprintf('m_dot_region0_Mg  %.6e \n',m_dot_region0_Mg);

    m_dot_region2_Mg = 0;
    m_dot_region2_CO2 = - m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO2.molar_mass;
    m_dot_region2_CO = m_dot_region1_Mg / mw.Mg.molar_mass *mw.CO.molar_mass - abs(m_dot_region1_CO);
    m_dot_region2_MgO = m_dot_region1_Mg / mw.Mg.molar_mass *mw.MgO.molar_mass - abs(m_dot_region1_MgO);
    m_dot_region2_total = m_dot_region2_Mg + m_dot_region2_CO2+ m_dot_region2_CO + m_dot_region2_MgO;

    %fprintf(' m_dot_region0_Mg: %.3e \n', m_dot_region0_Mg);
    %fprintf(' Y_Mg_frac0: %.3e \n', Y_Mg_frac0);
    n_dot_region0_Mg = abs(m_dot_region0_Mg) / mw.Mg.molar_mass;
    n_dot_region0_CO = abs(m_dot_region0_CO) / mw.CO.molar_mass;

    n_frac_region0 = n_dot_region0_Mg + n_dot_region0_CO;

    n_frac_region0_Mg = n_dot_region0_Mg / n_frac_region0;
    n_frac_region0_CO = n_dot_region0_CO / n_frac_region0;
    % 计算区域1和区域2的总质量流率 (归一化)
    %m_dot_region1_Mg = m_dot_region0_Mg - m_dot_surface_Mg;
    % 关联二区的质量流率占比与质量分数;
    %Y_CO2_surf = 0.0;
    %Y_CO2_frac1 = 0.0;
    %Y_Mg_frac2 = 0.0;
    %Y_MgO_surf = 0 ;
    
    %Y_Mg_flame = 0 ;
    %Y_CO2_flame = 0 ;


    % 计算区域0、1和区域2的物质输运佩克莱数
    Pe_m_region0 = m_dot_region0_total / (4 * pi * r_c * D_ox_gas);
    Pe_m_region1 = Pe_m_region1 ;
    Pe_m_region2 = m_dot_region2_total / (4 * pi * r_c * rho_D_gas);

    Pe_cp_region0 = m_dot_region0_total * cp_gas / (4 * pi * r_c * k_ox_gas);
    %Pe_cp_region0 = m_dot_region0_total  / (4 * pi * r_c * ratio_k_cp_ox);
    Pe_cp_region1 = m_dot_region1_total * cp_gas / (4 * pi * r_c * k_gas);
    %Pe_cp_region1 = m_dot_region1_total / (4 * pi * r_c * ratio_k_cp_gas);
    Pe_cp_region2 = m_dot_region2_total * cp_gas / (4 * pi * r_c * k_gas);
    %Pe_cp_region2 = m_dot_region2_total / (4 * pi * r_c * ratio_k_cp_gas);
    %fprintf('m_dot_region0_total  %.6e \n',m_dot_region0_total);
    %fprintf('m_dot_region1_total  %.6e \n',m_dot_region1_total);
    %fprintf('m_dot_region2_total  %.6e \n',m_dot_region2_total);

    %fprintf('Pe_m_region0  %.6e \n',Pe_m_region0);
   % fprintf('Pe_m_region1  %.6e \n',Pe_m_region1);
   % fprintf('Pe_m_region2  %.6e \n',Pe_m_region2);
    
    
    %fprintf('Pe_cp_region0  %.6e \n',Pe_cp_region0);
    %fprintf('Pe_cp_region1  %.6e \n',Pe_cp_region1);
    %fprintf('Pe_cp_region2  %.6e \n',Pe_cp_region2);
    % 无量纲远场边界
    %  r_inf_nd = 5.0 ;  % 无量纲远场边界 (远场边界是颗粒半径的5倍)
    r_inf = r_p * params.r_inf_particle;
    r_inf_nd = r_inf / r_c;
    %r_inf_nd = 10 ;
    %r_inf_nd = params.r_inf_nd; 
    % 初始化残差向量 (34个方程: 8个新方程 + 原有26个方程)
    F = zeros(32, 1);  % 减少到32个方程
    % --- 组0: 核心表面边界条件 (r = r_c) ---
     
    % BC 2: 表面Mg质量分数由克-克方程在沸点下的情况决定
    A = 4.32e5;
    B = 1.77e4;
    P_Mg_surface = A * exp(-B / T_c);
    V_frac_Mg = P_Mg_surface / pressure;
    V_frac_CO = 1 - V_frac_Mg;
    Y_Mg_theory = V_frac_Mg * mw.Mg.molar_mass / (V_frac_Mg * mw.Mg.molar_mass + V_frac_CO * mw.CO.molar_mass);
    %Y_Mg_surf  =  Y_mg_theory;
    %Y_CO_surf = 1 - Y_Mg_surf;
    
    % BC 3: CO在核心表面反应消耗
   % F(3) = relaxation * (m_dot_CO_reaction + m_dot_region0_CO);  % CO消耗与流量平衡
    
    % --- 组1: 氧化层外表面边界条件 (r = r_p) ---
    
    % 计算区域0在r_p处的组分分布
    exp_arg_region0_Y = Pe_m_region0 * (1 - 1/r_p_nd);
    Y_Mg_region0_surf = Y_Mg_core + (Y_Mg_core - Y_Mg_frac0) * (exp(exp_arg_region0_Y) - 1);
    Y_CO_region0_surf = Y_CO_core + (Y_CO_core - Y_CO_frac0) * (exp(exp_arg_region0_Y) - 1);
    
    % 温度分布
    
    exp_arg_region0_T = Pe_cp_region0 * (1 - 1/r_p_nd);
    T_region0_surf = Y_T_core_nd + (Y_T_core_nd + B0_nd) * (exp(exp_arg_region0_T) - 1);

    F(1) = relaxation * (Y_T_core_nd - 1.0)*10;  % BC 1: 核心表面温度为沸点(无量纲为1.0)
    F(2) = relaxation * (Y_Mg_core - Y_Mg_theory)*5;
    F(3) = relaxation * (Y_T_surf_nd - T_region0_surf)*10;  % 温度连续
    F(4) = relaxation * (Y_Mg_region0_surf - Y_Mg_surf)*10;  % Mg连续
    F(5) = relaxation * (Y_CO_region0_surf - Y_CO_surf)*10;  % CO连续
    F(6) = relaxation * Y_MgO_surf*10 ;  % 氧化层内无MgO气态组分
    F(7) = relaxation * Y_CO2_surf ;  % 氧化层内无CO2气态组分

    %fprintf('  Pe_m_region0: %.3e \n', Pe_m_region0);
    %fprintf('  (Y_Mg_core - Y_Mg_frac0): %.3e \n', (Y_Mg_core - Y_Mg_frac0));
    %fprintf('  Y_Mg_frac0: %.3e \n', Y_Mg_frac0);
    %fprintf('  (1 - 1/r_p_nd): %.3e \n', (1 - 1/r_p_nd));
    %fprintf('  (exp(exp_arg_region0_Y) - 1): %.3e \n', (exp(exp_arg_region0_Y) - 1));
    %fprintf('  Y_Mg_region0_surf: %.3e \n', Y_Mg_region0_surf);
  
    % 计算核心表面能量平衡
    % 蒸发潜热贡献
    H_reac_f = params.reaction_heats.flame_reac_H;  % [J/kg]
    L_v = mw.Mg.L_evap_Mg / mw.Mg.molar_mass;
    Q_evap = m_dot_region1_Mg * L_v;
    Q_flame = m_dot_region1_Mg * abs(H_reac_f);
    %fprintf('  m_dot_region1_Mg: %.3e \n', m_dot_region1_Mg);
    %fprintf('  Q_evap: %.3e, Q_flame: %.3e  \n', Q_evap,Q_flame);
    %fprintf('  ratio: %.3e \n', Q_evap/Q_flame);
    q_conv = params.k_gas / r_p * (T_amb -T_p);
    heat_convection = q_conv * (4 * pi * r_p^2);
    heat_convection = 0 ;


    % 表面反应热
    H_reac_CO = params.reaction_heats.surface_reac_H;
    Q_reac_CO = abs(m_dot_surface_Mg) * (-H_reac_CO);
    %Q_reac_CO = 0 ;
    % 辐射
    emissivity_oxide = params.emissivity ;
    Q_rad_ox = emissivity_oxide * params.sigma * (4 * pi * r_c^2) * ((T_p)^4 - (T_c)^4) ;
    
    % 计算B0参数
    Q_sens_region0 = m_dot_region0_total * cp_gas * T_c;
    A = Q_evap - Q_reac_CO - Q_rad_ox -heat_convection ;
    A_porosity = A  * porosity ; 
    B0 = A_porosity - Q_sens_region0;
    B0_calc_nd = B0 / (m_dot_region0_total * cp_gas* T_c);
    %fprintf('  B0_nd: %.3e \n', B0_nd);
    %fprintf('  B0_calc_nd: %.3e \n', B0_calc_nd);
    %fprintf(' MG表面 Q_cod_region0: %.3e \n', Q_cod_region0);
    F(8) = relaxation * (B0_nd - B0_calc_nd)*0.01;

    % --- 组2: 颗粒表面边界条件 (r = r_p) ---
    % 注意：不再强制表面温度等于沸点，而是通过连续性条件确定
    % F(10) = relaxation * (Y_T_surf_nd - 1.0) * 2;  % BC: 表面温度等于沸点(无量纲为1.0)
    % 放松表面温度约束，用较小的权重确保不会偏离太远

    F(9) = relaxation * (Y_T_surf_nd - T_p_nd)*10;  % 弱约束：表面温度接近初始设定值
    %fprintf(' Y_T_surf_nd: %.3e \n', Y_T_surf_nd*T_c);
    %fprintf(' T_p_nd: %.3e \n', T_p_nd*T_c);
    %F(11) = relaxation * Y_CO_surf;      % BC: 表面CO浓度为零
    exp_arg_T_region1 = Pe_cp_region1 * ( 1/r_p_nd - 1/r_f_nd);

    T_flame_region1 = Y_T_surf_nd+ (Y_T_surf_nd + B1_nd) * (exp(exp_arg_T_region1) - 1);
    T_flame = T_flame_region1 * T_c;
    
    %fprintf('  T_flame: %.3e \n', T_flame);
    % 辐射热贡献
    emissivity = params.emissivity;
    sigma = params.sigma ;

            %fprintf('  > Q_rad = %.6e\n', Q_rad );
    q_conv = params.k_gas / r_p * (T_amb -T_p);
    heat_convection = q_conv * (4 * pi * r_p^2);
    %heat_convection = 0 ;

    m_dot_Mg_1 = m_dot_region1_Mg ;
    H_reac_f = params.reaction_heats.flame_reac_H;  % [J/kg]
    Q_reac_f_to_face = abs(m_dot_Mg_1)  * (-H_reac_f) * k;
    %Q_reac_f_to_face = 0 ;

    H_dep_MgO = params.H_dep_MgO;
    Q_dep_MgO = abs(m_dot_region1_MgO)* (-H_dep_MgO);

    % Q_rad  绝对值
    Q_rad_f = emissivity * sigma * (4 * pi * r_p^2) * (T_flame^4 - T_p^4);
    %Q_rad_f = 0;
    Q_rad_amb = params.emissivity * params.sigma * (4 * pi * r_p^2) * (T_amb^4 - T_p^4);
    Q_rad_amb = 0 ;
    Q_rad = Q_rad_f+ Q_rad_amb ;
    
    H_reac_CO = params.reaction_heats.surface_reac_H;
    Q_reac_CO = abs(m_dot_surface_CO) * (-H_reac_CO);
    Q_reac_CO = 0 ;       


    % 显热贡献
    Q_sens_1 = m_dot_region1_total * cp_gas * T_p;
    %%% 左侧导热 = 颗粒表面（0区域）左侧  导热梯度项 * 孔隙率
    Q_conv_left = (A_porosity + m_dot_region0_total * cp_gas * ( T_p - T_c ) ) ;
    B =(Q_rad  + Q_reac_CO + Q_dep_MgO + heat_convection + Q_reac_f_to_face- ...
            Q_rad_ox - Q_conv_left) ;
            %fprintf(' 颗粒表面 B : %.3e \n', B );
    % 计算B1 (物理量)
    B1 = B  - Q_sens_1;
    %B1 = X1 /  (m_dot_region1_total * cp_gas );
    %fprintf(' 颗粒表面 B : %.3e \n', B );
    % 转换为无量纲B1
    B1_calc_nd = B1 / (m_dot_region1_total * cp_gas* T_c) ;
    %fprintf(' 颗粒表面 Q_conv_right: %.3e: Q_sens_region1:%.3e\n', Q_cond_right, Q_sens_1);
    %fprintf('  B1_nd: %.3e \n', B1_nd);
    %fprintf('  B1_calc_nd: %.3e \n', B1_calc_nd);
    % 约束B1参数
    F(10) = relaxation * (B1_nd - B1_calc_nd)*0.01;  % 原来是F(14)，现在是F(12)


    exp_arg_Y_region1 = Pe_m_region1 * (1/r_p_nd - 1/r_f_nd);
    Y_CO2_region1_flame = Y_CO2_surf + (Y_CO2_surf - Y_CO2_frac1) * (exp(exp_arg_Y_region1) - 1);
    Y_Mg_region1_flame = Y_Mg_surf + (Y_Mg_surf - Y_Mg_frac1) * (exp(exp_arg_Y_region1) - 1);
    Y_CO_region1_flame = Y_CO_surf + (Y_CO_surf - Y_CO_frac1) * (exp(exp_arg_Y_region1) - 1);
    Y_MgO_region1_flame = Y_MgO_surf + (Y_MgO_surf - Y_MgO_frac1) * (exp(exp_arg_Y_region1) - 1);
     %fprintf(' 初始猜测 exp_arg_Y_nd: %.3e \n', exp(exp_arg_Y_nd));
    %fprintf(' 初始猜测 Y_CO2_frac2: %.3e \n', Y_CO2_frac2);s's's
    %fprintf('  Y_Mg_surf: %.3e \n', Y_Mg_surf);
    %fprintf('  Y_Mg_frac1: %.3e \n', Y_Mg_surf);
    %fprintf('  Y_Mg_surf: %.3e \n', Y_Mg_surf);

    % --- 组3: 远场边界条件 (r = r_inf) ---
    
    % 远场温度
    exp_arg_T_far = Pe_cp_region2 * (1/r_f_nd - 1/r_inf_nd);
    %exp_arg_T_far = Pe_cp_region2 * ( 1/r_f_nd);
    exp_arg_T_region1 = Pe_cp_region1 * ( 1/r_p_nd - 1/r_f_nd);
    T_flame_region1 = Y_T_surf_nd+ (Y_T_surf_nd + B1_nd) * (exp(exp_arg_T_region1) - 1);

    Y_T_inf_nd = Y_T_flame_nd + (Y_T_flame_nd + B2_nd) * (exp(exp_arg_T_far) - 1);

    F(11) = relaxation * (Y_T_inf_nd - T_amb_nd);  % 增强远场温度约束
    %fprintf(' Y_T_inf_nd: %.3e \n', Y_T_inf_nd * T_c);
    % 远场组分
    exp_arg_Y_nd = Pe_m_region2 * (1/r_f_nd - 1/r_inf_nd);
    %exp_arg_Y_nd = Pe_m_region2 * (1/r_f_nd );
    Y_CO_inf = Y_CO_flame + (Y_CO_flame - Y_CO_frac2) * (exp(exp_arg_Y_nd) - 1);
    Y_MgO_inf = Y_MgO_flame + (Y_MgO_flame - Y_MgO_frac2) * (exp(exp_arg_Y_nd) - 1);
    Y_Mg_inf = Y_Mg_flame + (Y_Mg_flame - Y_Mg_frac2) * (exp(exp_arg_Y_nd) - 1);
    Y_CO2_inf = Y_CO2_flame + (Y_CO2_flame - Y_CO2_frac2) * (exp(exp_arg_Y_nd) - 1);
    %fprintf(' 初始猜测 exp_arg_Y_nd: %.3e \n', exp(exp_arg_Y_nd));
    %fprintf(' 初始猜测 Y_CO2_frac2: %.3e \n', Y_CO2_frac2);s's's
    %fprintf('  Y_MgO_inf: %.3e \n', Y_MgO_inf);
    %Y_CO_inf=0;
    %Y_MgO_inf=0;
    %Y_Mg_inf=0;
    %Y_CO2_inf=1;
    F(12) = relaxation * Y_CO_inf;       % 增强CO约束
    F(13) = relaxation * Y_MgO_inf;      % 增强MgO约束
    F(14) = relaxation * Y_Mg_inf;      % 强化远场Mg约束
    % F(15) = relaxation * (Y_CO2_inf-1) ;  % 强化远场CO2约束
    F(15) = relaxation * (Y_CO2_inf-1)*20;  % 强化远场CO2约束
    %fprintf('  Y_CO2_inf: %.3e \n', Y_CO2_inf);
    %fprintf('  Y_Mg_flame: %.3e \n', Y_Mg_flame);
    % --- 组4: 火焰面界面条件 (r = r_f) ---
    %fprintf('  T_flame_region1: %.3e \n', T_flame_region1 * T_c);
    %fprintf('  Y_T_flame_nd: %.3e \n', Y_T_flame_nd);
    %fprintf('  T_flame_region1: %.3e \n', T_flame_region1 *T_c);
    %fprintf('  Y_T_flame_nd: %.3e \n', Y_T_flame_nd *T_c);
    F(16) = relaxation * (T_flame_region1-Y_T_flame_nd);  % 增强温度连续性
    F(17) = relaxation * (Y_CO2_region1_flame-Y_CO2_flame)*10; % 增强组分连续
    %fprintf(' 求解过程中一区计算CO2浓度: %.3e \n', Y_CO2_region1_flame);
    %fprintf(' 求解过程中火焰面CO2浓度: %.3e \n', Y_CO2_flame);
    F(18) = relaxation *(Y_Mg_region1_flame-Y_Mg_flame); % 增强组分连续
    F(19) = relaxation * (Y_MgO_region1_flame-Y_MgO_flame); % 增强组分连续
    F(20) = relaxation * (Y_CO_region1_flame-Y_CO_flame); % 增强组分连续
    % 火焰面组分和为1
    F(21) = relaxation * (Y_CO_flame + Y_MgO_flame + Y_Mg_flame + Y_CO2_flame - 1.0); % 增强组分和约束
    F(22) = relaxation * (Y_CO2_region1_flame+Y_CO_region1_flame+Y_Mg_region1_flame+Y_MgO_region1_flame-1)*10; % 增强组分和约束
    % 反应物耗尽
    F(23) = relaxation * Y_Mg_region1_flame*10; % 强化Mg在火焰面耗尽条件
    F(24) = relaxation * Y_CO2_flame *20  ; % 强化火焰面CO2约束
    
    % 火焰面能量守恒
        
    % 物理火焰温度和半径
    T_f = T_flame;
    %fprintf(' 求解过程中右侧火焰面温度: %.3e K\n', Y_MgO_region1_flame);
    %fprintf(' 求解过程中zuo侧火焰面温度: %.3e K\n', T_flame_region1*T_p);
    % 火焰反应热贡献 (物理量)
    %m_dot_Mg_1 = m_dot_region1_total * Y_Mg_frac1;m_dot_region1_Mg
    H_reac_f = params.reaction_heats.flame_reac_H;  % [J/kg]
    %Q_reac_f_to_flame = abs(m_dot_region1_Mg)  * (-H_reac_f) * 0.5;
    Q_reac_f_to_flame = abs(m_dot_region1_Mg)  * (-H_reac_f) * (1-k);
    % 辐射热贡献 (物理量)
    emissivity = params.emissivity;
    sigma = params.sigma;
    Q_rad_f_p = emissivity * sigma * (4 * pi * r_f^2) * (T_f^4 - T_p^4);
    Q_rad_f_amb = emissivity * sigma * (4 * pi * r_f^2) * (T_f^4 - T_amb^4);
    Q_rad_f_total = Q_rad_f_p + Q_rad_f_amb;
    %fprintf(' m_dot_region1_Mg: %.3e \n', m_dot_region1_Mg);
    Q_sens_2 = m_dot_region2_total * cp_gas * T_f;
    Q_cond_left = m_dot_region1_total * cp_gas * (T_f - T_p) + B;
    C = -(Q_reac_f_to_flame- Q_rad_f_total - Q_cond_left);
    %%%%   Q_cond_right_1 = Q_cond_right + m_dot_region2_total * cp_gas * ( T_f - T_p );

    % 计算B2 (物理量)
    B2 = C - Q_sens_2;
    B2_calc_nd = B2 / (m_dot_region2_total * cp_gas* T_c);
    % 转换为无量纲B2
    %B2_calc_nd = B2 / T_p;
    %fprintf(' 火焰面 C : %.3e \n', C);
    % 约束B2参数
    F(25) = relaxation * (B2_nd - B2_calc_nd)*0.01; % 原来是F(29)
    
    %%%%%%远场组分守恒
    F(26) = relaxation * (Y_CO2_inf + Y_CO_inf + Y_Mg_inf + Y_MgO_inf - 1) ; 
    
    % 化学计量通量平衡 (无量纲)
    molar_ratio_Mg = Y_Mg_frac1 * m_dot_region1_total / mw.Mg.molar_mass;
    molar_ratio_CO2 = Y_CO2_frac2 * m_dot_region2_total / mw.CO2.molar_mass;
    molar_ratio_CO_1 = Y_CO_frac1 * m_dot_region1_total / mw.CO.molar_mass;
    molar_ratio_CO_2 = Y_CO_frac2 * m_dot_region2_total / mw.CO.molar_mass;
    molar_ratio_CO_total = molar_ratio_CO_2 - molar_ratio_CO_1;
    
    % 原有方程
    F(27) = relaxation * (molar_ratio_Mg - molar_ratio_CO_total); %
    F(28) = relaxation * (molar_ratio_Mg + molar_ratio_CO2); % 
    F(29) = relaxation * (Y_Mg_region0_surf + Y_CO_region0_surf - 1 )*30; % 
 
    %F(30) = relaxation * (Y_Mg_region0_surf+Y_CO_region0_surf - 1) * 20;  
    F(30) = relaxation * (Y_Mg_core + Y_CO_core - 1)*10 ; % 原来是F(34)
    F(31) = relaxation * (Y_Mg_frac0_cal- Y_Mg_frac0);
    F(32) = relaxation * (Y_CO_frac0_cal- Y_CO_frac0);
    
            
    
end

function [coeffs_0, coeffs_1, coeffs_2] = extract_solution_coeffs(X)
    % 提取无量纲系数（扩展版本）

    % 区域0: 核心表面到氧化层表面
    coeffs_0.Y_T_core = X(3);       % 核心表面温度(无量纲)
    %fprintf('  coeffs_0.Y_T_core: %.3e \n', coeffs_0.Y_T_core);
    coeffs_0.m_dot_B0_cp_ratio = X(4); % 无量纲B0参数
    %fprintf('  coeffs_0.m_dot_B0_cp_ratio: %.3e \n', coeffs_0.m_dot_B0_cp_ratio);
    coeffs_0.Y_Mg_core = X(5);      % Mg核心表面质量分数
    coeffs_0.Y_Mg_frac = X(6);      % Mg在区域0质量流率中的比例
    coeffs_0.Y_CO_core = X(7);      % CO核心表面质量分数
    coeffs_0.Y_CO_frac = X(8);      % CO在区域0质量流率中的比例
    
    % 区域1: 氧化层表面到火焰面
    coeffs_1.Y_T_surf = X(9);       % 表面温度(无量纲)
    %fprintf('  coeffs_1.Y_T_surf: %.3e \n', coeffs_1.Y_T_surf);
    coeffs_1.m_dot_B1_cp_ratio = X(10); % 无量纲B1参数
    %fprintf('  coeffs_1.m_dot_B1_cp_ratio: %.3e \n', coeffs_1.m_dot_B1_cp_ratio);
    coeffs_1.Y_Mg_surf = X(11);     % 表面Mg质量分数
    coeffs_1.Y_Mg_frac = X(12);     % Mg在区域1总质量流率中的比例
    coeffs_1.Y_CO2_surf = X(13);    % 表面CO2质量分数
    coeffs_1.Y_CO2_frac = X(14);    % CO2在区域1总质量流率中的比例
    coeffs_1.Y_CO_surf = X(15);     % 表面CO质量分数
    coeffs_1.Y_CO_frac = X(16);     % CO在区域1总质量流率中的比例
    coeffs_1.Y_MgO_surf = X(17);    % 表面MgO质量分数
    coeffs_1.Y_MgO_frac = X(18);    % MgO在区域1总质量流率中的比例
    
    % 区域2: 火焰面到远场
    coeffs_2.Y_T_flame = X(19);     % 火焰面温度(无量纲)
    %fprintf('  coeffs_2.Y_T_flame: %.3e \n', coeffs_2.Y_T_flame);
    coeffs_2.m_dot_B2_cp_ratio = X(20); % 无量纲B2参数
    %fprintf('  coeffs_2.m_dot_B2_cp_ratio: %.3e \n', coeffs_2.m_dot_B2_cp_ratio);
    coeffs_2.Y_Mg_flame = X(21);    % 火焰面Mg质量分数
    coeffs_2.Y_Mg_frac = X(22);     % 区域2 Mg在总质量流率中的比例
    coeffs_2.Y_CO2_flame = X(23);   % 火焰面CO2质量分数
    coeffs_2.Y_CO2_frac = X(24);    % 区域2 CO2在总质量流率中的比例
    coeffs_2.Y_CO_flame = X(25);    % 火焰面CO质量分数
    coeffs_2.Y_CO_frac = X(26);     % 区域2 CO在总质量流率中的比例
    coeffs_2.Y_MgO_flame = X(27);   % 火焰面MgO质量分数
    coeffs_2.Y_MgO_frac = X(28);    % 区域2 MgO在总质量流率中的比例
   
        % fprintf(' coeffs_0.Y_T_core: %.3e \n', coeffs_0.Y_T_core);
         fprintf(' coeffs_0.Y_Mg_frac: %.3e \n', coeffs_0.Y_Mg_frac);
        % fprintf(' coeffs_1.Y_Mg_surf: %.3e \n', coeffs_1.Y_Mg_surf);
        % fprintf(' coeffs_1.Y_MgO_frac: %.3e \n', coeffs_1.Y_MgO_frac);
        % fprintf(' coeffs_2.Y_T_flame: %.3e \n', coeffs_2.Y_T_flame);
        % fprintf(' coeffs_0.Y_T_core: %.3e \n', coeffs_0.Y_T_core);
        % fprintf(' coeffs_2.Y_Mg_flame: %.3e \n', coeffs_2.Y_Mg_flame);
        % fprintf(' coeffs_2.Y_Mg_frac: %.3e \n', coeffs_2.Y_Mg_frac);
        % fprintf(' coeffs_2.Y_CO2_flame: %.3e \n', coeffs_2.Y_CO2_flame);
        % fprintf(' coeffs_2.Y_CO2_frac: %.3e \n', coeffs_2.Y_CO2_frac);
        % 
end

function value = evaluate_solution(r, coeffs_1, coeffs_2, variable, r_f, Pe_m, Le, r_p)
    % 评估解在指定位置的值，使用一致的无量纲方法
    % 
    % 无量纲解析解形式：
    % 区域1: Y(r) = (Y_surf - Y_frac) * (exp(Pe_m * ((1/r_p_nd) - (1/r_nd))) - 1) + Y_surf
    % 区域2: Y(r) = (Y_flame - Y_frac) * (exp(Pe_m * ((1/r_f_nd) - (1/r_nd))) - 1) + Y_flame
    
    % 计算无量纲径向坐标
    r_nd = r / r_p;
    r_f_nd = r_f / r_p;
    
    % 计算区域1和区域2的总质量流率比例
    m_dot_region1_total_nd = coeffs_1.Y_Mg_frac + coeffs_1.Y_CO2_frac + coeffs_1.Y_CO_frac + coeffs_1.Y_MgO_frac;
    m_dot_region2_total_nd = coeffs_2.Y_Mg_frac + coeffs_2.Y_CO2_frac + coeffs_2.Y_CO_frac + coeffs_2.Y_MgO_frac;
    
    % 计算各区域的佩克莱数
    Pe_m_region1_nd = Pe_m;
    Pe_m_region2_nd = Pe_m * (m_dot_region2_total_nd / m_dot_region1_total_nd);
    
    % 估计热传导佩克莱数
    cp_gas = 1000;  % 简化假设
    k_gas = 0.1;    % 简化假设
    
    % 气体物性的无量纲参数     co2  60 j/mol/k
    cp_gas = 1400;  % 气体比热估计值 [J/(kg·K)]
    % co2  2000K
    k_gas = 0.102682;  % 气体导热系数估计值 [W/(m·K)]
    
    Pe_cp_region1_nd = Pe_m_region1_nd * cp_gas / k_gas;
    Pe_cp_region2_nd = Pe_m_region2_nd * cp_gas / k_gas;
    
    % 确定所在区域并评估解
    if r <= r_f
        % 区域1: Y(r) = (Y_surf - Y_frac) * (exp(Pe_m * ((1/r_p_nd) - (1/r_nd))) - 1) + Y_surf
        
        % 计算exp项
        exp_arg_nd = Pe_m_region1_nd * (1 - 1/r_nd);
        exp_term_nd = safe_exp_with_log(exp_arg_nd) - 1;
        
        if strcmp(variable, 'T')
            % 区域1温度
            exp_arg_T = Pe_cp_region1_nd * (1/r_nd - 1);
            log_exp_term_T = exp_arg_T;  % 存储对数值

            % 当需要实际值时才转换
            if log_exp_term_T > 30
                % 大值使用对数算术
                T_nd(i) = coeffs_1.Y_T_surf + (coeffs_1.Y_T_surf - coeffs_1.m_dot_B1_cp_ratio) * exp(log_exp_term_T);
            else
                % 正常范围使用标准计算
                exp_term_T = safe_exp_with_log(exp_arg_T);
                T_nd(i) = coeffs_1.Y_T_surf + (coeffs_1.Y_T_surf - coeffs_1.m_dot_B1_cp_ratio) * exp_term_T;
            end
        elseif strcmp(variable, 'Mg')
            value = (coeffs_1.Y_Mg_surf - coeffs_1.Y_Mg_frac) * exp_term_nd + coeffs_1.Y_Mg_surf;
        elseif strcmp(variable, 'CO2')
            value = (coeffs_1.Y_CO2_surf - coeffs_1.Y_CO2_frac) * exp_term_nd + coeffs_1.Y_CO2_surf;
        elseif strcmp(variable, 'CO')
            value = (coeffs_1.Y_CO_surf - coeffs_1.Y_CO_frac) * exp_term_nd + coeffs_1.Y_CO_surf;
        elseif strcmp(variable, 'MgO')
            value = (coeffs_1.Y_MgO_surf - coeffs_1.Y_MgO_frac) * exp_term_nd + coeffs_1.Y_MgO_surf;
        end
    else
        % 区域2: Y(r) = (Y_flame - Y_frac) * (exp(Pe_m * ((1/r_f_nd) - (1/r_nd))) - 1) + Y_flame
        
        % 计算exp项
        exp_arg_nd = Pe_m_region2_nd * (1/r_f_nd - 1/r_nd);
        exp_term_nd = safe_exp_with_log(exp_arg_nd) - 1;
        
        if strcmp(variable, 'T')
            % 区域2温度
            value = (coeffs_2.Y_T_flame - coeffs_2.m_dot_B2_cp_ratio) * exp_term_nd + coeffs_2.Y_T_flame;
        elseif strcmp(variable, 'Mg')
            value = (coeffs_2.Y_Mg_flame - coeffs_2.Y_Mg_frac) * exp_term_nd + coeffs_2.Y_Mg_flame;
        elseif strcmp(variable, 'CO2')
            value = (coeffs_2.Y_CO2_flame - coeffs_2.Y_CO2_frac) * exp_term_nd + coeffs_2.Y_CO2_flame;
        elseif strcmp(variable, 'CO')
            value = (coeffs_2.Y_CO_flame - coeffs_2.Y_CO_frac) * exp_term_nd + coeffs_2.Y_CO_flame;
        elseif strcmp(variable, 'MgO')
            value = (coeffs_2.Y_MgO_flame - coeffs_2.Y_MgO_frac) * exp_term_nd + coeffs_2.Y_MgO_flame;
        end
    end
end

% 使用对数表示处理指数计算
function exp_term = safe_exp_with_log(exp_arg)
    if exp_arg > 700
        % 直接使用对数形式计算
        log_result = exp_arg;
        exp_term = exp(log_result) - 1;
    elseif exp_arg < -30
        % 对于非常小的指数参数，使用泰勒展开近似
        exp_term = exp_arg;  % 当exp(x)≈1+x时，exp(x)-1≈x
    else
        % 正常范围内直接计算
        exp_term = exp(exp_arg) - 1;
    end
end

function rate_info = package_rate_info(m_dot_total_region1, r_f, T_f, params, r_p, r_c, coeffs_0, coeffs_1, coeffs_2,T_c)
    % 将计算结果封装为与原BVP方法兼容的结构
    % 所有输入均为物理量





    
    mw = params.materials;
    porosity = params.material_properties.oxide_porosity;
    % 基本信息
    rate_info.m_dot_region1_total = m_dot_total_region1;             % 总质量流率(蒸发速率)
    rate_info.r_flame = r_f;             % 火焰位置
    rate_info.T_flame = T_f;             % 火焰温度
    rate_info.r_f_to_r_p_ratio = r_f / r_p; % 火焰半径与颗粒半径比值



    m_dot_region1_total = m_dot_total_region1;
    m_dot_region1_Mg = coeffs_1.Y_Mg_frac * m_dot_region1_total;
    m_dot_region1_CO = coeffs_1.Y_CO_frac * m_dot_region1_total;
    m_dot_region1_MgO = coeffs_1.Y_MgO_frac * m_dot_region1_total;
    m_dot_region1_CO2 = coeffs_1.Y_CO2_frac * m_dot_region1_total;
    
    % 表面反应部分处理
    
    K_CO_surface = params.reaction_pre_exponential * exp(-params.reaction_activation_energy / (params.R_u * T_c));
    Y_Mg_region0_core = coeffs_0.Y_Mg_core;
    Y_CO_region0_core = coeffs_0.Y_CO_core;
    M_mix_surface = Y_Mg_region0_core * mw.Mg.molar_mass + Y_CO_region0_core * mw.CO.molar_mass;
    A_surface = 4 * pi * r_c^2 ;
    pressure = params.ambient_pressure  / 101325; 
    m_dot_surface_Mg = A_surface * Y_CO_region0_core * pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass ...
                        *  K_CO_surface;
                        
    m_dot_surface_CO = m_dot_surface_Mg * mw.CO.molar_mass / mw.Mg.molar_mass;
    m_dot_surface_MgO = m_dot_surface_Mg * mw.MgO.molar_mass / mw.Mg.molar_mass;



   
    
    % 计算各组分物理质量流率
    m_dot_region0_Mg = m_dot_region1_Mg ;
    m_dot_region0_CO = -(abs(m_dot_region1_CO) - m_dot_surface_CO);
    m_dot_region0_total = m_dot_region0_Mg + m_dot_region0_CO;
    
    % 反应产物质量
    m_dot_MgO_reaction_surface = m_dot_surface_MgO;
    m_dot_Mg_reaction_surface = m_dot_surface_Mg;
    m_dot_C_reaction_surface = m_dot_surface_CO / mw.CO.molar_mass * mw.C.molar_mass;
    
    % 区域2各组分的质量流率 - 2区总质量流率用反应算
    m_dot_region2_CO2 =  - m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO2.molar_mass;
    m_dot_region2_Mg = 0;
    m_dot_region2_CO = m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO.molar_mass - abs(m_dot_region1_CO);
    m_dot_region2_MgO = m_dot_region1_Mg / mw.Mg.molar_mass * mw.MgO.molar_mass - abs(m_dot_region1_MgO);
    m_dot_region2_total = m_dot_region2_CO2 + m_dot_region2_Mg + m_dot_region2_CO + m_dot_region2_MgO;
    %fprintf('m_dot_region0_total: %.6e ,m_dot_region1_total: %.6e,m_dot_region2_total: %.6e\n', m_dot_region0_total,m_dot_region1_total,m_dot_region2_total);
    fprintf('m_dot_region0_Mg: %.6e ,m_dot_region1_Mg: %.6e,m_dot_region2_Mg: %.6e ,m_dot_Mg_reaction_surface: %.6e\n', m_dot_region0_Mg,m_dot_region1_Mg,m_dot_region2_Mg,m_dot_Mg_reaction_surface);

    % 输出各组分的总体质量流率 (用于ODE求解)  带符号
    rate_info.dmdt_mg = - m_dot_region1_Mg - m_dot_surface_Mg;    % 
    rate_info.dmdt_mg_total =  - m_dot_region1_Mg - m_dot_surface_Mg;
    rate_info.dmdt_mg_surface_reaction = m_dot_surface_Mg;
    rate_info.dmdt_MgO_surface_reaction = m_dot_MgO_reaction_surface;
    rate_info.dmdt_C_surface_reaction = m_dot_C_reaction_surface;
    rate_info.dmdt_Mg_surface_reaction = m_dot_Mg_reaction_surface;

    rate_info.dmdt_CO2 = - m_dot_region2_CO2;  % CO2质量流率
    rate_info.dmdt_CO = m_dot_region1_CO;    % CO质量流率 (使用区域1的值)

    rate_info.dmdt_mgo = - m_dot_region1_MgO + m_dot_MgO_reaction_surface;  % MgO质量流率
    rate_info.dmdt_mgo_in = - m_dot_region1_MgO;
    rate_info.dmdt_c = m_dot_C_reaction_surface;      % C质量流率
    
    % 保存更多详细信息
    rate_info.coeffs_region0 = coeffs_0;   % 添加区域0的系数
    rate_info.coeffs_region1 = coeffs_1;
    rate_info.coeffs_region2 = coeffs_2;
    rate_info.m_dot_region1_total = m_dot_region1_total;
    rate_info.m_dot_region2_total = m_dot_region2_total;
    rate_info.m_dot_Mg_reaction = m_dot_Mg_reaction_surface;
    
    % 保存各区域的组分质量流率

    rate_info.m_dot_Mg_region0 = m_dot_region0_Mg;
    rate_info.m_dot_CO_region0 = m_dot_region0_CO;
    rate_info.m_dot_region0_total = m_dot_region0_total;

    rate_info.m_dot_Mg_region1 = m_dot_region1_Mg;
    rate_info.m_dot_CO2_region1 = m_dot_region1_CO2;
    rate_info.m_dot_CO_region1 = m_dot_region1_CO;
    rate_info.m_dot_MgO_region1 = m_dot_region1_MgO;
    
    rate_info.m_dot_Mg_region2 = m_dot_region2_Mg;
    rate_info.m_dot_CO2_region2 = m_dot_region2_CO2;
    rate_info.m_dot_CO_region2 = m_dot_region2_CO;
    rate_info.m_dot_MgO_region2 = m_dot_region2_MgO;

    % 添加区域0相关信息
    %m_dot_region0_total = coeffs_0.Y_Mg_frac * m_dot_region0_total + coeffs_0.Y_CO_frac * m_dot_region0_total;
    %m_dot_region0_Mg = coeffs_0.Y_Mg_frac * m_dot_region0_total;
    %m_dot_region0_CO = coeffs_0.Y_CO_frac * m_dot_region0_total;
    
    % 更新输出结构体
    %rate_info.m_dot_Mg_core = m_dot_region0_Mg;
    %rate_info.m_dot_CO_core = m_dot_region0_CO;
    %rate_info.m_dot_CO_reaction = m_dot_region1_CO;
end

function [X_nd] = convert_to_dimensionless(X, r_p, r_c, T_p, T_c, m_dot_total)
    % 将物理量转换为无量纲量 - 基于核心半径r_c和核心温度T_c的无量纲化
    % X - 原始物理量向量
    % r_p - 颗粒半径
    % r_c - 核心半径
    % T_p - 氧化层表面温度
    % T_c - 金属核心温度
    % m_dot_total - 总质量流率
    
    % 初始化无量纲向量
    X_nd = zeros(size(X));
    
    % 全局参数无量纲化
    X_nd(1) = X(1) / r_c;         % 火焰半径无量纲化 r_f_nd = r_f/r_c
    X_nd(2) = X(2);               % 佩克莱数(已是无量纲)
    
    % 温度相关参数无量纲化
    % 假设：X中某些元素可能代表温度
    % 例如，如果X(19)表示火焰温度，则应该将其无量纲化为X_nd(19) = X(19)/T_c
    
    % 所有其他参数都已经是无量纲形式，直接复制
    for i = 3:length(X)
        X_nd(i) = X(i);
    end
end

function [X] = convert_from_dimensionless(X_nd, r_p, r_c, T_p, T_c, m_dot_total)
    % 将无量纲量转换回物理量 - 基于核心半径r_c和核心温度T_c的无量纲化
    % X_nd - 无量纲向量
    % r_p - 颗粒半径
    % r_c - 核心半径
    % T_p - 氧化层表面温度
    % T_c - 金属核心温度
    % m_dot_total - 总质量流率
    
    % 初始化物理量向量
    X = zeros(size(X_nd));
    
    % 全局参数转回物理量
    X(1) = X_nd(1);         % 火焰半径 r_f = r_f_nd * r_c
    X(2) = X_nd(2);               % 佩克莱数保持不变
    
    % 温度相关参数从无量纲化转回
    % 例如，如果X_nd(19)表示无量纲火焰温度，则应该转换为X(19) = X_nd(19)*T_c
    
    % 所有其他参数都已经是无量纲形式，直接复制
    for i = 3:length(X_nd)
        X(i) = X_nd(i);
    end
end

    function visualize_initial_guess(X0, r_p, r_f, r_inf, T_p, T_amb, params, physicalModel)
        % 可视化初始猜测的组分和温度分布
        % 输入:
        %   X0 - 初始猜测向量 (无量纲形式)
        %   r_p - 粒子半径
        %   r_f - 火焰半径
        %   r_inf - 远场边界
        %   T_p - 表面温度
        %   T_amb - 环境温度
        
        % 打印重要系数信息，帮助调试
        fprintf('===================== 初始猜测分析 =====================\n');
        fprintf('温度参数(无量纲):\n');
        fprintf('  Y_T_surf = %.6f, B1 = %.6f\n', X0(3), X0(4));
        fprintf('  Y_T_flame = %.6f, B2 = %.6f\n', X0(13), X0(14));
        
        fprintf('组分参数(无量纲):\n');
        fprintf('  Y_Mg_surf = %.6f, Y_Mg_frac1 = %.6f\n', X0(5), X0(6));
        fprintf('  Y_CO2_surf = %.6f, Y_CO2_frac1 = %.6f\n', X0(7), X0(8));
        fprintf('  Y_CO_surf = %.6f, Y_CO_frac1 = %.6f\n', X0(9), X0(10));
        fprintf('  Y_MgO_surf = %.6f, Y_MgO_frac1 = %.6f\n', X0(11), X0(12));
        fprintf('  Y_Mg_flame = %.6f, Y_Mg_frac2 = %.6f\n', X0(15), X0(16));
        fprintf('  Y_CO2_flame = %.6f, Y_CO2_frac2 = %.6f\n', X0(17), X0(18));
        fprintf('  Y_CO_flame = %.6f, Y_CO_frac2 = %.6f\n', X0(19), X0(20));
        fprintf('  Y_MgO_flame = %.6f, Y_MgO_frac2 = %.6f\n', X0(21), X0(22));
        fprintf('=======================================================\n');
        mw = params.materials;
        % 提取系数
        r_f_nd = X0(1);  % 无量纲火焰半径
        r_f = r_f_nd * r_p;
        Pe_m = X0(2);    % 佩克莱数
        
        % 温度参数 (无量纲)
        Y_T_surf_nd = X0(3);
        B1_nd = X0(4);
        Y_T_flame_nd = X0(13);
        B2_nd = X0(14);
        
        % 区域1组分参数 (无量纲)
        Y_Mg_surf = X0(5);   Y_Mg_frac1 = X0(6);
        Y_CO2_surf = X0(7);  Y_CO2_frac1 = X0(8);
        Y_CO_surf = X0(9);   Y_CO_frac1 = X0(10);
        Y_MgO_surf = X0(11); Y_MgO_frac1 = X0(12);
        
        % 区域2组分参数 (无量纲)
        Y_Mg_flame = X0(15);  Y_Mg_frac2 = X0(16);
        Y_CO2_flame = X0(17); Y_CO2_frac2 = X0(18);
        Y_CO_flame = X0(19);  Y_CO_frac2 = X0(20);
        Y_MgO_flame = X0(21); Y_MgO_frac2 = X0(22);
        
        % 计算区域1和区域2的总质量流率比例
        m_dot_region1_total = Pe_m * 4 * pi * params.rho_D_gas * r_p;
        m_dot_region1_mg = m_dot_region1_total * Y_Mg_frac1;
        m_dot_region2_total = m_dot_region1_mg / mw.Mg.molar_mass * mw.CO2.molar_mass;
        Pe_m_region1_nd = Pe_m;
        Pe_m_region2_nd = m_dot_region2_total / (4 * pi * r_p * params.rho_D_gas);
        m_dot_region2_total_nd = Y_Mg_frac2 + Y_CO2_frac2 + Y_CO_frac2 + Y_MgO_frac2;
        

        % 简化的热传导参数
        cp_gas = physicalModel.get_gas_properties(T_p).Cp_gas;  % 气体比热估计值 [J/(kg·K)]
        k_gas = physicalModel.get_gas_properties(T_p).k_gas;  % 气体导热系数估计值 [W/(m·K)]
            
        % 气体物性的无量纲参数     co2  60 j/mol/k
        cp_gas = 1400;  % 气体比热估计值 [J/(kg·K)]
        % co2  2000K
        k_gas = 0.102682;  % 气体导热系数估计值 [W/(m·K)]
        
        Pe_cp_region1_nd= m_dot_region1_total * cp_gas / (4 * pi * r_p * k_gas);
        Pe_cp_region2_nd = m_dot_region2_total * cp_gas / (4 * pi * r_f * k_gas);
        
        % 生成径向网格 (无量纲)
        r_f_nd = r_f / r_p;
        r_inf_nd = r_inf / r_p;
        r_pts1_nd = linspace(1.0, r_f_nd, 200);       % 区域1增加到200个点
        r_pts2_nd = linspace(r_f_nd, r_inf_nd, 200);  % 区域2增加到200个点
        r_pts_nd = [r_pts1_nd, r_pts2_nd(2:end)];        % 合并 (无量纲)
        
        % 物理径向网格
        r_pts = r_pts_nd * r_p;
        
        % 初始化存储数组
        n_pts = length(r_pts_nd);
        T_nd = zeros(n_pts, 1);    % 无量纲温度
        Y_Mg = zeros(n_pts, 1);    % Mg质量分数
        Y_CO2 = zeros(n_pts, 1);   % CO2质量分数
        Y_CO = zeros(n_pts, 1);    % CO质量分数
        Y_MgO = zeros(n_pts, 1);   % MgO质量分数
        
        % 计算各点的温度和组分质量分数 (全部使用无量纲计算)
        for i = 1:n_pts
            r_nd_i = r_pts_nd(i);
            
            % 判断区域
            if r_nd_i <= r_f_nd  % 区域1: 颗粒表面到火焰面
                % 温度 (无量纲)
                exp_arg_region1_nd_T = Pe_cp_region1_nd * (1/r_nd_i - 1 );
                exp_term_region1_nd_T = exp(exp_arg_region1_nd_T) - 1;
                T_nd(i) = (Y_T_surf_nd + B1_nd) * exp_term_region1_nd_T + Y_T_surf_nd;
                
                % 组分
                exp_arg_region1_nd_Y = Pe_m_region1_nd * (1/r_nd_i - 1 );
                exp_term_region1_nd_Y = exp(exp_arg_region1_nd_Y) - 1;
                Y_Mg(i) = (Y_Mg_surf - Y_Mg_frac1) * exp_term_region1_nd_Y + Y_Mg_surf;
                Y_CO2(i) = (Y_CO2_surf - Y_CO2_frac1) * exp_term_region1_nd_Y + Y_CO2_surf;
                Y_CO(i) = (Y_CO_surf - Y_CO_frac1) * exp_term_region1_nd_Y + Y_CO_surf;
                Y_MgO(i) = (Y_MgO_surf - Y_MgO_frac1) * exp_term_region1_nd_Y + Y_MgO_surf;
                
            else  % 区域2: 火焰面到远场
                % 温度 (无量纲)
                exp_arg_region2_nd_T = Pe_cp_region2_nd * (1/r_nd_i - 1/r_f_nd);
                exp_term_region2_nd_T = exp(exp_arg_region2_nd_T) - 1;
                T_nd(i) = (Y_T_flame_nd + B2_nd) * exp_term_region2_nd_T + Y_T_flame_nd;
                
                % 组分
                exp_arg_region2_nd_Y = Pe_m_region2_nd * (1/r_nd_i - 1/r_f_nd);
                exp_term_region2_nd_Y = exp(exp_arg_region2_nd_Y) - 1;
                Y_Mg(i) = (Y_Mg_flame - Y_Mg_frac2) * exp_term_region2_nd_Y + Y_Mg_flame;
                Y_CO2(i) = (Y_CO2_flame - Y_CO2_frac2) * exp_term_region2_nd_Y + Y_CO2_flame;
                Y_CO(i) = (Y_CO_flame - Y_CO_frac2) * exp_term_region2_nd_Y + Y_CO_flame;
                Y_MgO(i) = (Y_MgO_flame - Y_MgO_frac2) * exp_term_region2_nd_Y + Y_MgO_flame;
            end
            
            % 确保组分质量分数在[0,1]范围内
            Y_Mg(i) = max(0, min(1, Y_Mg(i)));
            Y_CO2(i) = max(0, min(1, Y_CO2(i)));
            Y_CO(i) = max(0, min(1, Y_CO(i)));
            Y_MgO(i) = max(0, min(1, Y_MgO(i)));
        end
        
        % 物理温度 (从无量纲转换)
        T = T_nd * T_p;
        
        % 创建一个新图形窗口
        figure('Name', '初始猜测分布', 'Position', [100, 100, 1200, 800]);
        
        % 温度分布(顶部左图)
        subplot(2, 2, 1);
        % 在绘图前对数据进行样条插值平滑
        r_fine = linspace(min(r_pts), max(r_pts), 500);
        T_fine = interp1(r_pts, T, r_fine, 'spline');
        plot(r_fine*1e6, T_fine, 'k-', 'LineWidth', 2);
        hold on;
        plot([r_f, r_f]*1e6, [min(T), max(T)], 'r--', 'LineWidth', 1.5);  % 火焰面位置
        text(r_f*1e6, max(T)*0.9, '火焰面', 'Color', 'r', 'FontWeight', 'bold');
        
        % 添加标注点
        plot(r_c*1e6, Y_T_core_nd*T_c, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '核心表面温度');
        plot(r_p*1e6, T_p, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', '氧化层表面温度');
        plot(r_f*1e6, coeffs_2.Y_T_flame*T_c, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '火焰温度');
        plot(r_boundary*1e6, T_amb, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', '环境温度');
        
        hold off;
        xlabel('径向位置 (μm)');
        ylabel('温度 (K)');
        title('初始猜测: 温度分布');
        legend('Location', 'best');
        grid on;
        
        % 无量纲温度(顶部右图)
        subplot(2, 2, 2);
        plot(r_pts_nd, T_nd, 'k-', 'LineWidth', 2);
        hold on;
        plot([r_f_nd, r_f_nd], [0, max(T_nd)*1.1], 'r--', 'LineWidth', 1.5);
        text(r_f_nd, max(T_nd)*0.9, '火焰面', 'Color', 'r', 'FontWeight', 'bold');
        hold off;
        xlabel('无量纲径向位置 (r/r_p)');
        ylabel('无量纲温度 (T/T_p)');
        title('初始猜测: 无量纲温度分布');
        grid on;
        
        % 组分分布(底部左图)
        subplot(2, 2, 3);
        plot(r_pts*1e6, Y_Mg, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mg');
        hold on;
        plot(r_pts*1e6, Y_CO2, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CO2');
        plot(r_pts*1e6, Y_CO, 'g-', 'LineWidth', 1.5, 'DisplayName', 'CO');
        plot(r_pts*1e6, Y_MgO, 'm-', 'LineWidth', 1.5, 'DisplayName', 'MgO');
        plot([r_f, r_f]*1e6, [0, 1], 'k--', 'LineWidth', 1.5);  % 火焰面位置
        
        % 检查组分和是否接近1
        Y_sum = Y_Mg + Y_CO2 + Y_CO + Y_MgO;
        plot(r_pts*1e6, Y_sum, 'k:', 'LineWidth', 2, 'DisplayName', '组分和');
        
        xlabel('径向位置 (μm)');
        ylabel('质量分数');
        title('初始猜测: 组分分布');
        legend('Location', 'best');
        ylim([0, 1.05]);
        grid on;
        
        % 保存图像
        saveas(gcf, 'initial_guess_visualization.png');
        fprintf('初始猜测可视化已保存为: initial_guess_visualization.png\n');
        
        % 检查组分和偏离1的程度
        max_deviation = max(abs(Y_sum - 1));
        fprintf('组分和最大偏差: %.4f\n', max_deviation);
        
        % 检查在火焰面的连续性
        r_f_idx = find(r_pts >= r_f, 1);
        if ~isempty(r_f_idx) && r_f_idx > 1
            fprintf('火焰面处各物理量:\n');
            fprintf('  左侧温度: %.2f K, 右侧温度: %.2f K, 差异: %.2f K\n', T(r_f_idx-1), T(r_f_idx), T(r_f_idx)-T(r_f_idx-1));
            fprintf('  左侧组分和: %.4f, 右侧组分和: %.4f\n', Y_sum(r_f_idx-1), Y_sum(r_f_idx));
            fprintf('  组分质量分数跳变(左-右):\n');
            fprintf('    Mg: %.4f, CO2: %.4f\n', Y_Mg(r_f_idx-1)-Y_Mg(r_f_idx), Y_CO2(r_f_idx-1)-Y_CO2(r_f_idx));
            fprintf('    CO: %.4f, MgO: %.4f\n', Y_CO(r_f_idx-1)-Y_CO(r_f_idx), Y_MgO(r_f_idx-1)-Y_MgO(r_f_idx));
        end
    end
    function analyze_initial_guess(X0, r_p)
        % 分析初始猜测的合理性
        r_f = X0(1)*r_p;
        Pe_m = X0(2);
        
        fprintf('初始猜测分析：\n');
       
        % 检查火焰半径是否合理
        if r_f <= r_p
            fprintf('  警告：火焰半径小于等于颗粒半径！r_f/r_p = %.3f\n', r_f/r_p);
        elseif r_f > 5*r_p
            fprintf('  警告：火焰半径过大！r_f/r_p = %.3f\n', r_f/r_p);
        else
            fprintf('  火焰半径似乎合理: r_f/r_p = %.3f\n', r_f/r_p);
        end
        
        % 检查佩克莱数是否合理
        if Pe_m < 0
            fprintf('  警告：佩克莱数为负！Pe_m = %.3e\n', Pe_m);
        elseif Pe_m < 1e-6
            fprintf('  警告：佩克莱数过小！Pe_m = %.3e\n', Pe_m);
        elseif Pe_m > 100
            fprintf('  警告：佩克莱数过大！Pe_m = %.3e\n', Pe_m);
        else
            fprintf('  佩克莱数似乎合理: Pe_m = %.3e\n', Pe_m);
        end
        
        % 注释掉可视化初始猜测分布代码
        % figure;
        % bar(X0);
        % title('初始猜测值分布');
        % xlabel('参数编号'); ylabel('值');
        % saveas(gcf, 'initial_guess_analysis.png');
    end

    function stop = monitor_progress(x, optimValues, state, solve_env)
        % 监控求解过程
        stop = false;
        
        persistent residual_history;
        persistent current_stage;
        persistent has_printed_header;
        monitor_geniration = 500 ;
        if strcmp(state, 'init')
            % 初始化残差历史
            if isempty(current_stage) || current_stage ~= solve_env.stage_number
                current_stage = solve_env.stage_number;
                residual_history = struct('iter', [], 'total', [], ...
                    'surface_bc', [], 'farfield_bc', [], 'flame_interface', [], ...
                    'chemistry', [], 'energy', []);
                has_printed_header = false; % 重置标题打印状态
            end
            
            if isfield(solve_env, 'iter_count')
                solve_env.iter_count = 0;  % 初始化求解环境的迭代计数
            end
            return;
        end
        
        if strcmp(state, 'iter')
            % 更新迭代计数
            if isfield(solve_env, 'iter_count')
                solve_env.iter_count = solve_env.iter_count + 1;
            end
            
            % 检查参数是否合理
            r_f = x(1);
            Pe_m = x(2);
            
            % 计算残差
            F = equations_system(x, solve_env);
            max_residual = max(abs(F));
            [max_val, max_idx] = max(abs(F));
            
            % 标识当前阶段
            stage_name = '';
            if solve_env.stage_number == 1
                stage_name = '简化模型';
            elseif solve_env.stage_number == 2
                stage_name = '中度简化模型';
            else
                stage_name = '完整模型';
            end
            
            % 检查残差是否增大
            persistent last_max_residual;
            persistent residual_growth_count;
            
            if isempty(last_max_residual)
                last_max_residual = max_residual;
                residual_growth_count = 0;
            end
            
            if optimValues.iteration > 5 && max_residual > last_max_residual * 1.5  % 残差增大超过50%
                residual_growth_count = residual_growth_count + 1;
                fprintf('  警告: 残差显著增大 (%.3e -> %.3e), 增大次数: %d\n', ...
                        last_max_residual, max_residual, residual_growth_count);
                
                if residual_growth_count >= 3  % 连续3次残差增大
                    fprintf('  检测到残差持续增大，终止计算\n');
                    stop = true;
                    return;
                end
            else
                residual_growth_count = 0;  % 重置计数器
            end
            
            last_max_residual = max_residual;
            
            % 记录残差历史
            res_surface = norm(F(1:5));          % 表面边界条件残差
            res_farfield = norm(F(6:10));        % 远场边界条件残差
            res_flame = norm(F(11:16));          % 火焰界面条件残差
            res_energy = norm(F(17));            % 能量平衡残差
            res_chemistry = norm(F(18:22));      % 化学计量残差
            res_total = norm(F);                 % 总残差
            
            residual_history.iter(end+1) = optimValues.iteration;
            residual_history.total(end+1) = res_total;
            residual_history.surface_bc(end+1) = res_surface;
            residual_history.farfield_bc(end+1) = res_farfield;
            residual_history.flame_interface(end+1) = res_flame;
            residual_history.energy(end+1) = res_energy;
            residual_history.chemistry(end+1) = res_chemistry;
        

         %visualize_residual = obj.params.visualize_residual
         %if visualize_residual
            % 输出残差信息
            if mod(optimValues.iteration,monitor_geniration) == 0 || optimValues.iteration == 1
                % 第一行：阶段信息、迭代次数和最大残差
                fprintf('阶段%d(%s) - 迭代%d: 最大残差 = %.3e (方程%d)\n', solve_env.stage_number, stage_name, optimValues.iteration, max_val, max_idx);
                
                % 只在第一次迭代时打印标题行
                if ~has_printed_header || optimValues.iteration == 0
                    % 首先输出方程序号行，确保对齐
                    fprintf('%12s', ' '); % 为了对齐，左侧添加空格
                    for i = 1:length(F)
                        fprintf('%12d', i);
                    end
                    fprintf('\n');
                    
                 
                    
                    has_printed_header = true; % 标记已打印标题
                end
                
                % 输出所有方程的残差值
                fprintf('%12d', optimValues.iteration);
                for i = 1:length(F)
                    fprintf('%12.3e', F(i));
                end
                fprintf('\n');
            end
         %end
        end
        
        if strcmp(state, 'done')
            % 保存残差数据
            save(['residual_history_stage', num2str(current_stage), '.mat'], 'residual_history');
        end
    end

    % 在文件末尾添加解质量评估函数
    function [is_acceptable, quality_score] = assess_solution_quality(X, fval, solve_env)
        % 使用最大残差(而非总残差)
        max_residual = max(abs(fval));
        [max_val, max_idx] = max(abs(fval));
        % 提取关键物理量
        r_f_nd = X(1);
        Pe_m = X(2);
        Y_T_surf = X(3);
        Y_T_flame = X(13);
        T_p = solve_env.T_p;
        T_c = solve_env.T_c;
        r_c = solve_env.r_c;
        r_p = solve_env.r_p;
        r_f = r_f_nd * r_c;
        r_f_nd = r_f / r_p;
        %fprintf('T_surf: %.5e\n', Y_T_surf*1366);
        %fprintf('T_flame: %.5e\n', Y_T_flame*1366);
        %fprintf('表面温差: %.5f K\n', (Y_T_surf-1)*1366);
        fprintf('r_f_nd(相对于颗粒表面的): %.5e\n', r_f_nd);
        % 基本物理合理性检查
        physics_ok = ( r_f_nd > 1.0 )&& (r_f_nd < 10.0);
        %physics_ok =  (r_f_nd < 10.0);
        
        % 残差阈值 - 根据阶段调整
        if solve_env.stage_number == 1
            residual_threshold = 0.5;
        elseif solve_env.stage_number == 2
            residual_threshold = 0.5;
        else
            residual_threshold = 0.5;
        end
        
        % 残差合格检查 - 使用最大残差
        residual_ok = max_residual < residual_threshold;
        %fprintf('max_residual: %.5e\n', max_residual);
        % 计算质量分数 (0-100)
        quality_score = min(100, max(0, 100 * (1 - max_residual/residual_threshold)));
        
        % 综合判断 - 必须同时满足物理合理性和残差要求
        is_acceptable = physics_ok && residual_ok;
        
        % 打印评估结果
        fprintf('  解质量评估: 最大残差=%.3e, 物理合理=%d, 残差合格=%d  ,最大残差 = %.3e (方程%d)\n', ...
            max_residual, physics_ok, residual_ok,max_val,max_idx);
        fprintf('  质量分数: %.1f/100, 是否可接受: %d\n', quality_score, is_acceptable);
    end

    function rate_info = create_failed_result_struct()
        rate_info = struct();
        rate_info.success = false;
        % 添加所有必要字段
        rate_info.dmdt_mg = 0;
        rate_info.dmdt_mgo = 0;
        rate_info.dmdt_c = 0;
        rate_info.m_dot = 0;
        rate_info.r_flame = 0;
        rate_info.T_flame = 0;
        rate_info.r_f_to_r_p_ratio = 0;
        rate_info.m_dot_CO = 0;
    end

    function plot_radial_distributions(r_p, r_f, r_inf, Pe_m_region1, coeffs_0, coeffs_1, coeffs_2, T_p, T_c, params, physicalModel,pState,cp_gas...
                        ,k_gas,rho_D_gas,porosity,tortuosity,D_ox_gas,ratio_k_cp_ox,ratio_k_cp_gas,k_ox,k_ox_gas)
        % 使用持久变量保存图形句柄
        persistent fig_handle;
        
        % 检查是否已存在图形或已被关闭
        if isempty(fig_handle) || ~ishandle(fig_handle)
            fig_handle = figure('Name', '径向分布', 'Position', [100, 100, 1200, 800]);
        else
            % 使用已有图形并清空
            figure(fig_handle);
            clf;
        end

        % 获取环境温度和核心半径
        k_solid = params.k_solid;
        T_amb = params.ambient_temperature;
        r_c = pState.r_c;  % 使用传入的pState参数获取核心半径
        
        % 定义温度接近阈值（当温度与环境温度差异小于1%时认为达到远场）
        temp_threshold = 0.01;  % 1%的温度差异阈值

        % 生成径向网格（使用更多点以获得更精确的远场边界）
        % 注意: 所有无量纲量都是基于金属核心半径r_c和金属核心温度T_c
        r_c_nd = 1.0;               % 无量纲核心半径
        r_p_nd = r_p / r_c;         % 无量纲颗粒表面半径
        r_f_nd = r_f / r_c;         % 无量纲火焰半径
        r_inf_nd = r_inf / r_c;     % 无量纲远场边界
        
        % 温度无量纲化
        T_p_nd = T_p / T_c;        % 无量纲氧化层表面温度
        T_amb_nd = T_amb / T_c;    % 无量纲环境温度
        
        % 使用更密集的网格来捕捉温度变化
        r_pts0_nd = linspace(r_c_nd, r_p_nd, 500);    % 区域0: 核心表面到氧化层表面
        r_pts1_nd = linspace(r_p_nd, r_f_nd, 300);    % 区域1: 氧化层表面到火焰面
        r_pts2_nd = linspace(r_f_nd, r_inf_nd, 500);  % 区域2: 火焰面到远场边界
        
        % 合并所有区域的网格，确保不重复点
        r_pts_nd = [r_pts0_nd, r_pts1_nd(2:end), r_pts2_nd(2:end)];
        
        % 物理径向网格
        r_pts = r_pts_nd * r_c;  % 转换为物理半径 (基于r_c)
        
        % 初始化存储数组
        n_pts = length(r_pts_nd);
        T_nd = zeros(n_pts, 1);     % 无量纲温度
        Y_Mg = zeros(n_pts, 1);     % Mg质量分数
        Y_CO2 = zeros(n_pts, 1);    % CO2质量分数
        Y_CO = zeros(n_pts, 1);     % CO质量分数
        Y_MgO = zeros(n_pts, 1);    % MgO质量分数
    
        mw = params.materials;
        %D_ox_gas = D_eff *10 ;
        %D_ox_gas = D_eff ;
        % 计算区域1和区域2的总质量流率比例
        m_dot_region1_total = Pe_m_region1 * 4 * pi * r_c* rho_D_gas;
        m_dot_region1_Mg = coeffs_1.Y_Mg_frac * m_dot_region1_total;
        m_dot_region1_CO = coeffs_1.Y_CO_frac * m_dot_region1_total;
        m_dot_region1_MgO = coeffs_1.Y_MgO_frac * m_dot_region1_total;
        m_dot_region1_CO2 = coeffs_1.Y_CO2_frac * m_dot_region1_total;
        


        % 为区域0估计核心表面条件
        Y_T_core_nd = coeffs_0.Y_T_core;      % 无量纲核心表面温度（假设为沸点）
        B0_nd = coeffs_0.m_dot_B0_cp_ratio;           % 热量比例参数估计值
        Y_Mg_core = coeffs_0.Y_Mg_core;        % 核心表面Mg浓度估计
        Y_Mg_frac0 = coeffs_0.Y_Mg_frac;      % Mg在区域0流率中的估计比例
        Y_CO_core =coeffs_0.Y_CO_core;        % 核心表面CO浓度估计
        Y_CO_frac0 = coeffs_0.Y_CO_frac;      % CO在区域0流率中的估计比例

        Y_Mg_surf = coeffs_1.Y_Mg_surf ;
        Y_Mg_frac1= coeffs_1.Y_Mg_frac ;
        Y_CO_surf = coeffs_1.Y_CO_surf ;
        %exp_arg_Y_region1 = Pe_m_region1 * (1/r_c_nd - 1/r_f_nd);
        %Y_Mg_region1_flame = Y_Mg_surf + (Y_Mg_surf - Y_Mg_frac1) * (exp(exp_arg_Y_region1) - 1);

        pressure = params.ambient_pressure / 101325;
        K_CO_surface = params.reaction_pre_exponential * exp(-params.reaction_activation_energy / (params.R_u * T_c));
        M_mix_surface = Y_Mg_core * mw.Mg.molar_mass + Y_CO_core * mw.CO.molar_mass;
        A_surface = 4 * pi * r_c^2 ;
        m_dot_surface_Mg = A_surface * Y_CO_core * pressure / (params.R_u * T_c) * mw.Mg.molar_mass * M_mix_surface / mw.CO.molar_mass ...
                            *  K_CO_surface;
        m_dot_surface_CO = m_dot_surface_Mg * mw.CO.molar_mass / mw.Mg.molar_mass;
        m_dot_surface_MgO = m_dot_surface_CO * mw.MgO.molar_mass / mw.CO.molar_mass;

    
        % 计算各组分物理质量流率
        m_dot_region0_Mg = m_dot_region1_Mg ;
        m_dot_region0_CO = -(abs(m_dot_region1_CO) - m_dot_surface_CO);
        m_dot_region0_total = m_dot_region0_Mg + m_dot_region0_CO;
        
        
        % 区域2各组分的质量流率 - 2区总质量流率用反应算
        m_dot_region2_CO2 =  - m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO2.molar_mass;
        m_dot_region2_Mg = 0;
        m_dot_region2_CO = m_dot_region1_Mg / mw.Mg.molar_mass * mw.CO.molar_mass - abs(m_dot_region1_CO);
        m_dot_region2_MgO = m_dot_region1_Mg / mw.Mg.molar_mass * mw.MgO.molar_mass - abs(m_dot_region1_MgO);
        m_dot_region2_total = m_dot_region2_CO2 + m_dot_region2_Mg + m_dot_region2_CO + m_dot_region2_MgO;

        % 计算佩克莱数
        % 为区域0(氧化层内)估算佩克莱数
        Pe_m_region0 = m_dot_region0_total / (4 * pi * r_c * D_ox_gas);
        Pe_m_region1 = m_dot_region1_total / (4 * pi * r_c * rho_D_gas);
        Pe_m_region2 = m_dot_region2_total / (4 * pi * r_c * rho_D_gas);
        
        % 热传导佩克莱数
        Pe_cp_region0 = m_dot_region0_total * cp_gas / (4 * pi * r_c * k_ox_gas);
        %Pe_cp_region0 = m_dot_region0_total  / (4 * pi * r_c * ratio_k_cp_ox);
        %Pe_cp_region0 = m_dot_region0_total  / (4 * pi * r_c * ratio_k_cp_ox);
        Pe_cp_region1 = m_dot_region1_total * cp_gas / (4 * pi * r_c * k_gas);
        %Pe_cp_region1 = m_dot_region1_total / (4 * pi * r_c * ratio_k_cp_gas);
        Pe_cp_region2 = m_dot_region2_total * cp_gas / (4 * pi * r_c * k_gas);
        %Pe_cp_region2 = m_dot_region2_total / (4 * pi * r_c * ratio_k_cp_gas);
        
        % 记录有效远场边界索引
        effective_boundary_idx = n_pts;
        T_amb_nd = T_amb / T_c;  % 无量纲环境温度
        
        %fprintf(' 图像中 Y_Mg_region1_flame: %.3e \n', Y_Mg_region1_flame);
        %fprintf(' 图像中 Y_Mg_flame: %.3e \n', coeffs_2.Y_Mg_flame);
        % 计算各点的温度和组分质量分数
        for i = 1:n_pts
            r_nd_i = r_pts_nd(i);
            
            % 判断区域
            if r_nd_i <= r_p_nd  % 区域0: 核心表面到氧化层表面
                % 温度分布 - 使用简化的指数分布模型
                exp_arg_T = Pe_cp_region0 * (1 - 1/r_nd_i);
                exp_term_T = exp(exp_arg_T) - 1;
                T_nd(i) = Y_T_core_nd + (Y_T_core_nd + B0_nd) * exp_term_T;
                
                % 组分分布
                exp_arg_Y = Pe_m_region0 * (1 - 1/r_nd_i);
                exp_term_Y = exp(exp_arg_Y) - 1;
                Y_Mg(i) = Y_Mg_core + (Y_Mg_core - Y_Mg_frac0) * exp_term_Y;
                Y_CO(i) = Y_CO_core + (Y_CO_core - Y_CO_frac0) * exp_term_Y;
                Y_CO2(i) = 0;  % 氧化层内无CO2
                Y_MgO(i) = 0;  % 氧化层内无气态MgO
                
            elseif r_nd_i <= r_f_nd  % 区域1: 氧化层表面到火焰面
                % 温度
                exp_arg_T = Pe_cp_region1 * (1/r_p_nd - 1/r_nd_i);
                exp_term_T = exp(exp_arg_T) - 1;
                T_nd(i) = coeffs_1.Y_T_surf + (coeffs_1.Y_T_surf + coeffs_1.m_dot_B1_cp_ratio) * exp_term_T;
                
                % 组分
                exp_arg_Y = Pe_m_region1 * (1/r_p_nd - 1/r_nd_i);
                exp_term_Y = exp(exp_arg_Y) - 1;
                Y_Mg(i) = coeffs_1.Y_Mg_surf + (coeffs_1.Y_Mg_surf - coeffs_1.Y_Mg_frac) * exp_term_Y;
                Y_CO2(i) = coeffs_1.Y_CO2_surf + (coeffs_1.Y_CO2_surf - coeffs_1.Y_CO2_frac) * exp_term_Y;
                Y_CO(i) = coeffs_1.Y_CO_surf + (coeffs_1.Y_CO_surf - coeffs_1.Y_CO_frac) * exp_term_Y;
                Y_MgO(i) = coeffs_1.Y_MgO_surf + (coeffs_1.Y_MgO_surf - coeffs_1.Y_MgO_frac) * exp_term_Y;
            else  % 区域2: 火焰面到远场
                % 温度
                exp_arg_T = Pe_cp_region2 * (1/r_f_nd - 1/r_nd_i);
                exp_term_T = exp(exp_arg_T) - 1;
                T_nd(i) = coeffs_2.Y_T_flame + (coeffs_2.Y_T_flame + coeffs_2.m_dot_B2_cp_ratio) * exp_term_T;
                
                % 组分
                exp_arg_Y = Pe_m_region2 * (1/r_f_nd - 1/r_nd_i);
                exp_term_Y = exp(exp_arg_Y) - 1;
                Y_Mg(i) = coeffs_2.Y_Mg_flame + (coeffs_2.Y_Mg_flame - coeffs_2.Y_Mg_frac) * exp_term_Y;
                Y_CO2(i) = coeffs_2.Y_CO2_flame + (coeffs_2.Y_CO2_flame - coeffs_2.Y_CO2_frac) * exp_term_Y;
                Y_CO(i) = coeffs_2.Y_CO_flame + (coeffs_2.Y_CO_flame - coeffs_2.Y_CO_frac) * exp_term_Y;
                Y_MgO(i) = coeffs_2.Y_MgO_flame + (coeffs_2.Y_MgO_flame - coeffs_2.Y_MgO_frac) * exp_term_Y;
                
                % 检查是否达到环境温度
                %if abs((T_nd(i) - T_amb_nd)/T_amb_nd) < temp_threshold
                %    effective_boundary_idx = i;
                %    fprintf('找到有效远场边界: r/r_p = %.2f (物理位置: %.2e μm)\n',...
                %            r_pts_nd(i), r_pts(i)*1e6);
                %    break;  % 找到远场边界后停止计算
                %end
            end
        end
        exp_arg_Y_nd = Pe_m_region2 * (1/r_f_nd - 1/r_inf_nd);
        Y_MgO_inf = coeffs_2.Y_MgO_flame + (coeffs_2.Y_MgO_flame - coeffs_2.Y_MgO_frac) * (exp(exp_arg_Y_nd) - 1);
        %fprintf('  Y_MgO_inf: %.3e \n', Y_MgO_inf);


        exp_arg_T = Pe_cp_region1 * (1/r_p_nd - 1/r_f_nd);
        exp_term_T = exp(exp_arg_T) - 1;
        T_flame_1 = coeffs_1.Y_T_surf + (coeffs_1.Y_T_surf + coeffs_1.m_dot_B1_cp_ratio) * exp_term_T;
        fprintf('  T_flame_1: %.3e \n', T_flame_1*T_c);
        fprintf('  T_flame: %.3e \n', coeffs_2.Y_T_flame*T_c);


        % 截断到有效远场边界
        r_pts_effective = r_pts(1:effective_boundary_idx);
        T_effective = T_nd(1:effective_boundary_idx) * T_c;  % 转换为物理温度（基于金属核心温度）
        Y_Mg_effective = Y_Mg(1:effective_boundary_idx);
        Y_CO2_effective = Y_CO2(1:effective_boundary_idx);
        Y_CO_effective = Y_CO(1:effective_boundary_idx);
        Y_MgO_effective = Y_MgO(1:effective_boundary_idx);
        
        % 温度分布
        subplot(2, 1, 1);
        % 在绘图前对数据进行样条插值平滑
        r_fine = linspace(min(r_pts_effective), max(r_pts_effective), 500);
        T_fine = interp1(r_pts_effective, T_effective, r_fine, 'spline');
        plot(r_fine*1e6, T_fine, 'k-', 'LineWidth', 2);
        hold on;
        
        % 标记区域边界
        plot([r_c, r_c]*1e6, [min(T_effective), max(T_effective)], 'b--', 'LineWidth', 1.5);
        text(r_c*1e6, max(T_effective)*0.8, '核心表面', 'Color', 'b', 'FontWeight', 'bold');
        
        plot([r_p, r_p]*1e6, [min(T_effective), max(T_effective)], 'g--', 'LineWidth', 1.5);
        text(r_p*1e6, max(T_effective)*0.85, '氧化层表面', 'Color', 'g', 'FontWeight', 'bold');
        
        plot([r_f, r_f]*1e6, [min(T_effective), max(T_effective)], 'r--', 'LineWidth', 1.5);
        text(r_f*1e6, max(T_effective)*0.9, '火焰面', 'Color', 'r', 'FontWeight', 'bold');
        
        % 标记远场边界
        r_boundary = r_pts_effective(end);
        plot([r_boundary, r_boundary]*1e6, [min(T_effective), max(T_effective)], 'm--', 'LineWidth', 1.5);
        text(r_boundary*1e6, max(T_effective)*0.75, '远场边界', 'Color', 'm', 'FontWeight', 'bold');
        
        % 添加标注点
        plot(r_c*1e6, Y_T_core_nd*T_c, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', '核心表面温度');
        plot(r_p*1e6, T_p, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', '氧化层表面温度');
        plot(r_f*1e6, coeffs_2.Y_T_flame*T_c, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '火焰温度');
        plot(r_boundary*1e6, T_amb, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', '环境温度');
        
        % 添加颗粒尺寸文本信息到图上
        text_x = min(r_pts_effective)*1e6;
        text_y = max(T_effective)*0.6;
        text_str = sprintf('镁核半径径: %.2f μm\n氧化层厚度: %.2f μm', r_c*1e6, (r_p-r_c)*1e6);
        text(text_x, text_y, text_str, 'FontSize', 10, 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'k');
        
        xlabel('径向位置 (μm)');
        ylabel('温度 (K)');
        title('温度分布 (三区域模型)');
        legend('Location', 'best');
        grid on;
        
        % 组分分布
        subplot(2, 1, 2);
        plot(r_pts_effective*1e6, Y_Mg_effective, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Mg');
        hold on;
        plot(r_pts_effective*1e6, Y_CO2_effective, 'b-', 'LineWidth', 1.5, 'DisplayName', 'CO2');
        plot(r_pts_effective*1e6, Y_CO_effective, 'g-', 'LineWidth', 1.5, 'DisplayName', 'CO');
        plot(r_pts_effective*1e6, Y_MgO_effective, 'm-', 'LineWidth', 1.5, 'DisplayName', 'MgO');
        
        % 添加区域边界线
        plot([r_c, r_c]*1e6, [0, 1.0], 'b--', 'LineWidth', 1.5);
        plot([r_p, r_p]*1e6, [0, 1.0], 'g--', 'LineWidth', 1.5);
        plot([r_f, r_f]*1e6, [0, 1.0], 'r--', 'LineWidth', 1.5);
        plot([r_boundary, r_boundary]*1e6, [0, 1.0], 'm--', 'LineWidth', 1.5);
        
        % 检查组分和是否接近1
        Y_sum_effective = Y_Mg_effective + Y_CO2_effective + Y_CO_effective + Y_MgO_effective;
        plot(r_pts_effective*1e6, Y_sum_effective, 'k:', 'LineWidth', 2, 'DisplayName', '组分和');
        
        xlabel('径向位置 (μm)');
        ylabel('质量分数');
        title('组分分布 (三区域模型)');
        legend('Location', 'best');
        ylim([0, 1.05]);
        grid on;
        
        % 添加区域标签
        region_text = {'区域0: 核心-氧化层', '区域1: 氧化层-火焰面', '区域2: 火焰面-远场'};
        region_positions = [mean([r_c, r_p])*1e6, mean([r_p, r_f])*1e6, mean([r_f, r_boundary])*1e6];
        for i = 1:3
            text(region_positions(i), 0.05, region_text{i}, 'HorizontalAlignment', 'center', 'Color', 'k', 'FontWeight', 'bold');
        end
        
        % 强制立即更新显示
        drawnow;
        
        % 保存最新图像
        saveas(fig_handle, 'radial_distributions.png');
        fprintf('径向分布图已保存为: radial_distributions.png\n');
    end
