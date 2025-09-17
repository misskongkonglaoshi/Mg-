import numpy as np
from .chemical_control import ChemicalControl
from .diffusion_control import DiffusionControl
from .thermo_reader import ThermoReader
import os

class PhysicalModel:
    def __init__(self, params):
        self.params = params
        # 物理常数
        self.epsilon = 0.330  # 颗粒发射率
        self.sigma = 5.67e-8  # 斯特藩-玻尔兹曼常数 (W/m²·K⁴)
        self.Ru = 8314.5  # 普适气体常数 (J/(kmol·K))
        self.gama = 1.4
        self.nu = 1e-5
        self.pr = 0.71
        self.sc = 0.6
        
        # 反应相关参数
        self.h_reaction = -601.6e3  # Mg + CO2 -> MgO + C 的反应焓 (J/mol)
        self.M_mg = 24.305  # Mg的摩尔质量 (g/mol)
        self.MgO = 40.304
        self.M_co2 = 44.01  # CO2的摩尔质量 (g/mol)
        self.A_pre = 1.376e-4  # 指前因子 (s-1)
        self.E_a = 132.4e3  # 活化能 (J/mol)
        
        # 阶段温度阈值
        self.T_melting = 923  # Mg熔点 (K)
        self.T_ignition = 1366  # Mg点火温度 (K)
        
        # 物质特性
        self.Cp_mg = 1020  # Mg比热容 (J/kg·K)
        self.h_fusion = 8.7e5  # Mg熔化潜热 (J/kg)
        self.rho_mg = 1740  # Mg密度 (kg/m³)

        # 气体特性
        self.gas_properties = {
            "O2": {"M_w": 32.00, "sigma_d": 3.467e-10, "epsilon_min": 106.7},
            "CO2": {"M_w": 44.01, "sigma_d": 3.941e-10, "epsilon_min": 195.2},
            "H2O": {"M_w": 18.01528, "sigma_d": 2.641e-10, "epsilon_min": 809.1},
            "N2": {"M_w": 28.0134, "sigma_d": 3.798E-10, "epsilon_min": 71.4},
            "MgO": {"M_w": 40.304, "sigma_d": 3.941e-10, "epsilon_min": 195.2},
            "C": {"M_w": 12.011, "sigma_d": 3.941e-10, "epsilon_min": 195.2}
        }

        # 初始化热力学数据读取器
        thermo_file = os.path.join(os.path.dirname(__file__), '..', '..', 'thermo.dat')
        self.thermo_reader = ThermoReader(thermo_file)
        
        # 初始氧化层参数
        self.initial_oxide_ratio = 0.02  # 初始氧化层占颗粒半径的比例
        
        # 颗粒状态跟踪
        self.particle_state = {
            "metal_radius": 0.0,         # 金属核心半径 (m)
            "total_radius": 0.0,         # 总半径 (包括氧化层) (m)
            "oxide_thickness": 0.0,      # 氧化层厚度 (m)
            "initial_oxide_thickness": 0.0,  # 初始氧化层厚度 (m)
            "temperature": 300.0,        # 颗粒温度 (K)
            "reaction_area_factor": 1.0,  # 有效反应面积系数
            "mass_mg": 0.0,             # 金属Mg质量 (kg)
            "mass_mgo": 0.0,            # MgO质量 (kg)
            "mass_c": 0.0,              # C质量 (kg)
            "melted_fraction": 0.0,      # 熔融金属的质量分数
            "energy_accumulated": 0.0    # 累积的能量 (J)
        }

        # 初始化控制器
        self.chemical_control = ChemicalControl()
        self.diffusion_control = DiffusionControl()

    def determine_stage(self, particle_temp):
        """确定当前的反应阶段"""
        if particle_temp < self.T_melting:
            return "preheating"
        elif particle_temp < self.T_ignition:
            return "melting"
        else:
            return "reaction"

    def calculate_energy_change(self, particle_temp, ambient_temp, particle_diameter, gas_velocity, ambient_composition):
        """计算总的能量变化"""
        stage = self.determine_stage(particle_temp)
        
        # 基础传热计算（所有阶段都有）
        q_conv = self.compute_convective_heat_transfer(particle_temp, ambient_temp, particle_diameter, gas_velocity)
        q_rad = self.compute_radiation(particle_temp, ambient_temp, particle_diameter)
        
        total_energy = q_conv + q_rad
        
        # 根据阶段增加额外的能量项
        if stage == "melting":
            q_melt = self.compute_melting(particle_temp)
            total_energy -= q_melt  # 熔化吸热
        
        elif stage == "reaction":
            q_reaction, _ = self.compute_surface_reaction(
                particle_temp, ambient_temp, particle_diameter, 
                gas_velocity, ambient_composition, self.params.ambient_pressure
            )
            total_energy += q_reaction  # 反应放热
        
        return total_energy, stage

    def calculate_temperature_change(self, particle_temp, ambient_temp, particle_diameter, 
                                  gas_velocity, ambient_composition, dt):
        """计算温度变化"""
        # 计算总的能量变化
        total_energy, stage = self.calculate_energy_change(
            particle_temp, ambient_temp, particle_diameter, gas_velocity, ambient_composition
        )
        
        # 计算质量
        particle_mass = self.calculate_particle_mass(particle_diameter)
        
        # 根据不同阶段处理温度变化
        if stage == "preheating":
            # 预热阶段，所有能量用于升温
            dT = (total_energy * dt) / (particle_mass * self.Cp_mg)
            
        elif stage == "melting":
            # 熔化阶段
            energy_for_heating, new_melted_fraction = self.compute_melting(
                particle_temp, total_energy * dt, dt
            )
            
            if energy_for_heating > 0:
                # 有剩余能量用于升温
                dT = energy_for_heating / (particle_mass * self.Cp_mg)
            else:
                # 所有能量用于熔化，温度保持不变
                dT = 0
                
        else:  # reaction stage
            # 反应阶段，所有能量用于升温
            dT = (total_energy * dt) / (particle_mass * self.Cp_mg)
            
        return dT, stage

    def calculate_particle_mass(self, particle_diameter):
        """计算颗粒质量"""
        volume = (4/3) * np.pi * (particle_diameter/2)**3
        return volume * self.rho_mg

    def compute_convective_heat_transfer(self, particle_temp, ambient_temp, particle_diameter, gas_velocity):
        # 对流换热计算，根据不同的流动regime采用不同的模型
        # 计算Knudsen数以确定流动regime
        return self.compute_convective_heat_transfer_continous(particle_temp, ambient_temp, particle_diameter,gas_velocity )
    
    def compute_convective_heat_transfer_continous(self, particle_temp, ambient_temp, particle_diameter,gas_velocity):
        # 连续流中的对流换热，使用Ranz-Marshall模型
        kg = self.calculate_thermal_conductivity(ambient_temp)  # 气体热导率 (W/m·K)，实际应根据具体条件计算
        re = self.calculate_reynolds_number(particle_diameter, gas_velocity, self.nu)
        nusselt = 2 + 0.6 * (re ** 0.5) * (self.pr ** (1/3))
        h = (kg * nusselt) / particle_diameter
        surface_area = np.pi * particle_diameter ** 2
        q_conv = h * surface_area * (ambient_temp - particle_temp)
        return q_conv
    
    def calculate_viscosity_and_thermal_conductivity(self,gas_composition,temperature):
        viscosity=0.0
        thermal_conductivity=0.0
        for species,mass_frac in gas_composition.items():
            M_w=self.gas_properties[species]["M_w"]
            sigma_d=self.gas_properties[species]["sigma_d"]
            epsilon_min=self.gas_properties[species]["epsilon_min"]

            T_star=(sigma_d*temperature)/epsilon_min
            Omega_u=self.calculate_collision_integral(T_star)
            mu=2.6693e-6*np.sqrt(M_w*temperature)/(sigma_d**2*Omega_u)
            thermal_conductivity=15/4*self.Ru/M_w*mu*(4/15*M_w*self.calculate_Cp(species,temperature)/self.Ru+1/3)
            viscosity+=mass_frac*mu
            thermal_conductivity+=mass_frac*thermal_conductivity
        return viscosity,thermal_conductivity
    
    def calculate_collision_integral(self,T_star):
        collision_integral=1.147*T_star**(-0.145)+(T_star+0.5)**(-2)
        return collision_integral
    
    def compute_radiation(self, particle_temp, ambient_temp,particle_diameter):
        # 辐射换热计算
        surface_area = np.pi * particle_diameter**2  # 假设颗粒为球形，直径为半径的两倍
        q_rad = self.epsilon * self.sigma * surface_area * (ambient_temp**4 - particle_temp**4)
        return q_rad
    
    def calculate_thermal_conductivity(self, ambient_temp):
        return 0.026 * (ambient_temp / 298.15) ** 0.75
    
    def calculate_mach_number(self, gas_velocity,gama, ambient_temp):
        return gas_velocity / np.sqrt(gama * self.R * ambient_temp)
    
    def calculate_reynolds_number(self, particle_diameter, gas_velocity, nu):
        # 计算雷诺数
        return (particle_diameter * gas_velocity) / nu
    
    def calculate_stanton_number(self):
        # 计算Stanton数，这里使用简化公式
        a = 0.9  # 热容系数
        c_star = 1.4  # 平均比热比
        s = 1.0  # 分子速度比
        erf_s = 0.8427  # error function of s
        erfc_s = 1 - erf_s
        ierfc_s = 0.5 * np.exp(-s**2)  # integral of complementary error function
        st = (1/8 * a * c_star + 1/c_star * (ierfc_s + 0.5 * s * erfc_s)) / (s * erfc_s)
        return st
    
    def get_coeffs(self, species, T):
        """获取对应温度下的NASA系数"""
        return self.thermo_reader.get_coeffs(species, T)

    def calculate_Cp(self, species, T):
        """使用NASA七系数多项式计算比热容 (J/mol·K)"""
        return self.thermo_reader.calculate_Cp(species, T, self.Ru)

    def calculate_enthalpy(self, species, T):
        """使用NASA七系数多项式计算标准摩尔焓 (J/mol)"""
        return self.thermo_reader.calculate_H(species, T, self.Ru)

    def calculate_entropy(self, species, T):
        """使用NASA七系数多项式计算标准摩尔熵 (J/mol·K)"""
        return self.thermo_reader.calculate_S(species, T, self.Ru)

    def calculate_gibbs(self, species, T):
        """计算标准摩尔吉布斯自由能 (J/mol)
        
        G° = H° - TS°
        """
        H = self.calculate_enthalpy(species, T)
        S = self.calculate_entropy(species, T)
        return H - T * S

    def calculate_reaction_enthalpy(self, T):
        """计算反应焓变 (J/mol)"""
        # Mg + CO2 -> MgO + C
        # 计算生成物焓
        H_products = self.calculate_enthalpy("MgO", T) + self.calculate_enthalpy("C", T)
        
        # 计算反应物焓
        H_reactants = self.calculate_enthalpy("Mg", T) + self.calculate_enthalpy("CO2", T)
        
        # 计算反应焓变
        return H_products - H_reactants

    def compute_surface_reaction(self, particle_temp, ambient_temp, particle_diameter, 
                               gas_velocity, ambient_composition, ambient_pressure):
        """计算表面反应速率和热量"""
        # 获取当前的反应面积系数
        reaction_area_factor = self.particle_state["reaction_area_factor"]
        
        # 计算扩散控制下的反应速率
        r_diff = self.diffusion_control.compute_diffusion_controlled_rate(
            particle_temp, ambient_temp, particle_diameter,
            gas_velocity, ambient_composition, ambient_pressure,
            reaction_area_factor, self.nu
        )
        
        # 计算化学动力学控制下的反应速率
        r_chem = self.chemical_control.compute_chemical_controlled_rate(
            particle_temp, particle_diameter, ambient_composition,
            ambient_pressure, reaction_area_factor
        )


        # 取较小值作为实际反应速率
        reaction_rate = min(r_diff, r_chem)
        
        # 计算反应热
        
        reaction_heat = reaction_rate * self.chemical_control.calculate_reaction_enthalpy(particle_temp)

        return reaction_heat, reaction_rate
    
    def compute_melting(self, particle_temp, total_energy, dt):
        """
        计算熔化过程
        
        Args:
            particle_temp: 当前温度 (K)
            total_energy: 当前时间步的总能量输入 (J)
            dt: 时间步长 (s)
            
        Returns:
            tuple: (实际用于升温的能量, 更新后的熔融分数)
        """
        if particle_temp < self.T_melting:
            return total_energy, self.particle_state["melted_fraction"]
            
        # 如果温度达到熔点但未完全熔化
        if self.particle_state["melted_fraction"] < 1.0:
            # 计算完全熔化所需的能量
            metal_mass = self.particle_state["mass_mg"]
            total_fusion_energy = metal_mass * self.h_fusion
            
            # 计算当前未熔化部分所需的熔化能
            remaining_fusion_energy = total_fusion_energy * (1 - self.particle_state["melted_fraction"])
            
            # 累积能量
            self.particle_state["energy_accumulated"] += total_energy
            
            if self.particle_state["energy_accumulated"] >= remaining_fusion_energy:
                # 能量足够完成熔化
                self.particle_state["melted_fraction"] = 1.0
                # 返回剩余能量用于升温
                excess_energy = self.particle_state["energy_accumulated"] - remaining_fusion_energy
                self.particle_state["energy_accumulated"] = 0
                return excess_energy, 1.0
            else:
                # 能量不足以完成熔化，全部用于熔化过程
                new_melted_fraction = self.particle_state["melted_fraction"] + \
                                    (self.particle_state["energy_accumulated"] / total_fusion_energy)
                self.particle_state["melted_fraction"] = new_melted_fraction
                self.particle_state["energy_accumulated"] = 0
                return 0, new_melted_fraction
                
        # 已完全熔化，所有能量用于升温
        return total_energy, 1.0

    def calculate_dTdt(self, particle_temp, ambient_temp, ambient_pressure, particle_diameter, gas_velocity, ambient_composition):
        # 能量方程
        total_energy, stage = self.calculate_energy_change(
            particle_temp, ambient_temp, particle_diameter, gas_velocity, ambient_composition
        )
        particle_mass = self.calculate_particle_mass(particle_diameter)
        dTdt = total_energy / (particle_mass * self.Cp_mg)
        return dTdt, stage

    def initialize_particle(self, initial_diameter):
        """初始化颗粒状态"""
        # 计算初始氧化层厚度
        initial_oxide_thickness = initial_diameter * self.initial_oxide_ratio
        
        # 计算金属核心半径
        metal_radius = (initial_diameter - 2 * initial_oxide_thickness) / 2
        
        # 更新颗粒状态
        self.particle_state.update({
            "metal_radius": metal_radius,
            "total_radius": initial_diameter / 2,
            "oxide_thickness": initial_oxide_thickness,
            "initial_oxide_thickness": initial_oxide_thickness,
            "temperature": 300.0,  # 初始温度设为300K
            "reaction_area_factor": 1.0,
            "mass_mg": self.calculate_particle_mass(2 * metal_radius),  # 金属核心质量
            "mass_mgo": 0.0,  # 初始MgO质量
            "mass_c": 0.0,  # 初始C质量
            "melted_fraction": 0.0,  # 初始熔融分数为0
            "energy_accumulated": 0.0  # 初始累积能量为0
        })
        
        # 更新反应面积系数
        self.update_reaction_area_factor()

    def update_reaction_area_factor(self):
        """根据氧化层厚度更新有效反应面积系数
        
        随着氧化层厚度的增加，有效反应面积会减小
        """
        # 基于氧化层厚度计算有效反应面积系数
        oxide_thickness = self.particle_state["oxide_thickness"]
        initial_thickness = self.particle_state["initial_oxide_thickness"]
        
        if oxide_thickness <= initial_thickness:
            self.particle_state["reaction_area_factor"] = 1.0
        else:
            # 随着氧化层厚度增加，有效反应面积系数减小
            self.particle_state["reaction_area_factor"] = initial_thickness / oxide_thickness
    
    def calculate_effective_area(self):
        """计算有效反应面积
        
        Returns:
            float: 有效反应面积 (m²)
        """
        # 计算总表面积
        total_radius = self.particle_state["total_radius"]
        surface_area = np.pi * (2 * total_radius) ** 2
        
        # 应用有效反应面积系数
        return surface_area * self.particle_state["reaction_area_factor"]

    def update_particle_size(self, reacted_mg_mass):
        """更新颗粒尺寸"""
        if reacted_mg_mass <= 0:
            return
            
        # 更新质量
        self.particle_state["mass_mg"] -= reacted_mg_mass
        
        # 计算生成的MgO和C的质量
        produced_mgo_mass = reacted_mg_mass * (self.chemical_control.MgO / self.chemical_control.M_mg)
        produced_c_mass = reacted_mg_mass * (12.011 / self.chemical_control.M_mg)  # 12.011是C的摩尔质量
        
        self.particle_state["mass_mgo"] += produced_mgo_mass
        self.particle_state["mass_c"] += produced_c_mass
        
        # 更新金属核心半径
        if self.particle_state["mass_mg"] > 0:
            new_metal_volume = self.particle_state["mass_mg"] / self.rho_mg
            self.particle_state["metal_radius"] = (3 * new_metal_volume / (4 * np.pi))**(1/3)
        else:
            self.particle_state["metal_radius"] = 0
            
        # 计算氧化层体积
        oxide_volume = self.particle_state["mass_mgo"] / 3580  # 3580 kg/m³是MgO的密度
        carbon_volume = self.particle_state["mass_c"] / 2260   # 2260 kg/m³是C的密度
        total_product_volume = oxide_volume + carbon_volume
        
        # 计算新的总半径（考虑金属核心和产物层）
        if total_product_volume > 0:
            total_volume = (4/3) * np.pi * self.particle_state["metal_radius"]**3 + total_product_volume
            self.particle_state["total_radius"] = (3 * total_volume / (4 * np.pi))**(1/3)
        else:
            self.particle_state["total_radius"] = self.particle_state["metal_radius"]
            
        # 更新氧化层厚度
        self.particle_state["oxide_thickness"] = (
            self.particle_state["total_radius"] - self.particle_state["metal_radius"]
        )
        
        # 更新反应面积系数
        self.update_reaction_area_factor()

    def simulate_time_step(self, dt, ambient_temp, ambient_pressure, gas_velocity, ambient_composition):
        """模拟一个时间步长"""
        # 获取当前状态
        particle_temp = self.particle_state["temperature"]
        particle_diameter = 2 * self.particle_state["total_radius"]
        
        # 计算温度变化
        dT, stage = self.calculate_temperature_change(
            particle_temp, ambient_temp, particle_diameter,
            gas_velocity, ambient_composition, dt
        )
        
        # 更新温度
        new_temp = particle_temp + dT
        self.particle_state["temperature"] = new_temp
        
        # 如果处于反应阶段，计算反应进度
        if stage == "reaction":
            _, reaction_rate = self.compute_surface_reaction(
                new_temp, ambient_temp, particle_diameter,
                gas_velocity, ambient_composition, ambient_pressure
            )
            
            # 计算反应的Mg质量
            reacted_mg_mass = reaction_rate * dt
            
            # 更新颗粒尺寸和组分
            self.update_particle_size(reacted_mg_mass)
        
        return stage

    def get_particle_state(self):
        """获取当前颗粒状态"""
        return self.particle_state.copy()