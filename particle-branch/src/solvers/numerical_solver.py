import numpy as np

class NumericalSolver:
    def __init__(self, params, model):
        self.params = params
        self.model = model
        self.stages_history = []  # 记录每个时间步的阶段
        
        # 初始化结果存储
        self.results = {var: [] for var in self.params.record_variables}
        
    def solve(self):
        """求解颗粒燃烧过程"""
        # 初始化时间数组
        time = np.arange(0, self.params.total_time, self.params.time_step)
        n_steps = len(time)
        
        # 初始化颗粒状态
        self.model.initialize_particle(self.params.particle_diameter)
        
        # 记录初始状态
        self._record_state(0.0)
        
        # 时间步进求解
        for i in range(1, n_steps):
            current_time = time[i]
            
            # 模拟单个时间步
            state = self.model.simulate_time_step(
                dt=self.params.time_step,
                ambient_temp=self.params.ambient_temperature,
                ambient_pressure=self.params.ambient_pressure,
                gas_velocity=self.params.gas_velocity,
                ambient_composition=self.params.ambient_gas_composition
            )
            
            # 记录当前状态
            self._record_state(current_time)
            
            # 记录阶段
            self.stages_history.append(
                self.model.determine_stage(state["temperature"])
            )
            
            # 检查是否需要输出
            if i % self.params.output_interval == 0:
                self._print_progress(i, n_steps, state)
        
        return self._prepare_results()
    
    def _record_state(self, current_time):
        """记录当前状态到结果数组"""
        state = self.model.get_particle_state()
        
        # 记录基本变量
        self.results["time"].append(current_time)
        self.results["temperature"].append(state["temperature"])
        self.results["metal_radius"].append(state["metal_radius"])
        self.results["total_radius"].append(state["total_radius"])
        self.results["oxide_thickness"].append(state["oxide_thickness"])
        self.results["mass_mg"].append(state["mass_mg"])
        self.results["mass_mgo"].append(state["mass_mgo"])
        self.results["mass_c"].append(state["mass_c"])
        
        # 计算并记录能量项
        current_diameter = 2 * state["total_radius"]
        q_conv = self.model.compute_convective_heat_transfer(
            state["temperature"],
            self.params.ambient_temperature,
            current_diameter,
            self.params.gas_velocity
        )
        q_rad = self.model.compute_radiation(
            state["temperature"],
            self.params.ambient_temperature,
            current_diameter
        )
        q_reaction, reaction_rate = self.model.compute_surface_reaction(
            state["temperature"],
            self.params.ambient_temperature,
            current_diameter,
            self.params.gas_velocity,
            self.params.ambient_gas_composition,
            self.params.ambient_pressure
        )
        
        self.results["q_conv"].append(q_conv)
        self.results["q_rad"].append(q_rad)
        self.results["q_reaction"].append(q_reaction)
        self.results["reaction_rate"].append(reaction_rate)
    
    def _print_progress(self, current_step, total_steps, state):
        """打印模拟进度和关键参数"""
        progress = current_step / total_steps * 100
        print(f"\n模拟进度: {progress:.1f}%")
        print(f"当前时间: {current_step * self.params.time_step:.6f} s")
        print(f"颗粒温度: {state['temperature']:.1f} K")
        print(f"金属半径: {state['metal_radius']*1e6:.2f} μm")
        print(f"氧化层厚度: {state['oxide_thickness']*1e6:.2f} μm")
        print(f"反应阶段: {self.stages_history[-1]}")
    
    def _prepare_results(self):
        """准备返回结果"""
        # 将列表转换为numpy数组
        for key in self.results:
            self.results[key] = np.array(self.results[key])
        
        # 添加阶段历史
        self.results["stages"] = self.stages_history
        
        return self.results