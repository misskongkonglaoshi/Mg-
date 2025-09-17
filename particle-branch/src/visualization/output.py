import matplotlib.pyplot as plt
import numpy as np

class Output:
    @staticmethod
    def plot_results(results):
        """绘制模拟结果
        
        Args:
            results: 包含time, temperature, composition, energy_terms和stages的字典
        """
        time = results["time"]
        temperature = results["temperature"]
        composition = results["composition"]
        energy_terms = results["energy_terms"]
        stages = results["stages"]
        
        # 创建一个2x2的子图布局
        fig = plt.figure(figsize=(15, 12))
        
        # 1. 温度随时间的变化
        ax1 = plt.subplot(2, 2, 1)
        ax1.plot(time * 1000, temperature, 'r-', label='Temperature')
        ax1.set_xlabel('Time (ms)')
        ax1.set_ylabel('Temperature (K)')
        ax1.set_title('Temperature Evolution')
        ax1.grid(True)
        
        # 添加阶段分界线和标注
        stage_changes = Output._find_stage_changes(stages)
        for stage_time, stage_name in stage_changes:
            if stage_time is not None:
                ax1.axvline(x=stage_time * 1000, color='k', linestyle='--', alpha=0.5)
                ax1.text(stage_time * 1000, ax1.get_ylim()[1], f'{stage_name}',
                        rotation=90, verticalalignment='top')
        
        # 2. 组分随时间的变化
        ax2 = plt.subplot(2, 2, 2)
        for species, values in composition.items():
            ax2.plot(time * 1000, values, label=species)
        ax2.set_xlabel('Time (ms)')
        ax2.set_ylabel('Mole Fraction')
        ax2.set_title('Species Evolution')
        ax2.legend()
        ax2.grid(True)
        
        # 3. 能量项随时间的变化
        ax3 = plt.subplot(2, 2, 3)
        for term, values in energy_terms.items():
            if term != 'total':  # 总能量单独绘制
                ax3.plot(time * 1000, values, label=term)
        ax3.set_xlabel('Time (ms)')
        ax3.set_ylabel('Energy (W)')
        ax3.set_title('Energy Terms')
        ax3.legend()
        ax3.grid(True)
        
        # 4. 总能量随时间的变化
        ax4 = plt.subplot(2, 2, 4)
        ax4.plot(time * 1000, energy_terms['total'], 'k-', label='Total Energy')
        ax4.set_xlabel('Time (ms)')
        ax4.set_ylabel('Energy (W)')
        ax4.set_title('Total Energy Evolution')
        ax4.grid(True)
        
        plt.tight_layout()
        plt.show()
    
    @staticmethod
    def _find_stage_changes(stages):
        """找出阶段变化的时间点"""
        stage_changes = []
        current_stage = None
        
        for i, stage in enumerate(stages):
            if stage != current_stage:
                stage_changes.append((i, stage))
                current_stage = stage
        
        return stage_changes
    
    @staticmethod
    def save_results(results, filename):
        """保存模拟结果到文件"""
        np.savez(filename,
                 time=results["time"],
                 temperature=results["temperature"],
                 composition=results["composition"],
                 energy_terms=results["energy_terms"],
                 stages=results["stages"])