# 附录A：ICI计算栈：完整代码与SOP协议



> Converted from the supplied Word appendix. Review formulas and tables before publication.





## A.1 概述与使用说明

本附录提供ICI框架的完整计算实现，包括：

标准操作程序（SOP）：参数测量和估算的操作步骤

完整Python代码：ICI和R(t)的计算、不确定性传播、早期预警信号检测

跨域映射验证工具：KSG互信息检验和Procrustes对齐检验

可视化工具：R(t)时间轨迹图、ICI层级图、预警仪表盘

依赖库：

# 核心计算

numpy >= 1.24.0

scipy >= 1.10.0

pandas >= 2.0.0



# 信息论分析

scikit-learn >= 1.2.0

npeet >= 0.1.0          # KSG互信息估计



# 统计与贝叶斯

pymc >= 5.0.0           # 贝叶斯标定

arviz >= 0.15.0         # 贝叶斯诊断



# 可视化

matplotlib >= 3.7.0

seaborn >= 0.12.0

plotly >= 5.14.0        # 交互式仪表盘



# 时间序列

statsmodels >= 0.14.0

hurst >= 0.0.5          # Hurst指数计算



## A.2 标准操作程序（SOP）

### SOP-01：参数提取通用协议

文件编号：SOP-01

版本：3.0

适用域：所有ICI应用域

最后更新：2025年



目的：规范六参数的提取流程，确保跨域可比性



步骤1：域识别与映射规则选择

─────────────────────────────

1.1 确认目标系统所属域：

□ 生物域（细胞/神经系统/多细胞生物）

□ 历史域（帝国/大型政治组织）

□ 经济域（金融系统/宏观经济）

□ 生态域（生态系统/种群）

□ AI域（人工神经网络/智能体）



1.2 选择对应的域特定映射规则（附录D）



1.3 记录映射规则的版本号和最后验证日期



步骤2：数据质量评估

─────────────────────────────

2.1 对每个参数，确认数据来源层级：

□ 第一层（直接测量，误差±20-30%）

□ 第二层（间接推算，误差±40-60%）

□ 第三层（专家判断，误差±60-90%）



2.2 记录每个参数的：

- 数据来源（文献引用/数据库条目）

- 估算方法

- 点估计值

- 不确定性范围（68%置信区间）



步骤3：参数计算

─────────────────────────────

3.1 计算DCS乘积：

DCS = D × C × S

lg(DCS) = lg(D) + lg(C) + lg(S)



3.2 计算FWM乘积：

FWM = F × W × M

使用FWM_h = 7.52 × 10^11（贝叶斯标定值）



3.3 计算ICI：

ICI = k × lg(DCS) × (1 + sqrt(α × FWM/FWM_h))

其中 k = 1.259，α = 1.02 × 10^5



3.4 计算R(t)（域内相对值）：

R_abs = lg(FWM) - lg(DCS)

R_relative = R_abs - R_reference（以选定基准点为零）



步骤4：不确定性传播

─────────────────────────────

4.1 使用蒙特卡洛方法：

- 对每个参数的对数值从正态分布采样

- 采样次数：N = 10000

- 计算ICI和R(t)的分布



4.2 报告：

- 点估计（中位数）

- 68%置信区间（16th至84th百分位）

- 95%置信区间（2.5th至97.5th百分位）



步骤5：结果验证

─────────────────────────────

5.1 检查ICI值是否在域的合理范围内

5.2 检查R(t)的绝对值是否与已知系统的

健康/危机状态一致

5.3 与已有文献中的ICI估算进行比较

（如有）



步骤6：文档记录

─────────────────────────────

6.1 完整填写参数记录表（见SOP-01附表）

6.2 保存原始数据和计算脚本

6.3 生成标准报告（使用附录A的报告模板）

### SOP-02：R(t)时间序列监测协议

文件编号：SOP-02

版本：2.1

适用域：生态、经济、历史（时间序列数据充足时）



步骤1：时间序列构建

─────────────────────────────

1.1 确定时间步长：

- 生态域：月至年

- 经济域：日至月

- 历史域：年至十年



1.2 对每个时间点估算六参数

（使用SOP-01）



1.3 计算R(t)时间序列



步骤2：趋势分析

─────────────────────────────

2.1 计算R(t)的滚动均值（窗口=10个时间步）

2.2 使用Mann-Kendall检验趋势显著性

- tau < -0.3且p < 0.05：显著下降趋势

- tau > 0.3且p < 0.05：显著上升趋势



步骤3：临界慢化信号检测

─────────────────────────────

3.1 滚动方差（窗口=20个时间步）：

signal_1 = rolling_variance(R_t, window=20)

trend_1 = kendall_tau(signal_1)



3.2 滚动AR(1)自相关（窗口=20个时间步）：

signal_2 = rolling_AR1(R_t, window=20)

trend_2 = kendall_tau(signal_2)



3.3 滚动Hurst指数（窗口=40个时间步）：

signal_3 = rolling_hurst(R_t, window=40)

transition_3 = hurst_transition(signal_3)



步骤4：预警级别判定

─────────────────────────────

计算综合预警分数：

score = 0

if R_t_rolling_mean < -8.5: score += 2

if R_t_rolling_mean < -9.5: score += 2

if trend_1 > 0.4: score += 1    # 方差膨胀

if trend_2 > 0.4: score += 1    # 自相关拖尾

if transition_3: score += 2      # Hurst跃迁



预警级别：

score 0-1: 绿色（正常）

score 2-3: 黄色（关注）

score 4-5: 橙色（预警）

score 6+:  红色（危机）



## A.3 核心计算库

### A.3.1 主计算模块

"""

ici_core.py

ICI框架核心计算模块

版本：3.0.0

"""



import numpy as np

from scipy import stats

from dataclasses import dataclass, field

from typing import Optional, Tuple, Dict

import warnings





# ============================================================

# 全局常数

# ============================================================



ICI_CONSTANTS = {

'k': 1.259,                    # 比例常数（贝叶斯标定值）

'alpha': 1.02e5,               # 涌现项系数

'FWM_h': 7.52e11,             # 海马CA1归一化基准

'R_critical': -10.0,           # 临界R(t)值

'R_resilient': -9.0,           # 强韧稳态上界

'R_collapse': -11.0,           # 不可逆崩溃下界

}



# 各域的参数误差估计（对数标准差）

DOMAIN_UNCERTAINTY = {

'biological_lab':    {'D':0.02, 'C':0.04, 'S':0.13,

'F':0.09, 'W':0.11, 'M':0.17},

'biological_field':  {'D':0.06, 'C':0.11, 'S':0.22,

'F':0.26, 'W':0.17, 'M':0.35},

'ecological':        {'D':0.06, 'C':0.11, 'S':0.22,

'F':0.26, 'W':0.17, 'M':0.35},

'economic_modern':   {'D':0.04, 'C':0.06, 'S':0.17,

'F':0.22, 'W':0.13, 'M':0.26},

'historical_empire': {'D':0.17, 'C':0.22, 'S':0.30,

'F':0.35, 'W':0.30, 'M':0.39},

'ai_architecture':   {'D':0.04, 'C':0.02, 'S':0.04,

'F':0.22, 'W':0.26, 'M':0.30},

}





@dataclass

class ICIParameters:

"""ICI六参数容器，包含不确定性"""

D: float        # 分子/功能类型数

C: float        # 总量

S: float        # 加权吞吐量 (s^-1)

F: float        # 反馈回路数

W: float        # 归一化交互频率

M: float        # 加权记忆复杂度



# 不确定性（以对数标准差表示）

D_sigma: float = 0.10

C_sigma: float = 0.10

S_sigma: float = 0.15

F_sigma: float = 0.20

W_sigma: float = 0.20

M_sigma: float = 0.25



# 元数据

system_name: str = "unnamed"

domain: str = "unknown"

time_point: Optional[str] = None

data_sources: Dict = field(default_factory=dict)



def __post_init__(self):

"""验证参数合理性"""

for name, val in [('D', self.D), ('C', self.C), ('S', self.S),

('F', self.F), ('W', self.W), ('M', self.M)]:

if val <= 0:

raise ValueError(f"参数{name}必须为正值，当前值：{val}")

if val < 1e-50:

warnings.warn(f"参数{name}极小（{val:.2e}），"

f"请确认是否已正确归一化")



@classmethod

def from_domain_defaults(cls, domain: str, **kwargs) -> 'ICIParameters':

"""从域默认误差创建参数对象"""

uncertainty = DOMAIN_UNCERTAINTY.get(domain, DOMAIN_UNCERTAINTY['biological_field'])

defaults = {

f'{k}_sigma': v

for k, v in uncertainty.items()

}

defaults.update(kwargs)

defaults['domain'] = domain

return cls(**defaults)





@dataclass

class ICIResult:

"""ICI计算结果"""

ICI: float

ICI_lower: float    # 16th百分位

ICI_upper: float    # 84th百分位

ICI_lower_95: float # 2.5th百分位

ICI_upper_95: float # 97.5th百分位



R_t: float

R_t_lower: float

R_t_upper: float



baseline: float     # lg(DCS)

emergence: float    # sqrt(alpha * FWM/FWM_h)



alert_level: str    # 'green', 'yellow', 'orange', 'red'

params: ICIParameters



def summary(self) -> str:

"""生成简洁摘要"""

lines = [

f"系统：{self.params.system_name}",

f"域：{self.params.domain}",

f"",

f"ICI = {self.ICI:.1f} [{self.ICI_lower:.1f}, {self.ICI_upper:.1f}] (68%CI)",

f"     [{self.ICI_lower_95:.1f}, {self.ICI_upper_95:.1f}] (95%CI)",

f"",

f"R(t) = {self.R_t:.2f} [{self.R_t_lower:.2f}, {self.R_t_upper:.2f}] (68%CI)",

f"",

f"基线项 lg(DCS) = {self.baseline:.2f}",

f"涌现项 sqrt(α·FWM/FWM_h) = {self.emergence:.4f}",

f"",

f"预警级别：{self.alert_level.upper()}",

]

return "\n".join(lines)





def compute_ici(params: ICIParameters,

n_samples: int = 10000,

reference_Rt: Optional[float] = None) -> ICIResult:

"""

计算ICI值及其不确定性



参数

────

params : ICIParameters

六参数对象，包含点估计和不确定性

n_samples : int

蒙特卡洛采样数（默认10000）

reference_Rt : float, optional

R(t)的参考值（用于域内相对R(t)计算）

若为None，返回绝对R(t)



返回

────

ICIResult : 完整计算结果

"""

k     = ICI_CONSTANTS['k']

alpha = ICI_CONSTANTS['alpha']

FWM_h = ICI_CONSTANTS['FWM_h']



# 点估计计算

# ──────────

DCS      = params.D * params.C * params.S

FWM      = params.F * params.W * params.M

baseline = np.log10(DCS)

fwm_ratio = FWM / FWM_h

emergence = np.sqrt(alpha * fwm_ratio) if fwm_ratio > 0 else 0.0

ICI_point = k * baseline * (1 + emergence)

R_t_abs   = np.log10(FWM) - np.log10(DCS)

R_t_point = R_t_abs if reference_Rt is None else R_t_abs - reference_Rt



# 蒙特卡洛不确定性传播

# ────────────────────

rng = np.random.default_rng(42)



# 对各参数的对数值采样（正态分布）

log_D = rng.normal(np.log10(params.D), params.D_sigma, n_samples)

log_C = rng.normal(np.log10(params.C), params.C_sigma, n_samples)

log_S = rng.normal(np.log10(params.S), params.S_sigma, n_samples)

log_F = rng.normal(np.log10(params.F), params.F_sigma, n_samples)

log_W = rng.normal(np.log10(params.W), params.W_sigma, n_samples)

log_M = rng.normal(np.log10(params.M), params.M_sigma, n_samples)



# 转回线性空间

D_s = 10**log_D

C_s = 10**log_C

S_s = 10**log_S

F_s = 10**log_F

W_s = 10**log_W

M_s = 10**log_M



# 计算样本ICI

DCS_s      = D_s * C_s * S_s

FWM_s      = F_s * W_s * M_s

baseline_s = np.log10(DCS_s)

ratio_s    = np.clip(FWM_s / FWM_h, 0, None)

emerg_s    = np.sqrt(alpha * ratio_s)

ICI_s      = k * baseline_s * (1 + emerg_s)



# 计算样本R(t)

log_FWM_s = np.log10(np.clip(FWM_s, 1e-300, None))

log_DCS_s = np.log10(np.clip(DCS_s, 1e-300, None))

Rt_s      = log_FWM_s - log_DCS_s

if reference_Rt is not None:

Rt_s = Rt_s - reference_Rt



# 计算分位数

pcts = np.percentile(ICI_s, [2.5, 16, 50, 84, 97.5])

Rt_pcts = np.percentile(Rt_s, [16, 50, 84])



# 确定预警级别（使用R(t)的中位数）

Rt_median = Rt_pcts[1]

alert = _determine_alert_level(Rt_median)



return ICIResult(

ICI         = ICI_point,

ICI_lower   = pcts[1],

ICI_upper   = pcts[3],

ICI_lower_95 = pcts[0],

ICI_upper_95 = pcts[4],

R_t         = R_t_point,

R_t_lower   = Rt_pcts[0],

R_t_upper   = Rt_pcts[2],

baseline    = baseline,

emergence   = emergence,

alert_level = alert,

params      = params,

)





def _determine_alert_level(Rt: float) -> str:

"""根据R(t)值确定预警级别"""

if Rt > ICI_CONSTANTS['R_resilient']:

return 'green'

elif Rt > -9.5:

return 'yellow'

elif Rt > ICI_CONSTANTS['R_collapse']:

return 'orange'

else:

return 'red'





def compute_ici_timeseries(params_list: list,

reference_index: int = 0) -> list:

"""

计算ICI时间序列



参数

────

params_list : list of ICIParameters

各时间点的参数列表（时序排列）

reference_index : int

用作R(t)基准的时间点索引



返回

────

list of ICIResult

"""

# 计算基准R(t)

ref_result = compute_ici(params_list[reference_index])

ref_Rt     = ref_result.R_t



return [

compute_ici(p, reference_Rt=ref_Rt)

for p in params_list

]

### A.3.2 临界慢化检测模块

"""

ici_csd.py

临界慢化（Critical Slowing Down）检测模块

"""



import numpy as np

import pandas as pd

from scipy import stats

from typing import Tuple, Dict, Optional





def rolling_variance(series: np.ndarray,

window: int = 20) -> np.ndarray:

"""

计算滚动方差



参数

────

series : 时间序列（等间隔）

window : 滚动窗口大小



返回

────

滚动方差数组（前window-1个元素为NaN）

"""

s = pd.Series(series)

return s.rolling(window=window, min_periods=window).var().values





def rolling_ar1(series: np.ndarray,

window: int = 20) -> np.ndarray:

"""

计算滚动一阶自回归系数 AR(1)



使用OLS估计：x(t) = ρ·x(t-1) + ε

"""

n = len(series)

ar1 = np.full(n, np.nan)



for i in range(window - 1, n):

sub = series[i - window + 1: i + 1]

x   = sub[:-1] - sub[:-1].mean()

y   = sub[1:]  - sub[1:].mean()

if np.std(x) > 1e-10:

ar1[i] = np.dot(x, y) / np.dot(x, x)



return ar1





def rolling_hurst(series: np.ndarray,

window: int = 40,

min_window: int = 8) -> np.ndarray:

"""

计算滚动Hurst指数（R/S分析）



H = 0.5：随机游走

H > 0.5：持久性（趋势惯性）

H < 0.5：反持久性（均值回归）

"""

n    = len(series)

hurst_vals = np.full(n, np.nan)



for i in range(window - 1, n):

sub = series[i - window + 1: i + 1]

h   = _hurst_rs(sub)

if h is not None:

hurst_vals[i] = h



return hurst_vals





def _hurst_rs(series: np.ndarray) -> Optional[float]:

"""R/S分析计算单个Hurst指数"""

n = len(series)

if n < 8:

return None



lags   = range(4, n // 2 + 1)

rs_vals = []



for lag in lags:

chunks = [series[i:i+lag] for i in range(0, n - lag + 1, lag)]

rs_chunk = []

for chunk in chunks:

mean_c = np.mean(chunk)

devs   = np.cumsum(chunk - mean_c)

R      = np.max(devs) - np.min(devs)

S      = np.std(chunk, ddof=1)

if S > 1e-10:

rs_chunk.append(R / S)

if rs_chunk:

rs_vals.append((lag, np.mean(rs_chunk)))



if len(rs_vals) < 3:

return None



lags_log = np.log([r[0] for r in rs_vals])

rs_log   = np.log([r[1] for r in rs_vals])



slope, _, r, _, _ = stats.linregress(lags_log, rs_log)

return slope if abs(r) > 0.7 else None





def kendall_tau_trend(series: np.ndarray) -> Tuple[float, float]:

"""

计算Kendall's τ趋势统计量



返回

────

(tau, p_value)

"""

x   = np.arange(len(series))

tau, p = stats.kendalltau(x, series)

return tau, p





def detect_hurst_transition(hurst_series: np.ndarray,

threshold: float = 0.5,

min_delta: float = 0.2) -> bool:

"""

检测Hurst指数是否发生了从<0.5到>0.5的跃迁



参数

────

hurst_series : 滚动Hurst指数序列

threshold    : 跃迁阈值（默认0.5）

min_delta    : 最小跃迁幅度（默认0.2）



返回

────

bool：是否检测到显著跃迁

"""

valid = hurst_series[~np.isnan(hurst_series)]

if len(valid) < 10:

return False



mid   = len(valid) // 2

first = valid[:mid]

last  = valid[mid:]



mean_first = np.mean(first)

mean_last  = np.mean(last)



transition = (mean_first < threshold and

mean_last  > threshold and

(mean_last - mean_first) > min_delta)

return transition





def compute_csd_signals(Rt_series: np.ndarray,

window_var: int = 20,

window_ar1: int = 20,

window_hurst: int = 40) -> Dict:

"""

计算完整的临界慢化信号套件



返回

────

dict containing:

'variance'     : 滚动方差序列

'ar1'          : 滚动AR(1)序列

'hurst'        : 滚动Hurst指数序列

'tau_var'      : 方差趋势的Kendall τ

'p_var'        : 方差趋势的p值

'tau_ar1'      : AR(1)趋势的Kendall τ

'p_ar1'        : AR(1)趋势的p值

'hurst_trans'  : Hurst指数是否发生跃迁

'alert_score'  : 综合预警分数 (0-6)

'alert_level'  : 'green'/'yellow'/'orange'/'red'

"""

var_s   = rolling_variance(Rt_series, window=window_var)

ar1_s   = rolling_ar1(Rt_series, window=window_ar1)

hurst_s = rolling_hurst(Rt_series, window=window_hurst)



# 只对非NaN部分计算趋势

var_valid   = var_s[~np.isnan(var_s)]

ar1_valid   = ar1_s[~np.isnan(ar1_s)]

hurst_valid = hurst_s[~np.isnan(hurst_s)]



tau_var, p_var   = (kendall_tau_trend(var_valid)

if len(var_valid) >= 5 else (0, 1))

tau_ar1, p_ar1   = (kendall_tau_trend(ar1_valid)

if len(ar1_valid) >= 5 else (0, 1))

hurst_trans      = detect_hurst_transition(hurst_s)



# 综合预警分数

score = 0

Rt_recent = np.nanmean(Rt_series[-5:]) if len(Rt_series) >= 5 else Rt_series[-1]



if Rt_recent < -8.5:  score += 2

if Rt_recent < -9.5:  score += 2

if tau_var > 0.4 and p_var < 0.05:  score += 1

if tau_ar1 > 0.4 and p_ar1 < 0.05: score += 1

if hurst_trans:                      score += 2



alert_map = {0:'green', 1:'green', 2:'yellow', 3:'yellow',

4:'orange', 5:'orange', 6:'red', 7:'red', 8:'red'}

alert     = alert_map.get(min(score, 8), 'red')



return {

'variance'   : var_s,

'ar1'        : ar1_s,

'hurst'      : hurst_s,

'tau_var'    : tau_var,

'p_var'      : p_var,

'tau_ar1'    : tau_ar1,

'p_ar1'      : p_ar1,

'hurst_trans': hurst_trans,

'alert_score': score,

'alert_level': alert,

}

### A.3.3 跨域验证模块

"""

ici_crossdomain.py

跨域映射验证：KSG互信息检验 + Procrustes对齐

"""



import numpy as np

from scipy.spatial import procrustes

from scipy.stats import permutation_test

from typing import Tuple, Optional

import warnings



try:

from npeet.entropy_estimators import mi as ksg_mi

HAS_NPEET = True

except ImportError:

HAS_NPEET = False

warnings.warn("npeet未安装，KSG互信息将使用近似方法")





def ksg_mutual_information(X: np.ndarray,

Y: np.ndarray,

k: int = 5) -> float:

"""

KSG互信息估计（Kraskov-Stögbauer-Grassberger）



参数

────

X, Y : 形状为(n_samples, n_features)的数组

k    : 近邻数（默认5）



返回

────

互信息估计值（bits）

"""

if X.ndim == 1:

X = X.reshape(-1, 1)

if Y.ndim == 1:

Y = Y.reshape(-1, 1)



if HAS_NPEET:

# 使用精确的KSG估计

return ksg_mi(X.tolist(), Y.tolist(), k=k)

else:

# 近似方法（基于直方图，精度较低）

return _approx_mi(X, Y)





def _approx_mi(X: np.ndarray, Y: np.ndarray,

bins: int = 10) -> float:

"""互信息的直方图近似（备用方法）"""

if X.ndim > 1:

X = X[:, 0]

if Y.ndim > 1:

Y = Y[:, 0]



hist_XY, _, _ = np.histogram2d(X, Y, bins=bins)

hist_X  = hist_XY.sum(axis=1)

hist_Y  = hist_XY.sum(axis=0)



pXY = hist_XY / hist_XY.sum()

pX  = hist_X  / hist_X.sum()

pY  = hist_Y  / hist_Y.sum()



mi = 0.0

for i in range(bins):

for j in range(bins):

if pXY[i,j] > 1e-10 and pX[i] > 1e-10 and pY[j] > 1e-10:

mi += pXY[i,j] * np.log2(pXY[i,j] / (pX[i] * pY[j]))

return mi





def functional_equivalence_test(source_params: np.ndarray,

target_params: np.ndarray,

param_names: list,

threshold_bits: float = 0.5,

n_permutations: int = 1000

) -> dict:

"""

功能等价性检验（KSG互信息检验）



参数

────

source_params : 形状(n_timepoints, 6)的源域参数矩阵

target_params : 形状(n_timepoints, 6)的目标域参数矩阵

param_names   : 六个参数的名称列表

threshold_bits: 互信息阈值（默认0.5 bits）

n_permutations: 置换检验次数



返回

────

dict with per-parameter results and overall verdict

"""

assert source_params.shape[1] == 6, "需要六个参数列"

assert target_params.shape[1] == 6, "需要六个参数列"



results = {}

passed_count = 0



for i, name in enumerate(param_names):

src = np.log10(np.clip(source_params[:, i], 1e-300, None))

tgt = np.log10(np.clip(target_params[:, i], 1e-300, None))



# KSG互信息

mi_val = ksg_mutual_information(src, tgt)



# 置换检验获得p值

def statistic(x, y):

return ksg_mutual_information(x, y)



# 简化置换检验

null_mis = []

for _ in range(min(n_permutations, 200)):

perm_tgt = np.random.permutation(tgt)

null_mis.append(ksg_mutual_information(src, perm_tgt))



p_val  = np.mean(np.array(null_mis) >= mi_val)

passed = (mi_val >= threshold_bits) and (p_val < 0.05)



results[name] = {

'mutual_information': mi_val,

'p_value': p_val,

'passed': passed,

'threshold': threshold_bits,

}



if passed:

passed_count += 1



results['overall'] = {

'n_passed': passed_count,

'n_total': 6,

'verdict': 'PASS' if passed_count >= 4 else 'FAIL',

'note': f'{passed_count}/6个参数通过功能等价性检验',

}



return results





def topological_alignment_test(source_params: np.ndarray,

target_params: np.ndarray,

threshold_r2: float = 0.75,

n_permutations: int = 1000) -> dict:

"""

拓扑对齐检验（Procrustes分析）



参数

────

source_params : 形状(n_timepoints, 6)的源域参数矩阵

target_params : 形状(n_timepoints, 6)的目标域参数矩阵

threshold_r2  : R²阈值（默认0.75）

n_permutations: 置换检验次数



返回

────

dict with Procrustes statistics and verdict

"""

# 标准化（零均值、单位方差）

def standardize(X: np.ndarray) -> np.ndarray:

mu  = X.mean(axis=0)

std = X.std(axis=0)

std[std < 1e-10] = 1.0

return (X - mu) / std



# 对参数取对数后标准化

log_src = standardize(

np.log10(np.clip(source_params, 1e-300, None))

)

log_tgt = standardize(

np.log10(np.clip(target_params, 1e-300, None))

)



# Procrustes对齐

mtx1, mtx2, disparity = procrustes(log_src, log_tgt)



# 计算R²

ss_tot  = np.sum((log_src - log_src.mean(axis=0))**2)

ss_res  = disparity * ss_tot

r2      = 1 - ss_res / ss_tot



# 置换检验

null_r2 = []

for _ in range(n_permutations):

perm_idx = np.random.permutation(len(log_tgt))

_, _, d_null = procrustes(log_src, log_tgt[perm_idx])

r2_null = 1 - d_null * ss_tot / ss_tot

null_r2.append(r2_null)



p_val  = np.mean(np.array(null_r2) >= r2)

passed = (r2 >= threshold_r2) and (p_val < 0.05)



return {

'procrustes_r2'    : r2,

'disparity'        : disparity,

'p_value'          : p_val,

'threshold'        : threshold_r2,

'passed'           : passed,

'verdict'          : 'PASS' if passed else 'FAIL',

'note'             : (

f"Procrustes R²={r2:.3f}（阈值{threshold_r2}），"

f"p={p_val:.3f}"

),

}





def full_crossdomain_validation(source_params: np.ndarray,

target_params: np.ndarray,

param_names: list = None) -> dict:

"""

完整跨域验证（KSG + Procrustes）



两个检验都通过才认定映射有效。



返回

────

dict with combined verdict

"""

if param_names is None:

param_names = ['D', 'C', 'S', 'F', 'W', 'M']



feq = functional_equivalence_test(

source_params, target_params, param_names

)

topo = topological_alignment_test(

source_params, target_params

)



both_passed = (feq['overall']['verdict'] == 'PASS' and

topo['verdict'] == 'PASS')



return {

'functional_equivalence': feq,

'topological_alignment' : topo,

'combined_verdict'      : 'VALID' if both_passed else 'INVALID',

'summary': (

f"功能等价性：{feq['overall']['verdict']} "

f"({feq['overall']['n_passed']}/6参数通过)\n"

f"拓扑对齐：{topo['verdict']} "

f"(R²={topo['procrustes_r2']:.3f})\n"

f"综合判定：{'映射有效' if both_passed else '映射无效'}"

),

}

### A.3.4 贝叶斯标定模块

"""

ici_calibration.py

ICI常数（k, alpha）的贝叶斯MCMC标定

"""



import numpy as np

import pymc as pm

import arviz as az

from typing import List, Tuple





def calibrate_ici_constants(

reference_systems: List[dict],

n_samples: int = 2000,

n_chains: int = 4,

target_accept: float = 0.9,

) -> dict:

"""

使用MCMC标定k和alpha



参数

────

reference_systems : 参照系统列表，每个系统为dict，包含：

- 'name': 系统名称

- 'D', 'C', 'S', 'F', 'W', 'M': 参数点估计

- 'ICI_observed': 从神经科学/行为数据估算的ICI值

- 'ICI_sigma': ICI观测值的对数标准差



返回

────

dict with posterior summaries

"""

# 提取数据

DCS_vals = np.array([

np.log10(sys['D'] * sys['C'] * sys['S'])

for sys in reference_systems

])

FWM_ratio_vals = np.array([

sys['F'] * sys['W'] * sys['M'] / 7.52e11

for sys in reference_systems

])

ICI_obs  = np.array([sys['ICI_observed'] for sys in reference_systems])

ICI_sigma = np.array([

sys.get('ICI_sigma', 0.1) for sys in reference_systems

])



# PyMC模型

with pm.Model() as model:

# 先验：k ~ Normal(1.26, 0.1)（中心在当前估计值附近）

k_param = pm.Normal('k', mu=1.26, sigma=0.1)



# 先验：log_alpha ~ Normal(log(1.02e5), 0.5)

log_alpha = pm.Normal('log_alpha',

mu=np.log(1.02e5),

sigma=0.5)

alpha_param = pm.Deterministic('alpha', pm.math.exp(log_alpha))



# ICI模型预测

emergence  = pm.math.sqrt(alpha_param * FWM_ratio_vals)

ICI_pred   = k_param * DCS_vals * (1 + emergence)



# 似然（对数尺度）

log_ICI_obs  = np.log(ICI_obs)

log_ICI_pred = pm.math.log(pm.math.abs_(ICI_pred) + 1e-10)



obs = pm.Normal('obs',

mu=log_ICI_pred,

sigma=ICI_sigma,

observed=log_ICI_obs)



# MCMC采样

trace = pm.sample(

draws=n_samples,

chains=n_chains,

target_accept=target_accept,

return_inferencedata=True,

progressbar=True,

)



# 提取后验统计

k_post     = trace.posterior['k'].values.flatten()

alpha_post = trace.posterior['alpha'].values.flatten()



# 收敛诊断

rhat     = az.rhat(trace)

ess      = az.ess(trace)

converged = (

all(rhat.data_vars[v].values < 1.01

for v in ['k', 'log_alpha']) and

all(ess.data_vars[v].values > 400

for v in ['k', 'log_alpha'])

)



return {

'k_mean'    : float(np.mean(k_post)),

'k_std'     : float(np.std(k_post)),

'k_ci95'    : (float(np.percentile(k_post, 2.5)),

float(np.percentile(k_post, 97.5))),

'alpha_mean': float(np.mean(alpha_post)),

'alpha_std' : float(np.std(alpha_post)),

'alpha_ci95': (float(np.percentile(alpha_post, 2.5)),

float(np.percentile(alpha_post, 97.5))),

'converged' : converged,

'r_hat_max' : float(max(

rhat.data_vars[v].values.max()

for v in ['k', 'log_alpha']

)),

'trace'     : trace,

}



## A.4 可视化模块

"""

ici_visualization.py

ICI框架可视化工具

"""



import numpy as np

import matplotlib.pyplot as plt

import matplotlib.patches as mpatches

from matplotlib.colors import LinearSegmentedColormap

import seaborn as sns

from typing import List, Optional



# ICI配色方案

ICI_COLORS = {

'green'  : '#2E7D32',

'yellow' : '#F9A825',

'orange' : '#E65100',

'red'    : '#B71C1C',

'blue'   : '#1565C0',

'gray'   : '#546E7A',

}





def plot_rt_timeseries(time_points: list,

Rt_values: np.ndarray,

Rt_lower: Optional[np.ndarray] = None,

Rt_upper: Optional[np.ndarray] = None,

system_name: str = "",

events: Optional[dict] = None,

figsize: tuple = (12, 5)) -> plt.Figure:

"""

绘制R(t)时间轨迹图，含临界区间和不确定性带



参数

────

time_points  : 时间点列表

Rt_values    : R(t)点估计

Rt_lower/upper : 68%置信带（可选）

events       : 标注事件字典 {time: label}

"""

fig, ax = plt.subplots(figsize=figsize)



# 绘制临界区间背景

ax.axhspan(-9.0,  5.0,  alpha=0.08, color=ICI_COLORS['green'],

label='强韧稳态 (R > −9)')

ax.axhspan(-10.0, -9.0, alpha=0.12, color=ICI_COLORS['yellow'],

label='亚稳临界带 (−10 ≤ R < −9)')

ax.axhspan(-15.0, -10.0, alpha=0.12, color=ICI_COLORS['red'],

label='崩溃区 (R < −10)')



# 临界线

ax.axhline(-9.0,  color=ICI_COLORS['yellow'], lw=1.5,

ls='--', alpha=0.7)

ax.axhline(-10.0, color=ICI_COLORS['red'], lw=2.0,

ls='--', alpha=0.9, label='临界阈值 R = −10')



# 不确定性带

if Rt_lower is not None and Rt_upper is not None:

ax.fill_between(time_points, Rt_lower, Rt_upper,

alpha=0.25, color=ICI_COLORS['blue'],

label='68%置信带')



# 主曲线（按区间着色）

for i in range(len(Rt_values) - 1):

Rt_mid = (Rt_values[i] + Rt_values[i+1]) / 2

color  = (ICI_COLORS['green']  if Rt_mid > -9.0 else

ICI_COLORS['yellow'] if Rt_mid > -10.0 else

ICI_COLORS['red'])

ax.plot(time_points[i:i+2], Rt_values[i:i+2],

color=color, lw=2.5, solid_capstyle='round')



# 数据点标记

for i, (t, r) in enumerate(zip(time_points, Rt_values)):

color = (ICI_COLORS['green']  if r > -9.0 else

ICI_COLORS['yellow'] if r > -10.0 else

ICI_COLORS['red'])

ax.scatter(t, r, color=color, s=60, zorder=5,

edgecolors='white', linewidths=1.0)



# 事件标注

if events:

for t, label in events.items():

ax.axvline(t, color='black', lw=1.0, ls=':',

alpha=0.6)

ax.annotate(label,

xy=(t, ax.get_ylim()[1] * 0.95),

xytext=(5, 0),

textcoords='offset points',

fontsize=9, color='black',

va='top', rotation=45)



# 格式化

ax.set_xlabel('时间', fontsize=12)

ax.set_ylabel('R(t)', fontsize=12)

ax.set_title(f'R(t) 时间轨迹：{system_name}', fontsize=13)

ax.legend(loc='upper right', fontsize=9,

framealpha=0.9, ncol=2)

ax.set_ylim(min(Rt_values) - 1.5, 2.0)

ax.grid(True, alpha=0.3, ls='--')

ax.spines['top'].set_visible(False)

ax.spines['right'].set_visible(False)



plt.tight_layout()

return fig





def plot_csd_dashboard(time_points: list,

Rt_series: np.ndarray,

csd_signals: dict,

system_name: str = "",

figsize: tuple = (14, 10)) -> plt.Figure:

"""

绘制临界慢化仪表盘（4面板）



面板：R(t)轨迹 | 滚动方差 | 滚动AR(1) | 滚动Hurst指数

"""

fig, axes = plt.subplots(2, 2, figsize=figsize)

fig.suptitle(f'临界慢化预警仪表盘：{system_name}',

fontsize=14, y=0.98)



# 面板1：R(t)

ax = axes[0, 0]

ax.axhline(-9.0,  color=ICI_COLORS['yellow'], lw=1.5, ls='--')

ax.axhline(-10.0, color=ICI_COLORS['red'],    lw=2.0, ls='--')

ax.plot(time_points, Rt_series,

color=ICI_COLORS['blue'], lw=2)

ax.set_title('R(t) 轨迹', fontsize=11)

ax.set_ylabel('R(t)')

ax.grid(True, alpha=0.3)



# 面板2：滚动方差

ax = axes[0, 1]

var_s   = csd_signals['variance']

tau_var = csd_signals['tau_var']

color   = (ICI_COLORS['red'] if tau_var > 0.4

else ICI_COLORS['gray'])

ax.plot(time_points, var_s, color=color, lw=2)

ax.set_title(f'滚动方差  τ={tau_var:.2f}  '

f'p={csd_signals["p_var"]:.3f}',

fontsize=10)

ax.set_ylabel('方差')

ax.grid(True, alpha=0.3)



# 面板3：滚动AR(1)

ax = axes[1, 0]

ar1_s   = csd_signals['ar1']

tau_ar1 = csd_signals['tau_ar1']

color   = (ICI_COLORS['red'] if tau_ar1 > 0.4

else ICI_COLORS['gray'])

ax.plot(time_points, ar1_s, color=color, lw=2)

ax.axhline(1.0, color='black', lw=1, ls='--', alpha=0.5)

ax.set_title(f'滚动AR(1)  τ={tau_ar1:.2f}  '

f'p={csd_signals["p_ar1"]:.3f}',

fontsize=10)

ax.set_ylabel('AR(1)系数')

ax.set_ylim(-0.2, 1.2)

ax.grid(True, alpha=0.3)



# 面板4：滚动Hurst指数

ax = axes[1, 1]

hurst_s = csd_signals['hurst']

color   = (ICI_COLORS['red'] if csd_signals['hurst_trans']

else ICI_COLORS['gray'])

ax.plot(time_points, hurst_s, color=color, lw=2)

ax.axhline(0.5, color='black', lw=1.5, ls='--',

label='H=0.5（随机游走）')

ax.fill_between(time_points,

np.where(~np.isnan(hurst_s), hurst_s, 0.5),

0.5,

where=~np.isnan(hurst_s),

alpha=0.2, color=color)

trans_label = '⚠ 检测到Hurst跃迁' if csd_signals['hurst_trans'] else '未检测到跃迁'

ax.set_title(f'滚动Hurst指数  {trans_label}', fontsize=10)

ax.set_ylabel('H')

ax.set_ylim(0.0, 1.0)

ax.legend(fontsize=9)

ax.grid(True, alpha=0.3)



# 预警分数标注

score = csd_signals['alert_score']

level = csd_signals['alert_level']

color = ICI_COLORS.get(level, 'gray')

fig.text(0.5, 0.01,

f'综合预警分数：{score}/8  |  级别：{level.upper()}',

ha='center', fontsize=12,

color=color, fontweight='bold')



for ax in axes.flat:

ax.spines['top'].set_visible(False)

ax.spines['right'].set_visible(False)



plt.tight_layout(rect=[0, 0.04, 1, 1])

return fig





def plot_ici_spectrum(systems: list,

ici_values: list,

ici_errors: Optional[list] = None,

figsize: tuple = (10, 7)) -> plt.Figure:

"""

绘制跨系统ICI谱图



参数

────

systems    : 系统名称列表

ici_values : ICI值列表

ici_errors : ICI不确定性（可选）

"""

fig, ax = plt.subplots(figsize=figsize)



# 按ICI排序

order = np.argsort(ici_values)

sys_sorted = [systems[i] for i in order]

ici_sorted = [ici_values[i] for i in order]

err_sorted = ([ici_errors[i] for i in order]

if ici_errors else None)



# 颜色映射

def ici_color(val):

if val < 50:    return ICI_COLORS['gray']

if val < 5000:  return ICI_COLORS['blue']

if val < 8000:  return '#7B1FA2'

return '#B71C1C'



colors = [ici_color(v) for v in ici_sorted]



# 绘制水平条形图（对数轴）

y_pos = np.arange(len(sys_sorted))

bars  = ax.barh(y_pos, ici_sorted, color=colors,

alpha=0.8, height=0.7)



if err_sorted:

ax.errorbar(ici_sorted, y_pos,

xerr=err_sorted, fmt='none',

color='black', capsize=4, lw=1.5)



# 分层背景

ax.axvspan(0,    50,   alpha=0.04, color='gray',

label='初级层 (0–50)')

ax.axvspan(50,   5000, alpha=0.04, color=ICI_COLORS['blue'],

label='中级层 (50–5000)')

ax.axvspan(5000, 8000, alpha=0.04, color='purple',

label='高级层 (5000–8000)')

ax.axvspan(8000, 12000, alpha=0.04, color='red',

label='顶层 (8000+)')



ax.set_yticks(y_pos)

ax.set_yticklabels(sys_sorted, fontsize=10)

ax.set_xscale('log')

ax.set_xlabel('ICI值（对数轴）', fontsize=12)

ax.set_title('系统复杂性谱：ICI比较', fontsize=13)

ax.legend(loc='lower right', fontsize=9)

ax.grid(True, axis='x', alpha=0.3, which='both')

ax.spines['top'].set_visible(False)

ax.spines['right'].set_visible(False)



plt.tight_layout()

return fig



## A.5 完整使用示例

"""

example_analysis.py

完整分析示例：以历史罗马帝国为例

"""



from ici_core import ICIParameters, compute_ici, compute_ici_timeseries

from ici_csd import compute_csd_signals

from ici_visualization import (plot_rt_timeseries,

plot_csd_dashboard,

plot_ici_spectrum)

import numpy as np





def rome_analysis():

"""罗马帝国R(t)完整分析示例"""



print("=" * 60)

print("ICI分析：罗马帝国兴衰轨迹")

print("=" * 60)



# 定义各时间点参数

# ──────────────────

time_labels = [

"鼎盛期\n(公元100年)",

"三世纪危机\n(公元235年)",

"戴克里先改革\n(公元285年)",

"后期帝国\n(公元400年)",

"西罗马灭亡\n(公元476年)",

]



params_series = [

ICIParameters(

D=50, C=5e7, S=1e6, F=20, W=4.8e-13, M=1e3,

D_sigma=0.17, C_sigma=0.22, S_sigma=0.30,

F_sigma=0.35, W_sigma=0.30, M_sigma=0.39,

system_name="罗马帝国·鼎盛期",

domain='historical_empire',

),

ICIParameters(

D=65, C=5e7, S=8e5, F=12, W=2.4e-13, M=6e2,

D_sigma=0.17, C_sigma=0.22, S_sigma=0.30,

F_sigma=0.35, W_sigma=0.30, M_sigma=0.39,

system_name="罗马帝国·三世纪危机",

domain='historical_empire',

),

ICIParameters(

D=55, C=5e7, S=9e5, F=18, W=4.0e-13, M=9e2,

D_sigma=0.17, C_sigma=0.22, S_sigma=0.30,

F_sigma=0.35, W_sigma=0.30, M_sigma=0.39,

system_name="罗马帝国·戴克里先改革",

domain='historical_empire',

),

ICIParameters(

D=80, C=3e7, S=4e5, F=6, W=9.5e-14, M=3e2,

D_sigma=0.17, C_sigma=0.22, S_sigma=0.30,

F_sigma=0.35, W_sigma=0.30, M_sigma=0.39,

system_name="罗马帝国·后期",

domain='historical_empire',

),

ICIParameters(

D=90, C=1.5e7, S=2e5, F=3, W=3.2e-14, M=1e2,

D_sigma=0.17, C_sigma=0.22, S_sigma=0.30,

F_sigma=0.35, W_sigma=0.30, M_sigma=0.39,

system_name="罗马帝国·西罗马灭亡",

domain='historical_empire',

),

]



# 计算时间序列

# ─────────────

results = compute_ici_timeseries(params_series, reference_index=0)



Rt_vals  = np.array([r.R_t       for r in results])

Rt_lower = np.array([r.R_t_lower for r in results])

Rt_upper = np.array([r.R_t_upper for r in results])

ICI_vals = np.array([r.ICI       for r in results])



print("\n参数与结果汇总：")

print("-" * 60)

for i, (label, res) in enumerate(zip(time_labels, results)):

label_clean = label.replace('\n', ' ')

print(f"{label_clean:20s} | "

f"ICI={res.ICI:6.1f} | "

f"R(t)={res.R_t:+.2f} | "

f"预警：{res.alert_level.upper()}")



# 绘制R(t)轨迹图

# ───────────────

fig1 = plot_rt_timeseries(

time_points = list(range(len(time_labels))),

Rt_values   = Rt_vals,

Rt_lower    = Rt_lower,

Rt_upper    = Rt_upper,

system_name = "罗马帝国",

events      = {1: "三世纪危机", 2: "改革", 4: "灭亡"},

)

fig1.savefig('rome_rt_trajectory.pdf',

dpi=150, bbox_inches='tight')

print("\n已保存：rome_rt_trajectory.pdf")



# 临界慢化分析（需要更长时间序列才有意义，此处演示）

# ────────────────────────────────────────────────

if len(Rt_vals) >= 5:

csd = compute_csd_signals(Rt_vals)

print(f"\n临界慢化信号：")

print(f"  方差趋势 τ = {csd['tau_var']:.2f}  "

f"(p = {csd['p_var']:.3f})")

print(f"  AR(1)趋势 τ = {csd['tau_ar1']:.2f}  "

f"(p = {csd['p_ar1']:.3f})")

print(f"  Hurst跃迁检测：{csd['hurst_trans']}")

print(f"  综合预警分数：{csd['alert_score']}/8"

f"  级别：{csd['alert_level'].upper()}")



return results





def biological_spectrum_example():

"""生物系统ICI谱图示例"""



systems = [

"大肠杆菌",

"酿酒酵母",

"拟南芥",

"线虫(C.elegans)",

"水螅神经网",

"果蝇",

"小鼠皮层",

"大鼠海马CA1",

"猕猴前额叶",

"人类海马CA1",

]



# 书中各系统的ICI估算值

ici_values = [20.3, 28.4, 36.8, 89.4, 124.6,

1247, 5001, 6723, 7891, 10000]

ici_errors = [3.2, 4.5, 6.1, 15.2, 22.3,

210, 850, 1150, 1340, 1700]



fig = plot_ici_spectrum(systems, ici_values, ici_errors)

fig.savefig('biological_ici_spectrum.pdf',

dpi=150, bbox_inches='tight')

print("已保存：biological_ici_spectrum.pdf")





if __name__ == '__main__':

results = rome_analysis()

biological_spectrum_example()

print("\n所有分析完成。")



## A.6 参数记录标准表

ICI参数记录表 v3.0

═══════════════════════════════════════════════════════════



基本信息

────────

系统名称：______________________________________________

分析日期：______________________________________________

分析人员：______________________________________________

域类型：  □生物 □历史 □经济 □生态 □AI □其他_________

时间节点：______________________________________________



参数估算

────────

参数  | 点估计 | 对数误差σ | 数据层级 | 来源

──────|────────|──────────|──────────|──────────────────

D     |        |          | □1 □2 □3 |

C     |        |          | □1 □2 □3 |

S     |        |          | □1 □2 □3 |

F     |        |          | □1 □2 □3 |

W     |        |          | □1 □2 □3 |

M     |        |          | □1 □2 □3 |



数据层级：1=直接测量(±20-30%) 2=间接推算(±40-60%)

3=专家判断(±60-90%)



计算结果

────────

lg(DCS) = _____________

FWM/FWM_h = _____________

ICI（点估计） = _____________

ICI（68%CI）  = [________, ________]

ICI（95%CI）  = [________, ________]

R(t)（绝对值）= _____________

R(t)（域内相对）= _____________（基准：_______________）

预警级别：□绿色 □黄色 □橙色 □红色



验证检查

────────

□ 参数值均为正数

□ ICI在域的合理范围内

□ R(t)绝对值与已知案例一致

□ 不确定性已充分量化

□ 数据来源已完整记录



备注

────

_______________________________________________________

_______________________________________________________



数据文件路径（原始数据+计算脚本）：

_______________________________________________________



