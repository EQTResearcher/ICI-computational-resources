# SOP-02：R(t)时间序列监测协议



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

