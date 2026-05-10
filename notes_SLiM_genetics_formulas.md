# SLiM 软件遗传学公式、计算原理与随机建模过程全面总结

## 一、软件概述

本软件是 **SLiM**（Simulating Evolutionary and Ecological Dynamics），一个正向时间（forward-time）、个体水平（individual-based）的种群遗传学模拟器，支持 Wright-Fisher（WF）和非 WF（nonWF）模型，以及多物种（multi-species）模型，并可记录完整系谱（tree sequence recording）。

---

## 二、突变过程（Mutation）

### 2.1 突变事件计数——泊松过程

每次减数分裂中，新突变数量从泊松分布抽样：

$$N_{mut} \sim \text{Poisson}(\mu_{total})$$

其中 $\mu_{total}$ 为整条染色体上所有基因组元素区间的突变率加权总和（突变率图与基因组元素图取交集后的分段积分）。

### 2.2 突变率重参数化（Position Uniquing 校正）

从 SLiM 3.5 开始，同一配子中不允许同一位置出现两个突变（unique 操作）。为保证实现突变率等于用户请求率，对每段区间进行重参数化：

$$\mu_{adj} = -\ln(1 - \mu_{req})$$

原理：$P[\text{Poisson}(\lambda) \geq 1] = 1 - e^{-\lambda} = p$，故 $\lambda = -\ln(1-p)$。

### 2.3 突变位置抽样

突变位置从各子区间的加权离散分布抽取（Walker/Vose 别名法，`gsl_ran_discrete`），权重为区间率 × 区间长度。

### 2.4 适应度效应分布（DFE）

每类突变（`MutationType`）的选择系数 $s$ 从以下分布之一抽取：

| 类型 | 分布 | 参数 | GSL 实现 |
|------|------|------|---------|
| `f` | 固定值 | $s_0$ | — |
| `g` | Gamma | 均值 $\bar{s}$, 形状 $\alpha$ | `gsl_ran_gamma(rng, α, s̄/α)` |
| `e` | 指数 | 均值 $\bar{s}$ | `gsl_ran_exponential(rng, s̄)` |
| `n` | 正态 | 均值 $\mu$, 标准差 $\sigma$ | `gsl_ran_gaussian(rng, σ) + μ` |
| `w` | Weibull | 尺度 $\lambda$, 形状 $k$ | `gsl_ran_weibull(rng, λ, k)` |
| `p` | Laplace | 均值 $\mu$, 尺度 $b$ | `gsl_ran_laplace(rng, b) + μ` |
| `s` | 脚本自定义 | Eidos 代码 | — |

---

## 三、重组过程（Recombination）

### 3.1 交叉互换断点模型（Crossover Breakpoint Model）

断点数从泊松分布抽样：

$$N_{break} \sim \text{Poisson}(\lambda_{total})$$

用户指定的重组率 $p_i$ 为相邻碱基间发生交叉互换的概率，经重参数化校正（原因同突变）：

$$\lambda_i = -\ln(1 - p_i)$$

约束：所有区间重组率 $p_i \in [0, 0.5]$（0.5 代表完全独立；不允许"反连锁"）。

断点位置从 `gsl_ran_discrete` 按加权离散分布抽取，抽取后排序并去重（unique），保证同一区间最多有一个断点。

### 3.2 联合概率优化（Joint Poisson Optimization）

为了高效跳过"无突变且无断点"的常见情况，预计算联合概率：

$$P_0^{mut} = e^{-\mu_{total}}, \quad P_0^{break} = e^{-\lambda_{total}}$$

$$P_{\text{both 0}} = P_0^{mut} \cdot P_0^{break}$$

$$P_{\text{mut=0 only}} = P_0^{mut} \cdot (1 - P_0^{break})$$

$$P_{\text{break=0 only}} = (1 - P_0^{mut}) \cdot P_0^{break}$$

这三个区间的累积概率用于单次均匀随机数判断，实现高效分支。

### 3.3 双链断裂模型（DSB Gene Conversion Model）

可选替代模型，通过 `initializeGeneConversion()` 激活：

- DSB 点位置抽取同交叉互换断点
- 每个 DSB 向两侧随机延伸，形成基因转变片段（gene conversion tract）：

$$\text{extent}_1, \text{extent}_2 \sim \text{Geometric}\!\left(\frac{1}{L/2}\right)$$

其中 $L$ 为用户指定的平均片段长度，几何分布参数为 $1/(L/2)$。

- 每个 DSB 独立判断：
  - **非交叉互换概率**：$P(\text{noncrossover}) = \text{nonCrossoverFraction}$
  - **简单转变概率**：$P(\text{simple}) = \text{simpleConversionFraction}$
  - 复杂（非简单）转变产生杂合双链区（heteroduplex regions）

- **gBGC（GC 偏向性基因转变）**：`mismatch_repair_bias` 参数，偏向 G/C 碱基的错配修复，适用于核苷酸模型。

---

## 四、适应度计算（Fitness）

### 4.1 个体适应度——乘积模型

个体适应度 $W$ 为所有突变效应的乘积：

$$W = \prod_{i} w_i$$

对**二倍体**个体（单倍体组 H1, H2）：

- **杂合突变**（仅存在于一条单倍体组）：
  $$w_i = \max(0,\ 1 + h \cdot s)$$
  其中 $h$ = 显性系数（dominance coefficient），$s$ = 选择系数。

- **纯合突变**（两条单倍体组均含有）：
  $$w_i = \max(0,\ 1 + s)$$

- **半合（Hemizygous）突变**（如雄性 X 染色体上的突变）：
  $$w_i = \max(0,\ 1 + h_{hemi} \cdot s)$$
  其中 $h_{hemi}$ = hemizygous 显性系数。

两个值被预缓存为：

- `cached_one_plus_sel_ = max(0, 1+s)`
- `cached_one_plus_dom_sel_ = max(0, 1+h·s)`

### 4.2 适应度缩放（Fitness Scaling）

总适应度 = 个体突变组合效应 × 个体缩放因子 × 亚群缩放因子：

$$W_{final} = W_{genotype} \times f_{individual} \times f_{subpopulation}$$

允许通过 `fitnessScaling` 属性实现密度依赖或环境依赖的适应度缩放。

### 4.3 亲本采样（WF 模型）

在 WF 模型中，每个子代从亲本群体中依照适应度比例采样，使用 Vose-Alias O(1) 别名法（`gsl_ran_discrete`）实现。

---

## 五、种群遗传学统计量

### 5.1 核苷酸多样性 π

$$\hat{\pi} = \frac{\sum_j c_j(n - c_j)}{\binom{n}{2} \cdot L} = \frac{2\sum_j p_j(1-p_j)}{L}$$

$c_j$ = 第 $j$ 位点的衍生等位基因计数，$n$ = 单倍体组数，$L$ = 序列长度，$p_j = c_j/n$。  
等价于 Korunes & Samuk 2021 公式 1。

### 5.2 Watterson's θ

$$\hat{\theta}_W = \frac{k}{a_n \cdot L}, \quad a_n = \sum_{i=1}^{n-1} \frac{1}{i}$$

$k$ = 分离位点数（segregating sites），$n$ = 单倍体组数。

### 5.3 Tajima's D

$$D = \frac{(\hat{\pi} - \hat{\theta}_W / L) \cdot L}{\sqrt{e_1 k + e_2 k(k-1)}}$$

中间变量：

$$a_1 = \sum_{i=1}^{n-1}\frac{1}{i}, \quad a_2 = \sum_{i=1}^{n-1}\frac{1}{i^2}$$

$$b_1 = \frac{n+1}{3(n-1)}, \quad b_2 = \frac{2(n^2+n+3)}{9n(n-1)}$$

$$c_1 = b_1 - \frac{1}{a_1}, \quad c_2 = b_2 - \frac{n+2}{a_1 n} + \frac{a_2}{a_1^2}$$

$$e_1 = \frac{c_1}{a_1}, \quad e_2 = \frac{c_2}{a_1^2 + a_2}$$

### 5.4 FST（群体间分化）

使用 Hudson 式估计量：

$$F_{ST} = 1 - \frac{\overline{H}_S}{\overline{H}_T}$$

其中：

- $\bar{p} = (p_1 + p_2)/2$（两群体均值）
- $H_T = 2\bar{p}(1-\bar{p})$（总体期望杂合度）
- $H_S = p_1(1-p_1) + p_2(1-p_2)$（亚群期望杂合度之和）
- $\overline{H}_S$, $\overline{H}_T$ 为跨所有位点的均值

### 5.5 Dxy（绝对遗传分歧）

Nei 定义（未除以序列长度版本）：

$$D_{xy} = \frac{\sum_j \left[ c_{1j}(n_2 - c_{2j}) + c_{2j}(n_1 - c_{1j}) \right]}{2 n_1 n_2}$$

可通过 `normalize=T` 除以序列长度得到每碱基 Dxy。

### 5.6 连锁不平衡（LD）

- **D 系数**：$D = p_{12} - p_1 p_2$（$p_{12}$ = 同时含两突变的单倍型频率）
- **r²**：$r^2 = \dfrac{D^2}{p_1 p_2 (1-p_1)(1-p_2)}$
- **r（有符号）**：$r = \dfrac{D}{\sqrt{p_1 p_2 (1-p_1)(1-p_2)}}$

### 5.7 杂合度（Heterozygosity）

$$H = \frac{2\sum_j p_j(1-p_j)}{L}$$

### 5.8 近交负荷（Inbreeding Load, B）

Morton et al. (1956) 公式：

$$B = \sum_j q_j s_j - \sum_j q_j^2 s_j - 2\sum_j q_j(1-q_j) s_j h_j$$

其中 $q_j$ = 有害等位基因频率，$s_j = -\text{selectionCoeff}$（大于零），$h_j$ = 显性系数。

### 5.9 F_ROH（基于纯合延伸段的近交系数）

$$F_{ROH} = \frac{\sum_k L_k^{ROH}}{L_{total}}$$

ROH 段（runs of homozygosity）定义为相邻杂合位点之间的纯合区间，仅计入长度超过最小阈值（默认 1 Mb）的段；$L_{total}$ 为二倍体染色体总长度。

### 5.10 加性遗传方差（VA）

$$V_A = \text{Var}\!\left(\sum_{i \in \text{mutType}} s_i \cdot \text{dosage}_i\right)$$

即指定突变类型的"剂量加权选择系数之和"在个体间的方差。

### 5.11 等位基因频率谱（SFS）

- 将所有分离位点按衍生等位基因频率（或样本计数）分箱（histogram）
- 支持展开式（unfolded）和折叠式（folded，将对称频率合并）
- 输出密度（density）或计数（count）

---

## 六、核苷酸进化模型

### 6.1 Jukes-Cantor（JC）模型

4×4 突变矩阵（行=起始碱基，列=目标碱基，顺序 A/C/G/T）：

$$M_{JC} = \begin{pmatrix} 0 & \alpha & \alpha & \alpha \\ \alpha & 0 & \alpha & \alpha \\ \alpha & \alpha & 0 & \alpha \\ \alpha & \alpha & \alpha & 0 \end{pmatrix}, \quad 3\alpha \leq 1$$

所有替换率相等；`mmJukesCantor(alpha)`。

### 6.2 Kimura 两参数（K2P）模型

$$M_{K2P} = \begin{pmatrix} 0 & \beta & \alpha & \beta \\ \beta & 0 & \beta & \alpha \\ \alpha & \beta & 0 & \beta \\ \beta & \alpha & \beta & 0 \end{pmatrix}, \quad \alpha + 2\beta \leq 1$$

$\alpha$ = 转换率（transition: A↔G, C↔T），$\beta$ = 颠换率（transversion）；`mmKimura(alpha, beta)`。

### 6.3 通用 4×4 突变矩阵

`initializeGenomicElementType()` 中通过 `mutationMatrix` 参数指定完整的 4×4（或 16-element）概率矩阵，支持完全非对称、上下文依赖等复杂模型。

### 6.4 扩展为 256 元素矩阵

`mm16To256()` 将 4×4 矩阵扩展为 256 元素矩阵，支持三核苷酸背景依赖的突变模型（类似 SBS 模式）。

---

## 七、空间互作核函数（Spatial Interaction Kernels）

用于空间模型中的交配、竞争、合作等互作强度，以距离 $d$ 为自变量：

| 核函数 | 代码 | 公式 |
|--------|------|------|
| 固定强度 | `f` | $f_{max}$ |
| 线性衰减 | `l` | $f_{max} (1 - d/d_{max})$ |
| 指数衰减 | `e` | $f_{max} e^{-\lambda d}$ |
| 高斯（正态）| `n` | $f_{max} e^{-d^2/(2\sigma^2)}$ |
| Cauchy | `c` | $f_{max} / (1 + (d/\lambda)^2)$ |
| Student's t | `t` | $f_{max} / (1 + (d/\tau)^2/\nu)^{(\nu+1)/2}$ |

空间近邻搜索使用 k-d tree 数据结构，支持 1D/2D/3D 和周期边界条件。

---

## 八、随机模拟的建模过程

### 8.1 随机数生成基础设施

- 主 RNG：GSL Mersenne Twister（MT19937）
- 线程局部 RNG：OpenMP 并行时每个线程独立 RNG
- 自定义快速泊松：`Eidos_FastRandomPoisson` / `Eidos_FastRandomPoisson_PRECALCULATE`（预计算 $e^{-\lambda}$）

### 8.2 Wright-Fisher 模型每代流程

```
1. 计算亲本适应度 W_i（乘积模型 + 缩放因子）
2. 构建加权采样查找表（Vose-Alias，O(1)）
3. 生成子代：
   for each offspring:
     a. 按适应度采样父本/母本
     b. 减数分裂（meiosis）：
        - 抽 N_mut ~ Poisson(μ_total)
        - 抽 N_break ~ Poisson(λ_total) [or DSB model]
        - 抽突变位置（加权离散分布 + 去重）
        - 抽重组断点位置（加权离散分布 + 去重）
        - 两条亲本单倍体组按断点交换片段
        - 在各突变位置从 DFE 抽取 s，生成新突变对象
4. 突变注册/固定检测（freq=1 的突变转为替代 Substitution）
5. 子代替换亲代（离散世代）
```

### 8.3 非 WF（nonWF）模型

- 允许世代交叠；个体生存/繁殖通过用户脚本（`reproduction()`/`survival()` callbacks）控制
- 种群大小可动态变化（增删个体而非整代替换）
- 子代可为同一亲本的多个后代
- 死亡概率可基于适应度或自定义公式

### 8.4 树序列记录（Tree Sequence Recording）

- 记录每个新单倍体组的亲本关系（edges）
- 在模拟过程中周期性运行 Simplification，删除非祖先节点
- 最终输出 `.trees` 文件，兼容 tskit/msprime 工具链
- 支持连接到追溯模拟（coalescent simulation）以填充到 MRCA

### 8.5 突变运行优化（Mutation Run Caching）

- 每条单倍体组按位置分段为"突变运行"（mutation runs）
- 维护"非中性突变"缓存，跳过 $s=0$ 突变的适应度贡献
- 自适应调整突变运行数目（"mutation run experiments"），最小化内存和计算开销

### 8.6 多物种（Multi-species）社区模型

- `Community` 协调多个物种（`Species`）的生命周期排序
- 每个物种独立维护突变类型、染色体、亚群、重组率等
- 跨物种互作通过 `InteractionType` 实现（如捕食、竞争），但遗传信息不跨物种传递

---

## 总结表

| 模块 | 核心公式/方法 | 随机分布 |
|------|------------|---------|
| 突变数目 | $N_{mut} \sim \text{Poisson}(\mu_{total})$ | 泊松分布 |
| 突变位置 | 加权离散分布 | Walker 别名法 |
| 选择系数 | DFE：固定/Gamma/指数/正态/Weibull/Laplace | 各自对应分布 |
| 重组断点数 | $N_{break} \sim \text{Poisson}(\lambda_{total})$ | 泊松分布 |
| 重组断点位置 | 加权离散分布 | Walker 别名法 |
| 基因转变片段 | $\text{extent} \sim \text{Geometric}(1/(L/2))$ | 几何分布 |
| 适应度 | $W = \prod_i (1 + h_i s_i)$ 或 $(1 + s_i)$ | 乘积模型 |
| 亲本选择 | 适应度比例采样 | Vose-Alias O(1) |
| π | $2\sum p_j(1-p_j)/L$ | — |
| Watterson's θ | $k/(a_n L)$ | — |
| Tajima's D | $(\hat\pi - \hat\theta_W)/\sqrt{Var}$ | — |
| FST | $1 - \overline{H}_S/\overline{H}_T$ | — |
| Dxy | $\sum c_{1j}(n_2-c_{2j}) + \ldots\ /\ 2n_1n_2$ | — |
| LD D/r² | $D=p_{12}-p_1p_2$；$r^2=D^2/(p_1p_2(1-p_1)(1-p_2))$ | — |
| 近交负荷 B | Morton et al. 1956 | — |
| F_ROH | $\sum L_{ROH}/L_{total}$ | — |
| JC 矩阵 | 均等替换率 $\alpha$ | — |
| K2P 矩阵 | 转换 $\alpha$，颠换 $\beta$ | — |
| 空间核 | 固定/线性/指数/高斯/Cauchy/Student-t | — |
