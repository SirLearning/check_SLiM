# `calcInbreedingLoad(B)` 从随机数到最终结果的完整计算链

> 重点：模拟过程中低适应度个体的过滤参数作用原理

---

## 一、全局视角：三条平行的数据流

```
[随机数流 1] → 突变产生 → s 值确定 → 突变对象存活于单倍体组
[随机数流 2] → 亲本采样（适应度加权） → 等位基因频率演化
[随机数流 3] → 重组断点 → 单倍体组组合传递

三条流的交汇点：每代末的等位基因频率 q_j
最终汇入 B 公式：B = Σq_j·s_j - Σq_j²·s_j - 2Σq_j(1-q_j)·s_j·h_j
```

---

## 二、第一阶段：突变产生（随机数 → 突变对象 → s 值）

### 步骤 1：每代每条染色体，对每个子代配子，抽取突变数

```c
// chromosome.cpp，DrawMutationAndBreakpointCounts()
N_mut ~ Poisson(μ_total)
```

`μ_total` 是各基因组元素区间的 `rate × length` 之和（已做 uniquing 重参数化：`μ_adj = -ln(1-μ_req)`）。

**快速优化分支**（`_InitializeJointProbabilities`）：

```
预计算：P₀_mut  = exp(-μ_total)
预计算：P₀_break = exp(-λ_total)

抽一个 U(0,1)：
  U < P₀_mut × P₀_break          → 无突变无断点（最常见，直接跳过）
  else if U < P₀_mut             → 只有断点，无突变
  else if U < P₀_mut + (1-P₀_mut)×P₀_break → 只有突变，无断点
  else                           → 两者都有（需分别 Poisson 抽数）
```

### 步骤 2：抽取突变位置

```c
// chromosome.cpp，DrawSortedUniquedMutationPositions()
int interval = gsl_ran_discrete(rng_gsl, lookup_mutation_H_);  // O(1) 别名法
position = subrange.start + Eidos_rng_interval_uint64(rng_64, subrange_length);
```

多个位置经 `std::sort` + `std::unique` 去重（同一配子中每个位置只能有一个突变）。

### 步骤 3：从 DFE 抽取选择系数 s

```c
// mutation_type.cpp，DrawSelectionCoefficient()
switch (dfe_type_) {
  case kFixed:       s = dfe_parameters_[0];
  case kGamma:       s = gsl_ran_gamma(rng_gsl, shape, mean/shape);
  case kExponential: s = gsl_ran_exponential(rng_gsl, mean);
  case kNormal:      s = gsl_ran_gaussian(rng_gsl, sd) + mean;
  case kWeibull:     s = gsl_ran_weibull(rng_gsl, scale, shape);
  case kLaplace:     s = gsl_ran_laplace(rng_gsl, scale) + mean;
}
```

### 步骤 4：缓存适应度乘子

```c
// population.cpp，ValidateMutationFitnessCaches()
mut->cached_one_plus_sel_               = max(0.0, 1.0 + s)           // 纯合用
mut->cached_one_plus_dom_sel_           = max(0.0, 1.0 + h*s)         // 杂合用
mut->cached_one_plus_hemizygousdom_sel_ = max(0.0, 1.0 + h_hemi*s)    // 半合用
```

此时突变对象被添加到子代单倍体组的 `MutationRun` 中，加入全局 `mutation_registry_`（引用计数 refcount=1），等待自然选择的审判。

---

## 三、第二阶段：适应度计算——每个个体的"分数"

每代开始时，`Population::RecalculateFitness()` → `Subpopulation::UpdateFitness()` 为所有亲本计算适应度 $W_i$。

### 步骤 5：扫描单倍体组，计算个体适应度乘积

核心函数：`_Fitness_DiploidChromosome<>()` （`subpopulation.cpp`）

```
w = 1.0

对每条染色体的每个 MutationRun（按位置排序）：
  同时扫描 haplosome1 和 haplosome2 的突变列表（双指针归并）：

  ┌─ 情况 A：突变只在 haplosome1（杂合）─────────────────┐
  │  w *= mutation->cached_one_plus_dom_sel_              │  ← w = w × (1 + h·s)
  └──────────────────────────────────────────────────────┘

  ┌─ 情况 B：突变只在 haplosome2（杂合）─────────────────┐
  │  w *= mutation->cached_one_plus_dom_sel_              │  ← w = w × (1 + h·s)
  └──────────────────────────────────────────────────────┘

  ┌─ 情况 C：两条单倍体组都有同一突变（纯合）────────────┐
  │  w *= mutation->cached_one_plus_sel_                  │  ← w = w × (1 + s)
  └──────────────────────────────────────────────────────┘

  ┌─ 情况 D：haplosome1 为 null（半合）──────────────────┐
  │  w *= mutation->cached_one_plus_hemizygousdom_sel_    │  ← w = w × (1 + h_hemi·s)
  └──────────────────────────────────────────────────────┘

  如果 w ≤ 0.0，立即 return 0.0（短路优化）
```

**关键：非中性缓存优化**

`MutationRun` 维护 `nonneutral_mutations_` 缓存，只包含 $s \neq 0$ 的突变。当 $h=0.5$（加性）且 $s=0$ 时，`(1+h·s)=1.0`，乘以 1.0 无贡献，直接跳过，使纯中性模型以 $O(1)$ 计算适应度。

### 步骤 6：适应度缩放

```c
// subpopulation.cpp，UpdateFitness()
double fitness = individual->fitness_scaling_;    // 个体缩放因子（默认 1.0）
fitness *= FitnessOfParent_TEMPLATED(index, ...); // 上面的乘积 w
// fitnessEffect() 回调可再乘以额外因子
// subpop_fitness_scaling_ 代表亚群密度依赖因子
```

最终 `cached_fitness_UNSAFE_ = fitness`，**这是低适应度过滤的核心参数**。

---

## 四、第三阶段：亲本采样——低适应度个体的过滤机制

这是自然选择真正发生的地方。

### 步骤 7：构建加权查找表（Vose-Alias 法）

```c
// subpopulation.cpp，UpdateWFFitnessBuffers()
// 把每个个体的 cached_fitness_UNSAFE_ 写入 cached_parental_fitness_[] 数组

// 然后用 Walker-Vose Alias 方法预处理：
lookup_parent_ = gsl_ran_discrete_preproc(N, cached_parental_fitness_);
```

**Vose-Alias 方法原理**：

- 预处理 O(N)；单次抽样 O(1)
- 将 N 个概率值 $w_i / \sum w_j$ 组织成一个包含等高矩形的表
- 每次只需两次随机数：一次选列（整数），一次决定取本列还是别名列（浮点比较）

### 步骤 8：为每个子代采样亲本（关键过滤点）

```c
// population.cpp，EvolveSubpopulation()
for each child:
    parent1 = source_subpop.DrawParentUsingFitness(rng_state);
    parent2 = source_subpop.DrawParentUsingFitness(rng_state);
```

```c
// subpopulation.cpp（内联）
slim_popsize_t Subpopulation::DrawParentUsingFitness(Eidos_RNG_State *rng_state) {
    if (lookup_parent_)
        return (slim_popsize_t)gsl_ran_discrete(rng_gsl, lookup_parent_);
    else
        return (slim_popsize_t)(Eidos_rng_uniform_int_MT64(rng_64, parent_subpop_size_));
}
```

**低适应度个体的过滤原理**：

设总体 N=4，适应度 $W = [0.9, 1.0, 0.5, 0.1]$，$\sum W = 2.5$。

则采样概率分别为 $[0.36, 0.40, 0.20, 0.04]$。

适应度 0.1 的个体**并未被强制删除**，而是其被选为亲本的概率被压低至 4%。经过多代，其携带的有害等位基因的期望频率会持续下降，直至被遗传漂变消除（或在隐性情况下漂浮于低频）。

### 步骤 9：减数分裂 & 变异传递

被选中的亲本经过重组（断点已由随机数抽取）产生后代单倍体组，过程中父本的各突变随机分配到下一代。

---

## 五、第四阶段：突变频率演化——q_j 如何被过滤

### 步骤 10：引用计数追踪（频率的实现）

每次完成一代生殖后，`Population::MaintainMutationRegistry()` 对每个突变做全局引用计数统计：

```c
// population.cpp，TallyMutationReferencesAcrossPopulation()
// 对所有子代单倍体组扫描，每找到一个突变就 ++refcount
slim_refcount_t refcount = gSLiM_Mutation_Refcounts[mut_index];

// 等位基因频率：
double q_j = refcount / total_haplosome_count;  // 总单倍体组数（二倍体 × 2）
```

### 步骤 11：突变的命运判断

```c
// population.cpp，RemoveAllFixedMutations()
if (refcount == 0) {
    mutation->state_ = MutationState::kLostAndRemoved;    // 丢失（自然选择或漂变消除）
    remove_mutation = true;
}
else if (refcount == total_haplosome_count) {
    if (mutation_type->convert_to_substitution_) {
        // WF 模型默认 true → 固定的突变转为 Substitution 对象，从 registry 删除
        mutation->state_ = MutationState::kFixedAndSubstituted;
    }
    // nonWF 模型默认 false → 留在 registry 中（因为可能再次变为多态）
}
// 否则 0 < refcount < total → 仍为 Segregating（分离中），保持 kInRegistry 状态
```

**突变状态与 calcInbreedingLoad 中的 q 值对应关系**：

| 突变状态 | 在 calcInbreedingLoad 中的 q 值 |
|---------|-------------------------------|
| `kInRegistry`（分离中） | `refcount / total_haplosome_count` |
| `kLostAndRemoved`（丢失）| 0.0，已不在 subsetMutations() 中 |
| `kFixedAndSubstituted`（固定）| 1.0，不影响 B（B 只计算分离的有害突变）|

---

## 六、第五阶段：calcInbreedingLoad(B) 的精确计算步骤

到达用户调用 `calcInbreedingLoad()` 的时刻，上述所有过滤已作用了若干代，等位基因频率 $q_j$ 已经是自然选择 + 遗传漂变共同塑造的结果。

### 步骤 12：获取有害突变集合

```eidos
// slim_functions.cpp，gSLiMSourceCode_calcInbreedingLoad
muts = species.subsetMutations(chromosome=chromosome);  // 取 kInRegistry 中的突变
muts = muts[muts.selectionCoeff < 0.0];                 // 只保留有害突变（s < 0）
```

`subsetMutations()` 只返回当前仍在 `mutation_registry_` 中（`kInRegistry` 状态）的突变，即已被选择/漂变过滤后仍分离的突变。

### 步骤 13：获取等位基因频率

```eidos
q = haplosomes.mutationFrequenciesInHaplosomes(muts);
// 底层：refcount / tallied_haplosome_count（已缓存）

inHaplosomes = (q > 0);          // 过滤掉频率为 0 的（已丢失但尚未从列表移除的）
muts = muts[inHaplosomes];
q = q[inHaplosomes];
```

### 步骤 14：获取 s 和 h

```eidos
s = -muts.selectionCoeff;    // 注意：SLiM 的 selectionCoeff 为负，Morton 1956 中 s > 0
s[s > 1.0] = 1.0;            // 截断：s 不能超过 1（不能比致死还致死）
h = muts.mutationType.dominanceCoeff;
```

### 步骤 15：代入 Morton et al. 1956 公式

```eidos
B = (sum(q*s) - sum(q^2*s) - 2*sum(q*(1-q)*s*h))
```

展开含义：

$$B = \underbrace{\sum_j q_j s_j}_{\text{项①}} - \underbrace{\sum_j q_j^2 s_j}_{\text{项②}} - \underbrace{2\sum_j q_j(1-q_j) s_j h_j}_{\text{项③}}$$

**Morton 1956 公式各项的遗传学含义**：

| 项 | 含义 |
|---|---|
| $\sum q_j s_j$ | 若全部近交（F=1），每个有害等位基因都会暴露为纯合，期望适应度降低的上界 |
| $-\sum q_j^2 s_j$ | 减去已经是纯合的部分（无需近交就已表达了的） |
| $-2\sum q_j(1-q_j) s_j h_j$ | 减去杂合时已经通过显性系数部分表达的适应度代价 |

$B$ 代表**每个单倍体组上的致死当量**（lethal equivalents per haploid genome），反映了"若使一个个体完全近交（F=1，即使所有基因座纯合）时，其相对于随机交配后代额外增加的遗传死亡风险"。

---

## 七、低适应度个体的过滤参数作用原理全图

```
初始突变产生
     ↓
  s ~ DFE（如 Gamma 分布，均值 s̄ < 0）
  h = 显性系数（如 h=0.25 → 隐性有害）
     ↓
缓存：cached_one_plus_dom_sel_ = max(0, 1 + h·s)   ← 杂合适应度
     cached_one_plus_sel_     = max(0, 1 + s)       ← 纯合适应度
     ↓
每代 UpdateFitness()：
  个体适应度 W_i = ∏ over all mutations:
    杂合出现 → × (1 + h·s)    ← 部分降低
    纯合出现 → × (1 + s)      ← 完全降低
     ↓
UpdateWFFitnessBuffers()：
  cached_parental_fitness_[i] = W_i（个体"被采样权重"）
  gsl_ran_discrete_preproc() → Vose-Alias 查找表
     ↓
EvolveSubpopulation()：
  parent = gsl_ran_discrete(rng, lookup_parent_)
  采样概率 ∝ W_i
  ──────────────────────────────────────────────
  W_i 低 → 被选为亲本的概率低 → 其携带的有害等位基因频率 q_j 下降
  W_i = 0 → 永远不会被选为亲本（致死突变的纯合子彻底淘汰）
  ──────────────────────────────────────────────
     ↓
MaintainMutationRegistry()：
  按 refcount 统计 q_j = count / total_haplosome_count
  q_j = 0 → kLostAndRemoved（从 registry 删除）
     ↓
calcInbreedingLoad(haplosomes)：
  muts ← subsetMutations()（只有 kInRegistry 的分离突变）
  muts ← muts[s < 0]（过滤出有害突变）
  q ← mutationFrequenciesInHaplosomes(muts)
  s ← -muts.selectionCoeff
  h ← muts.mutationType.dominanceCoeff
  B = sum(q*s) - sum(q²*s) - 2*sum(q*(1-q)*s*h)
```

---

## 八、三个关键过滤机制总结

### 机制 1：适应度乘积截断（硬边界）

`max(0, 1+s)` 保证适应度非负。若某突变 $s \leq -1$（致死），则含该突变的纯合个体 $W=0$，永远不被选为亲本。其对应等位基因会在一代内从纯合子消失（但若 $h < 1$，杂合子仍可携带）。

### 机制 2：加权采样压制（核心过滤参数）

适应度 $W_i$ 直接决定 `cached_parental_fitness_[i]` 中的权重值。$W_i$ 越低，个体在每代被"过滤掉"（未被选为亲本）的概率越高。随代数积累，含高负荷有害突变个体的后代比例持续收缩，导致有害等位基因频率 $q_j$ 的期望值随选择压力下降。

| 参数 | 对过滤强度的影响 |
|------|----------------|
| $s$（选择系数绝对值）大 | 纯合时适应度降低更多，淘汰更快 |
| $h$（显性系数）大 | 杂合时也有明显适应度代价，隐性突变不易被过滤 |
| $N$（有效种群大小）小 | 遗传漂变作用强，有害突变可能随机固定（遗传漂变与选择的竞争） |
| 迁移率高 | 不同亚群间有害突变频率趋于平衡 |

### 机制 3：引用计数归零清除（最终淘汰）

当某突变在全群所有单倍体组中引用计数 `refcount=0`，该突变对象状态转为 `kLostAndRemoved`，从注册表中永久删除，不再参与后续计算（包括 `calcInbreedingLoad`）。

---

## 九、最终 B 值的解读

$B$ 的数值反映的是：经历了以上所有过滤轮次后，仍然**以分离状态存活**在群体中的有害突变的频率加权总负荷。它们是漂变、选择和显性系数共同博弈后的"幸存者"，其累积效应即为群体的近交负荷。

- **B 偏高**：群体中积累了较多的有害隐性突变（低 h、低 q），近交会使它们暴露为纯合，大幅降低后代适应度
- **B 偏低**：群体经历了有效的净化选择，有害突变已被大量清除，或群体本身具有高近交历史（隐性有害突变已在历史近交中被清除）
- **B 受 h 影响显著**：若 h 接近 0（完全隐性），项③趋近于 0，B 值主要由低频有害等位基因的频率和选择系数决定
