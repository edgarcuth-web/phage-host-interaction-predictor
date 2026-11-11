# 两套 rescue 估计量的对账

`docs/annotation-recall-audit.md`（PR #31）和
`docs/retrieval-rescue-benchmark-v0.md` + `docs/offline-retrieval-rescue-evaluator.md`
（PR #29、#30）都在谈"证据路线漏掉、另一条路线捞回"的阳性。
**它们定义的不是同一个量。** 本文件说明差别、各自的适用范围，以及哪些数字可比、哪些不可比。

本文件不改变任何已报告的数值，也不否定任何一方。

## 两者各做什么

| | benchmark v0 + evaluator | recall audit |
|---|---|---|
| 性质 | 协议规范 + 离线评估器实现 | 在已发表实测数据上的测量 |
| 是否训练模型 | 不训练（路线 G 标注为 *not yet trained*） | 不训练 |
| 回答的问题 | 怎么把 learned 与 evidence 的等预算比较做对 | evidence 路线本身漏掉多少真阳性 |
| 数据 | 冻结的 177 host × 94 phage，16,597 条审计 MLC 标签 | 公开的 402 × 94 序数矩阵；PhageHostLearn 215 × 132 |

benchmark v0 自己留了一处空缺，原文为：

> Do **not** label the annotation model "mechanism-only retrieval"… A genuine
> mechanistic-rule baseline needs a separately documented, source-supported
> receptor/RBP/defense compatibility policy and **measured coverage of those inputs**.

recall audit 的第四节和第十节就是那个 measured coverage：
前者给出 Escherichia 上单特征注释规则的清单（96 个 phage 中 75 个存在规则，
ABC_serotype 供给 32 个、LPS_type 22 个）和规则迁移的召回损失；
后者给出 Klebsiella 上 capsule 匹配的覆盖天花板
（precision 上限 0.632 / recall 0.282，68.1% 的阳性没有 K-locus 同一性 ≥0.90 的已感染邻居）。

所以两者互补：一方定义怎么跑比较，另一方量化 baseline A 配不配叫机制筛选。

## 四处不可比

### 一、预算单位：per-host 的 K 对全局池化的 N

benchmark v0 的预算是 per-host：每个 host 给 K 个槽位（primary K=2），
对应"一个病人分离株、做 K 次点板"。

recall audit 第十节的 rescue 用的是**全局池化 top-N**——在整个矩阵上取评分最高的
183 格（N 取 K-locus 规则在同一性 ≥0.90 处选中的格数），一株菌可以占用多个槽位。
那回答的是"整个库里该挑哪些格子做实验"，是库策展的问法。

**因此第十节的 2/348 不能当作 evaluator 的 `primary_k_rescue` 值。**
换成 per-host K=2，槽位总数更少，数字会不同。

recall audit 第五节的 recall@k 本身是 per-strain 的，与 benchmark 的排序方向一致
（给定一株菌，对 phage 排序）。只有第十节的 rescue 口径不同。

### 二、并列策略

benchmark v0 要求 membership 用声明过的确定性策略
（`score descending then phage_id ascending`），指标用并列排列的精确期望，
并明确要求两者不可混用。

recall audit 的脚本原先用 `np.argsort(-score)`，并列时由输入列顺序决定——未声明。
本 PR 已把 `step2c_recall_at_k.py`、`step3b_recall_robust.py`、`step5b_rescue.py`
改为与 benchmark 一致的策略（分数降序，并列时按标识符升序；
第三个脚本对全局格还额外按菌株标识符再按 phage 标识符升序）。

重跑后的实际影响：

| 受影响的数字 | 改前 | 改后 |
|---|---|---|
| prevalence 占 oracle，k=2 | 0.428 | 0.429 |
| nn_phylo 占 oracle，k=1 | 0.568 | 0.569 |
| Klebsiella rescue 计数 | 93 / 18 / 2 / 280 | 不变 |

只有并列最多的两条基线动了 0.1 个百分点，没有结论改变。
Klebsiella 不变的原因已单独验证：183 格的截断值为 0.901327，
严格大于它的有 182 格，截断处的并列组只有 1 格。
更高处确有大并列组（32 格和 23 格恰好等于 1.0000），
阈值若落在那里则membership会是任意的，所以这个策略必须声明。

**recall audit 的指标仍未做并列排列的精确期望**，只用了确定性 membership 策略。
所以它的 recall@k 与 evaluator 的 tie-expected Hit@K / Recall@K 是两个不同的估计量。

### 三、不完整 host 的处理

evaluator 的行为是 `Only complete-label hosts are scored`，
unknown host 保留在排序审计里但不计 rescue credit。

recall audit 按观测格计算，未排除测试覆盖不完整的菌株。
在 Escherichia 上影响很小（38,592 格中仅 157 格缺失，0.41%）。
**在 Klebsiella 上影响不小**：28,380 格中只有 12,970 格被测过（45.7%），
且缺失是块状的，每株已测格数从 46 到 86。
因此第十节的 per-strain 分母随各株的测试覆盖度变化，这是一处偏倚。

evaluator 的处理更严。若要把第十节的数字纳入 benchmark 口径，
需要先按完整标签面板筛选菌株。

### 四、标签端点

| | 端点 | 规模 |
|---|---|---|
| benchmark v0 | 审计过的 MLC > 0 / non-positive，41 对因缺失主读数被排除 | 177 host × 94 phage，16,597 条 |
| recall audit 第一至九节 | 公开矩阵的 0–4 序数分数，在 ≥1 处二值化（≥2/≥3/≥4 做敏感性） | 402 × 94，38,592 格 |
| recall audit 第十节 | PhageHostLearn 的 1:10 稀释点板二值标签 | 215 × 132，12,970 格已测 |

三者是三个端点。**不要把任一方的 Hit@K 或 Recall@K 直接与另一方对比。**

## recall audit 对 benchmark 的两点建议

### 一、报告 oracle 归一化行

benchmark v0 把 Hit@2 定为 primary estimand，Recall@K 作为次级指标，
但没有按每个 host 自身的阳性数归一化 Recall@K。

每个 host 的 recall@K 受该 host 阳性数封顶。在 Gaborieau 的完整矩阵上
（每株阳性数中位 16.5），完美排序的 oracle 为：

```
oracle recall@k:  k=1 0.115   k=2 0.205   k=3 0.283   k=5 0.419   k=10 0.630
```

所以一个 0.13 的 Recall@2 其实是可达上限的 63%，不是失败。
建议 evaluator 的输出里加一行 oracle recall@K
（每个 host 取 `min(K, n_positives) / n_positives`，再按 host 宏平均），
以及各方法占 oracle 的比例。没有这一行，Recall@K 会被系统性低估解读。

### 二、第十节不是路线 G 的一个实例

benchmark v0 写明：

> Untrained Evo2 vector cosine similarity is not a validated compatibility score.

recall audit 第十节用的正是未训练的中心化 ESM-2 cosine。
它的结论——**同一批 K-locus 序列的学习表示，作为相似度使用时 rescue 为 2/348，
且按 phage 衡量显著劣于原始 capsule 同一性（配对 Wilcoxon p=0.036）**——
界定的是"表示相似度单独能做到什么"的上界。

**按 benchmark 的分类法，那不是路线 G。** 路线 G 要求训练过的、
使用 phage 序列信息的兼容性模型；PhageHostLearn 实际拟合的是
RBP × locus 拼接后的 XGBoost，严格强于最近邻 cosine。
所以第十节既不实现也不否证路线 G，不应被读作"G 已经被测过了"。

第十节真正对 benchmark 有用的部分是另一条：宿主侧与噬菌体侧的覆盖天花板相近
（汇总 AUROC 0.750 对 0.716）且不互补，提示剩余信号在交互本身，
而不在任一边的边缘分布。这支持 benchmark 坚持 G 必须使用 phage 序列信息的要求。

## 两份文件都成立的共同结论

- 未观测对不是阴性。两边都不把 unknown 转成负标签。
- evidence 路线不是机制预言机。benchmark 在命名上拒绝这么叫，
  recall audit 量化了它漏多少。
- "learned-only positive" 不是发现了未知的因果基因。两边都明确写了这一条。
- 需要 prevalence 基线。benchmark 列为路线 P，recall audit 实测它在
  Escherichia 上已达 oracle 的 0.43–0.50，在 Vibrio 上只有 0.038。

## 下一步（本文件不实施）

1. 在 evaluator 的输出里加 oracle recall@K 行与占 oracle 比例。
2. 若要把 Klebsiella 纳入 benchmark 口径：先按完整标签面板筛菌株，
   再把 rescue 从全局池化 N 改为 per-host K。
3. recall audit 若要与 evaluator 的指标直接可比，需要实现并列排列的精确期望，
   而不只是确定性 membership 策略。
4. 任一方的 Hit@K / Recall@K 在跨端点引用时必须标注端点。
