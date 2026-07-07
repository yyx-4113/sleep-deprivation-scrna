# 福建省自然科学基金面上项目 — 实验标准操作流程 (SOP)

> 项目：急性睡眠剥夺下丘脑Pomc神经元转录调控机制及Pomc-Pcsk2正反馈环路的功能研究
> 版本：v1.0 | 2026-05-22
> 适用：内容一（表达验证）、内容二（因果验证）、内容三（机制初探）

---

## SOP-01: 急性睡眠剥夺小鼠模型建立

### 1.1 动物信息
- **品系**: C57BL/6J 雄性小鼠
- **周龄**: 8周龄（体重20-25g）
- **来源**: 福建医科大学实验动物中心 / 斯莱克实验动物
- **饲养条件**: 12h:12h 明暗循环（ZT0=开灯 07:00），22±1°C，50±10%湿度，自由饮食饮水
- **适应期**: 至少7天

### 1.2 分组设计
| 组别 | 缩写 | 处理 | 时间窗 | n |
|------|------|------|--------|---|
| 正常睡眠对照组 | A1 | 不干预，于ZT0-5间自然睡眠 | ZT0-5 | 6 |
| 睡眠剥夺组 | A2 | 轻柔触碰法剥夺睡眠5h | ZT0-5 | 6 |
| 恢复睡眠组 | A3 | SD 5h后允许自由睡眠2h | ZT0-7 | 6 |

总计: **18只小鼠**

### 1.3 睡眠剥夺操作流程
1. 于 ZT0（07:00，开灯时）开始
2. 实验者坐在鼠笼旁，观察小鼠行为
3. 当小鼠出现以下入睡迹象时，用软毛刷轻触其背部或轻轻摇晃笼子：
   - 静止不动 >10秒
   - 闭眼
   - 蜷缩姿势
4. **关键**: 触碰力度仅需唤醒小鼠，避免引起应激反应
5. 每笼同时剥夺 2-3 只小鼠（不要单独剥夺以降低社交隔离应激）
6. 连续进行 5 小时
7. **恢复组**: SD结束后让小鼠自由睡眠2h，取材

### 1.4 行为验证
- 全程录像，每30分钟记录一次运动活性
- SD期间应观察到持续的活动增加
- 记录被触碰次数作为睡眠压力指标（>30次/小时为有效SD）

---

## SOP-02: 组织取材与样本处理

### 2.1 设备与试剂
- 异氟烷麻醉机
- 冰台（预冷至-20°C）
- 手术器械（剪刀、镊子、止血钳，180°C烘烤4h灭活RNase）
- RNAlater / Trizol
- 1.5mL RNase-free EP管
- 干冰、液氮

### 2.2 取材流程

**血液收集**:
1. 异氟烷深度麻醉小鼠
2. 摘眼球取血：0.5-1.0 mL 全血 → EDTA抗凝管
3. 4°C, 3000×g 离心 15 分钟 → 取上清血浆
4. 分装50μL/管 → -80°C保存（用于ELISA）

**脑组织取材**:
1. 断头处死后快速剥离全脑，置于预冷的人工脑脊液（aCSF）中
2. 在冰台上使用解剖显微镜分离：
   - **下丘脑**: 视交叉至乳头体之间，深至第三脑室底，重点含弓状核区域
   - **脑干**: 第四脑室底区域
   - **前额叶皮层**: Bregma +3.0 至 +1.5 mm
3. 每个脑区一分为二：
   - 一半 → Trizol（RNA提取）
   - 一半 → RIPA裂解液+蛋白酶抑制剂（蛋白提取）
4. 下丘脑额外取2只/组 → 4% PFA固定（免疫荧光）
5. 液氮速冻后 → -80°C保存

---

## SOP-03: qRT-PCR（实时定量PCR）

### 3.1 引物序列

| 基因 | 正向引物 (5'-3') | 反向引物 (5'-3') | 产物(bp) |
|------|-----------------|-----------------|----------|
| Pomc | ACCTCACCACGGAGAGCA | GCGAGAGGTCGAGTTTGC | 120 |
| Pcsk2 | TGGCTGACACAGACCTACG | CCACTCGGTGTAGTTCTCCA | 135 |
| Scg2 | GAGCTGAGCGATGACTCAGA | TCCAAGTGTTGAGCCAGTCC | 108 |
| Cga | CTGGTCCCAGATGTGTCACA | TGTGTAGCCAAAGCGGAAG | 142 |
| Gnas | CTGCCAACAACTACGACTCG | TTGATGCCCTCTGACACCTT | 115 |
| Dbp | CCGTGGAGGTGCTAATGACC | CTCTTCATGGCGTAGAATGCG | 128 |
| Nr3c1 | AGCAGTGTGAAACCTGGAAGC | GATTTCAAAGGCAGTCTGGTGC | 98 |
| Malat1 | GGAGTCCGAGGACGAATGC | ATTGGCACGCACATTTCGC | 86 |
| Ccnd2 | GAGTGGGAACCTGGTAGTGTTG | CGCACTTCAATCTCCTCAGCA | 145 |
| Fos | CGGGTTTCAACGCCGACTA | TTGGCACTAGAGACGGACAGA | 105 |
| Jun | AACGACCTTCTACGACGATGC | TGAGTTGGCACCCACTGTT | 118 |
| Egr1 | TCGGCTCCTTTCCTCACTCA | CTCCAGCTTAGGGTAGTGTAGAG | 156 |
| Gapdh | AGGTCGGTGTGAACGGATTTG | TGTAGACCATGTAGTTGAGGTCA | 123 |
| Actb | GGCTGTATTCCCCTCCATCG | CCAGTTGGTAACAATGCCATGT | 154 |

### 3.2 操作流程
1. **RNA提取** (Trizol法):
   - 组织 + 500μL Trizol → 电动匀浆
   - 加100μL氯仿 → 剧烈振荡15s → 室温3min
   - 4°C, 12000×g, 15min → 取上层水相
   - 加250μL异丙醇 → 冰上10min
   - 4°C, 12000×g, 10min → 弃上清
   - 75%乙醇洗2次 → 晾干 → DEPC水溶解
   - Nanodrop测浓度: A260/A280 = 1.8-2.1

2. **反转录** (Takara PrimeScript RT kit):
   - 500ng total RNA → 10μL体系
   - 37°C 15min → 85°C 5s → 4°C
   - cDNA稀释5倍备用

3. **qPCR** (TB Green Premix Ex Taq II):
   - 10μL体系: 5μL TB Green + 0.4μL each primer (10μM) + 1μL cDNA + 3.2μL H₂O
   - 程序: 95°C 30s → [95°C 5s → 60°C 30s] × 40 cycles → 熔解曲线
   - 每样3个技术重复
   - 相对定量: ΔΔCt法，Gapdh和Actb双内参

### 3.3 预期结果
- SD组: Pomc ↓80%+, Pcsk2 ↓30-50%, Scg2 ↓, Dbp ↓
- 恢复组: Pomc反弹至基线以上
- Fos/Jun/Egr1: SD组可能上调（IEGs的急性应激响应）

---

## SOP-04: Western Blot

### 4.1 抗体信息

| 靶蛋白 | 货号/供应商 | 分子量 | 稀释比 |
|--------|-----------|--------|--------|
| POMC | Abcam ab227487 / CST #7074 | 31 kDa (前体) | 1:1000 |
| PC2 (PCSK2) | Abcam ab28540 | 70 kDa | 1:500 |
| SCG2 | Proteintech 12264-1-AP | 70 kDa | 1:1000 |
| DBP | Abcam ab181029 | 40 kDa | 1:1000 |
| β-actin | CST #4970 | 45 kDa | 1:5000 |

### 4.2 操作流程
1. **蛋白提取**:
   - RIPA + PMSF(1mM) + 蛋白酶抑制剂cocktail (1:100)
   - 组织 + 200μL裂解液 → 电动匀浆 → 冰上30min, 每10min涡旋
   - 4°C, 14000×g, 15min → 取上清
   - BCA法定量 → 调至等浓度

2. **电泳**:
   - 10-12% SDS-PAGE（POMC用15%以分离小分子）
   - 上样量: 30-50μg/孔
   - 80V 30min（积层胶）→ 120V 60-90min（分离胶）

3. **转膜**:
   - PVDF膜（甲醇激活15s）
   - 湿转: 300mA, 60-90min（根据分子量调整）
   - 丽春红染色确认转膜效率

4. **封闭与抗体孵育**:
   - 5%脱脂牛奶/TBST, 室温1h
   - 一抗4°C过夜
   - TBST洗 3×10min
   - HRP二抗 (1:5000-1:10000), 室温1h
   - TBST洗 3×10min

5. **显影**:
   - ECL化学发光液（A液:B液=1:1）
   - 化学发光成像仪曝光（自动/梯度曝光）
   - ImageJ定量灰度值

### 4.3 预期结果
- SD组 vs 对照: POMC ↓ (50-70%), PC2 ↓ (30-50%), SCG2 ↓

---

## SOP-05: 免疫荧光双标

### 5.1 试剂与抗体

| 试剂/抗体 | 供应商 | 货号 | 稀释比 |
|-----------|--------|------|--------|
| POMC 一抗 (兔源) | CST | #7074 | 1:200 |
| PC2 一抗 (小鼠源) | Santa Cruz | sc-376721 | 1:100 |
| 驴抗兔 Alexa Fluor 488 | Invitrogen | A21206 | 1:500 |
| 驴抗小鼠 Alexa Fluor 594 | Invitrogen | A21203 | 1:500 |
| DAPI 封片剂 | Vector Labs | H-1200 | — |

### 5.2 操作流程
1. **冰冻切片制备**:
   - 4% PFA灌注固定 → 取脑 → 4% PFA后固定4h → 30%蔗糖脱水（4°C过夜至沉底）
   - OCT包埋 → -20°C恒冷切片机, 14-16μm冠状切片
   - 贴于多聚赖氨酸包被载玻片 → -20°C保存

2. **免疫荧光染色**:
   - 室温复温20min → PBS洗 3×5min
   - 抗原修复: 柠檬酸钠缓冲液(pH 6.0), 95°C, 10min → 冷却至室温
   - 封闭: 5%驴血清 + 0.3% Triton X-100/PBS, 室温1h
   - 一抗混合液（POMC + PC2）4°C过夜
   - PBS洗 3×10min
   - 二抗混合液（避光）, 室温1h
   - PBS洗 3×10min（避光）
   - DAPI封片剂封片 → 4°C避光保存

3. **成像与分析**:
   - 共聚焦显微镜（或宽场荧光显微镜）
   - 20× 和 40× 拍摄弓状核区域（Bregma -1.4 至 -2.0 mm）
   - ImageJ/Fiji 定量:
     - 每个神经元 POMC/PC2 荧光强度
     - Pearson 相关系数（共定位分析，Coloc 2插件）
     - 每只小鼠统计 ≥50个POMC+神经元

### 5.3 预期结果
- SD组: POMC+PC2双阳性细胞数↓，共定位系数↓
- 恢复组: 部分恢复

---

## SOP-06: ELISA（酶联免疫吸附实验）

### 6.1 试剂盒

| 检测物 | 供应商 | 货号 | 检测范围 | 样本量 |
|--------|--------|------|----------|--------|
| α-MSH | Phoenix Pharmaceuticals | EK-043-01 | 0-25 ng/mL | 50 μL |
| ACTH | MD Bioproducts | M046006 | 0-1000 pg/mL | 25 μL |
| β-endorphin | Phoenix Pharmaceuticals | EK-022-14 | 0-100 ng/mL | 50 μL |

### 6.2 操作流程
1. 血浆样本从-80°C取出 → 冰上融解
2. 按试剂盒说明书配制标准品梯度
3. **α-MSH/β-endorphin**（竞争法ELISA）:
   - 50μL样品/标准品 + 25μL一抗 + 25μL生物素化肽 → 室温2h
   - 洗板5次 → 100μL SA-HRP → 室温1h
   - 洗板5次 → 100μL TMB底物 → 室温20-30min（避光）
   - 100μL 2N HCl终止 → 450nm读数
4. **ACTH**（夹心法ELISA）:
   - 25μL样品/标准品 + 100μL检测抗体 → 室温2h
   - 后续同上
5. 4参数logistic曲线拟合 → 计算浓度
6. 每样2个技术重复

### 6.3 注意事项
- 肽类激素易降解：全程冰上操作，加抑肽酶（Aprotinin 500 KIU/mL）保护
- 如信号低于检测限 → 考虑冻干浓缩或改用RIA（放射免疫法，更高灵敏度）

---

## SOP-07: 细胞培养与siRNA敲低

### 7.1 细胞系

**首选**: **GT1-7细胞**（小鼠下丘脑GnRH神经元细胞系，内源性表达Pomc）
- 来源: ATCC / MilliporeSigma (SCC296)
- 培养基: DMEM高糖 + 10% FBS + 1% Pen/Strep
- 培养条件: 37°C, 5% CO₂
- 传代: 1:4-1:6, 每3-4天

**备选**: **N2a细胞**（小鼠神经母细胞瘤）
- 来源: ATCC CCL-131
- 培养基: DMEM + 10% FBS + 1% Pen/Strep
- 需转染Pomc表达质粒后使用

### 7.2 siRNA设计

设计3条Pomc siRNA + 1条阴性对照 (scrambled):

| siRNA | 靶序列 (sense, 5'-3') | 靶区域 |
|-------|----------------------|--------|
| siPomc-1 | GCUUCUGACUCCUUCAUAA | CDS |
| siPomc-2 | GGACAUCUUUAGCCUCAGA | CDS |
| siPomc-3 | CCUCUGAGAAGUCUACAUA | CDS |
| siNC | UUCUCCGAACGUGUCACGU | (无靶标) |

- 合成: 吉玛基因 / GenePharma, 2'-OH modified
- 使用浓度: 50nM

### 7.3 转染流程

1. **Day 0**: 种板: GT1-7细胞, 2×10⁵/孔（6孔板）
2. **Day 1**（细胞汇合度60-70%）:
   - 管A: 5μL siRNA (20μM) + 250μL Opti-MEM
   - 管B: 5μL Lipofectamine RNAiMAX + 250μL Opti-MEM
   - 混合A+B → 室温15min → 滴加入孔
3. **Day 2**（转染后24h）: 换新鲜完全培养基
4. **Day 3**（转染后48h）: 收细胞
   - 一半 → Trizol（qPCR验证敲低效率）
   - 一半 → RIPA（WB验证蛋白水平）

### 7.4 敲低效率验证

- qPCR: Pomc mRNA 应在siPomc组降低 >70%（与siNC比较）
- WB: POMC蛋白应显著降低
- 选择敲低效率最高的siRNA进行后续RNA-seq和Rescue实验

### 7.5 靶基因qPCR Panel（20个基因）

在siPomc vs siNC细胞中检测:
- **核心环路**: Pcsk2, Scg2, Cga, Gnas
- **预测靶基因**: Malat1, Ccnd2, Dbp, Nr3c1
- **共表达TFs**: Fos, Jun, Junb, Egr1, Crem, Atf3
- **其他预测**: Rbm3, Rnaset2a, Tsc22d3, Mt1
- **内参**: Gapdh, Actb

### 7.6 Rescue实验

1. siPomc转染后24h → 转染POMC过表达质粒（pcDNA3.1-Pomc, 0.5μg/孔）
2. 对照组: siPomc + pcDNA3.1空载体
3. 转染后48h → 收细胞 → qPCR检测靶基因
4. 预期: POMC过表达后靶基因（特别是Pcsk2）表达恢复

---

## SOP-08: 3' mRNA-seq

### 8.1 实验设计

| 组别 | n | 总计 |
|------|---|------|
| siPomc | 3 | 3个生物学重复 |
| siNC | 3 | 3个生物学重复 |
| **合计** | **6** | **6个样本** |

### 8.2 文库构建与测序

1. 总RNA提取（Trizol法，同上）→ QC: RIN > 7.0（Bioanalyzer）
2. **QuantSeq 3' mRNA-seq**（Lexogen）:
   - 500ng total RNA起始
   -  oligo-dT priming → 只捕获polyA+ mRNA的3'端
   - 优点: 不需要polyA富集和mRNA打断，节省成本
3. 测序: NovaSeq 6000 / HiSeq X, PE150, 每个样本 5M reads
4. 送样: 福建伯豪 / 华大基因（报价约2500元/样本）

### 8.3 数据分析

- 比对: STAR → mm39参考基因组
- 定量: featureCounts
- 差异分析: DESeq2 (adjusted p < 0.05, |log2FC| > 0.5)
- GSEA: clusterProfiler（GO + KEGG）
- 与虚拟KO预测结果对比: 计算Spearman相关系数，绘制散点图

---

## SOP-09: ChIP-qPCR

### 9.1 抗体

| 抗体 | 供应商 | 货号 | 用量/反应 |
|------|--------|------|-----------|
| CREB (48H2) | CST | #9197 | 2-5 μg |
| c-Fos (9F6) | CST | #2250 | 2-5 μg |
| JunB (C37F9) | CST | #3753 | 2-5 μg |
| EGR1 (15F7) | CST | #4153 | 2-5 μg |
| Normal Rabbit IgG | CST | #2729 | 2-5 μg (阴性对照) |

### 9.2 靶区域引物（Pcsk2启动子）

基于JASPAR预测结果设计的候选ChIP-qPCR区域：

| 位点 | 相对TSS | 预测结合的TF | 正向引物 (5'-3') | 反向引物 (5'-3') | 产物(bp) |
|------|---------|-------------|-----------------|-----------------|----------|
| Site 1 (TRE) | -607 | c-Fos/JunB | GACTCAGAGCTGGGATGGAG | CCTAGCCCTGACTTGACTGC | 112 |
| Site 2 (CRE) | -291 | CREB1/ATF3 | TGACTCTGAGGCAGAGAGGA | TCTGGACAGGCCTTTGTTCA | 98 |
| Site 3 (Egr) | -144 | EGR1 | GAGGCTGGTGAGATGGTTCA | CCAGAGCTAGGCGAGAGAAA | 105 |
| Negative Ctrl | -3500 | (阴性区域) | GTCTGCTGAGCTGGTTACCC | CAGTGGCACAGTGACAGAGA | 125 |

### 9.3 操作流程

1. **交联**:
   - 下丘脑组织（每组n=3, ~15-20mg）→ 剪碎至1mm³
   - 1%甲醛/PBS, 室温交联10min（旋转）
   - 甘氨酸（终浓度0.125M）终止交联 → 室温5min
   - 冷PBS洗2次

2. **染色质制备**:
   - 裂解液（含PMSF+蛋白酶抑制剂）匀浆 → 冰上10min
   - 4°C, 3000×g, 5min → 弃上清，保留核沉淀
   - 核裂解液重悬 → 超声打断（Bioruptor / Covaris）:
     - 条件: 30s ON / 30s OFF, 15-20 cycles → 200-500bp片段
   - 取10μL检测片段大小（琼脂糖凝胶电泳）
   - 4°C, 14000×g, 10min → 取上清

3. **免疫沉淀**:
   - 取5-10μg染色质（按DNA量算）→ 稀释至500μL ChIP dilution buffer
   - 取1% 作为 "Input"
   - 加2-5μg抗体（目标抗体 / IgG阴性对照）→ 4°C旋转过夜
   - 加30μL Protein G磁珠 → 4°C旋转2h

4. **洗涤与洗脱**:
   - 依次用低盐→高盐→LiCl→TE buffer各洗1次, 每次5min
   - 洗脱液: 200μL (1% SDS + 0.1M NaHCO₃) → 室温15min（旋转）
   - 收集上清, 重复洗脱1次（总计400μL）

5. **解交联与DNA纯化**:
   - 加NaCl至终浓度0.2M → 65°C 4h或过夜
   - 加EDTA+Tris+蛋白酶K → 45°C 1h
   - 酚/氯仿抽提 → 乙醇沉淀 → TE溶解
   - QIAquick PCR Purification Kit纯化 → 50μL洗脱

6. **qPCR检测**:
   - 同SOP-03 qPCR条件
   - 计算 % Input = 2^(Ct_Input - Ct_IP) × 100
   - 相对富集 = % Input (抗体) / % Input (IgG)

### 9.4 预期结果
- SD组 vs 对照组:
  - CREB在Pcsk2启动子Site 2 (-291) 富集减少
  - c-Fos/JunB在Site 1 (-607 TRE) 富集减少
- IgG对照组应无明显富集（fold enrichment < 2）

---

## SOP-10: 双荧光素酶报告基因实验（可选）

### 10.1 质粒构建

1. **报告载体**: pGL4.10[luc2] (Promega)
   - 克隆Pcsk2启动子片段 (-2000~+100bp) 至luciferase上游
   - 构建截短体: -607 TRE位点单独（用于验证AP-1调控）

2. **TF表达质粒**:
   - pcDNA3.1-Fos, pcDNA3.1-Junb, pcDNA3.1-Creb1
   - 加C-terminal HA/Flag tag便于验证表达

3. **内参**: pRL-TK (Renilla luciferase)

### 10.2 操作流程

1. HEK293T细胞 → 96孔板, 1×10⁴/孔
2. 转染（脂质体法Lipofectamine 3000）:
   - 报告质粒(Pcsk2-promoter-luc): 100ng/孔
   - TF表达质粒: 50ng/孔
   - 内参(pRL-TK): 10ng/孔
   - 每组3孔重复
3. 转染后48h → Dual-Luciferase Reporter Assay
4. 计算: Firefly/Renilla 比值 → 归一化后比较

### 10.3 预期结果
- Fos+Junb共转染 → Pcsk2启动子活性↑ (vs 空载体)
- Creb1转染 → Pcsk2启动子活性↑
- 截短TRE位点 → AP-1共转染诱导效应消失

---

## 附录A: 关键试剂采购清单

| 类别 | 试剂/耗材 | 供应商 | 预估价格(元) |
|------|----------|--------|-------------|
| 动物 | C57BL/6J小鼠 (♂, 8w) | 斯莱克/维通利华 | 120/只 |
| RNA | Trizol reagent | Invitrogen | 800/100mL |
| RNA | TB Green Premix Ex Taq II | Takara | 1200/200rxn |
| RNA | PrimeScript RT kit | Takara | 1500/100rxn |
| 蛋白 | RIPA裂解液 | 碧云天 | 150/100mL |
| 蛋白 | BCA蛋白定量试剂盒 | 碧云天 | 300/kit |
| 蛋白 | ECL发光液 | Millipore | 800/100mL |
| 抗体 | POMC (CST #7074) | CST | 2800/100μL |
| 抗体 | PC2 (Abcam ab28540) | Abcam | 4200/100μg |
| 细胞 | Lipofectamine RNAiMAX | Invitrogen | 2500/1.5mL |
| ELISA | α-MSH EIA Kit | Phoenix Pharm | 5800/96 wells |
| ELISA | β-Endorphin EIA Kit | Phoenix Pharm | 5800/96 wells |
| ChIP | SimpleChIP Kit | CST #9003 | 3800/24 rxn |

---

## 附录B: 数据分析软件

| 用途 | 软件 | 版本 |
|------|------|------|
| qPCR分析 | Bio-Rad CFX Manager | ≥3.6 |
| WB定量 | ImageJ / Fiji | ≥1.54 |
| IF定量 | ImageJ / Fiji + Coloc 2 | ≥1.54 |
| ELISA | GraphPad Prism | ≥9.0 |
| 统计 | R + ggplot2 | ≥4.2 |
| RNA-seq | STAR + DESeq2 | latest |
| 图版制作 | Adobe Illustrator / Inkscape | — |

---

*SOP编写日期: 2026-05-22*
*下一步: 项目批准后, 在动物伦理审批通过后按SOP-01开始实验*
