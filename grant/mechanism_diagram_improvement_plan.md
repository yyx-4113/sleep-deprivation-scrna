# 机制假说图精修方案

## 当前状态
- 已有 `grant/mechanism_diagram_nsfc.png`（代码生成，v3版本）
- 12×8英寸，300dpi，白底
- 展示：睡眠剥夺触发 → 细胞核(Pomc转录↓) → 细胞质(POMC→Pcsk2→肽↓) → 下游表型

## 建议改进项

### 1. 加GRNBoost2重要性数值（强化计算预测证据）
目前环路箭头只标了"[X] 环路断裂"，建议在箭头上标注：
- Pomc → Pcsk2: importance=0.766 (Top 1靶基因)
- Pcsk2 → Pomc: importance=0.510 (Top 1调控因子)

### 2. 加Pcsk2启动子TF结合预测结果（新增前期数据）
在细胞核区域，在Pcsk2基因位点旁加注预测的TF结合位点：
- "CREB1: 105 sites; Fos: 69 sites; Junb: 216 sites"
- 或用小标注框显示 "-291 CRE, -607 TRE" 关键位点

### 3. 明确展示"双打击"概念
当前是单线路径，建议加两条汇聚箭头：
- 打击1: Pomc mRNA↓ → POMC蛋白↓ (转录层面)
- 打击2: Pcsk2↓ → PC2酶↓ → POMC加工↓ (翻译后加工层面)
- 两条汇聚到"功能性肽崩溃"

### 4. 加虚拟KO预测靶基因
在Pcsk2下方加一行小字标注预测受影响的下游基因：
"预测靶基因: Scg2(-0.32), Ccnd2(-0.38), Dbp(↓)"

### 5. 六基因生物标志物
在下游结局区域加注：
"AUC=0.931 | Dbp+Nr3c1+Pomc+Rbm3+Rnaset2a+Tsc22d3"

### 6. 样式优化
- 中文字体渲染（SimHei/微软雅黑代替默认sans-serif）
- 反馈环路箭头颜色：正常=绿色虚线，断裂=红色粗线+大X
- 尽量用 BioRender 风格配色（比当前matplotlib默认色更专业）

## 建议操作
1. 短期（提交前）：代码快速加上#1、#2、#4的标注文字
2. 中期（如中标后发表）：请人用BioRender/Illustrator重绘为发表级
