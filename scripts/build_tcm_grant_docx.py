"""
Build the Fujian TCM Administration Grant Application as a DOCX file.
"""
from docx import Document
from docx.shared import Pt, Cm, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import os

OUT_DIR = 'D:/药膳'
doc = Document()

# Page setup
for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

style = doc.styles['Normal']
font = style.font
font.name = 'SimSun'
font.size = Pt(12)
style.paragraph_format.line_spacing = 1.5

def add_heading_styled(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = 'SimHei'
        if level == 0: run.font.size = Pt(16)
        elif level == 1: run.font.size = Pt(15)
        elif level == 2: run.font.size = Pt(14)
        elif level == 3: run.font.size = Pt(13)
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(text, bold=False, indent=False, alignment=None, size=12):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.5
    if indent:
        p.paragraph_format.first_line_indent = Cm(0.74)
    run = p.add_run(text)
    run.font.name = 'SimSun'
    run.font.size = Pt(size)
    run.bold = bold
    if alignment:
        p.alignment = alignment
    return p

def add_blank():
    doc.add_paragraph()

def add_rich_para(segments, indent=False):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.5
    if indent:
        p.paragraph_format.first_line_indent = Cm(0.74)
    for seg in segments:
        text = seg[0]
        bold = seg[1] if len(seg) > 1 else False
        italic = seg[2] if len(seg) > 2 else False
        fn = seg[3] if len(seg) > 3 else 'SimSun'
        run = p.add_run(text)
        run.font.name = fn
        run.font.size = Pt(12)
        run.bold = bold
        run.italic = italic
    return p

def add_table_with_data(headers, rows):
    ncols = len(headers)
    nrows = len(rows) + 1
    table = doc.add_table(rows=nrows, cols=ncols, style='Table Grid')
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, h in enumerate(headers):
        cell = table.rows[0].cells[j]
        cell.text = ''
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(h)
        run.font.name = 'SimHei'
        run.font.size = Pt(10)
        run.bold = True
        shading = cell._element.get_or_add_tcPr()
        shd = shading.makeelement(qn('w:shd'), {qn('w:fill'): 'D9E2F3', qn('w:val'): 'clear'})
        shading.append(shd)
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            cell = table.rows[i+1].cells[j]
            cell.text = ''
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER
            run = p.add_run(str(val))
            run.font.name = 'SimSun'
            run.font.size = Pt(9)
    add_blank()

def add_ref(text):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.3
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(10)

# ================================================================
# TITLE
# ================================================================
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
title.paragraph_format.space_after = Pt(4)
run = title.add_run('福建省中医药管理局 — 中医药科研课题申请书')
run.font.name = 'SimHei'
run.font.size = Pt(18)
run.bold = True
add_blank()

# ---- Cover Info ----
add_heading_styled('封面信息', 1)

cover_data = [
    ['项目名称', '基于AI虚拟筛选与网络扰动的闽产药食同源食材助眠活性肽发现及机制研究'],
    ['申报学科', '中药药理学 / 中医药信息学'],
    ['研究期限', '2027年1月 — 2028年12月（2年）'],
    ['申请经费', '3万元'],
    ['项目负责人', '杨永新'],
    ['依托单位', '福建省第二人民医院'],
    ['协作单位', '福建医科大学药学院（拟协作，提供HT22细胞实验平台支持）'],
]
add_table_with_data(['项目', '内容'], cover_data)
doc.add_page_break()

# ================================================================
# 一、摘要
# ================================================================
add_heading_styled('一、摘要（限300字）', 1)

add_para(
    '药食同源食材中的生物活性肽是安全温和的睡眠改善候选分子，但传统筛选耗时昂贵且缺乏分子'
    '到转录组的系统药效预测手段。本项目以三种闽产药食同源食材（建瓯核桃、武夷桑葚、闽北黑芝麻）'
    '为对象，建立"虚拟酶切→AI筛选→分子对接→虚拟细胞验证→体外初筛"的助眠肽高效发现管道：'
    '（1）UniProt获取三种食材蛋白组，虚拟胃肠酶切构建理论肽段库；'
    '（2）BiLSTM深度学习结合随机森林预测神经保护/助眠活性肽；'
    '（3）以GABA_AR、OX_2R、MC4R等5个睡眠相关受体为靶点批量分子对接筛选；'
    '（4）整合GSE137665单细胞数据构建下丘脑Pomc神经元虚拟细胞，网络传播算法预测候选肽段对'
    '睡眠剥夺转录组的逆转效应；'
    '（5）委托第三方合成TOP3候选肽段，在HT22细胞中进行CCK-8和ROS活性初步验证。'
    '本研究将为闽产药食同源资源的深度开发提供分子层面科学依据。',
    indent=True
)

add_blank()
add_rich_para([
    ('关键词：', True, False, 'SimHei'),
    ('药食同源、生物活性肽、虚拟筛选、深度学习、睡眠', False, False, 'SimSun'),
])
add_blank()

# ================================================================
# 二、立项依据
# ================================================================
add_heading_styled('二、立项依据', 1)

add_heading_styled('2.1 药食同源与睡眠健康：中医药现代化的突破口', 2)
add_para(
    '"药食同源"是中医药理论的核心思想之一，《黄帝内经》云"空腹食之为食物，患者食之为药物"。'
    '国家卫生健康委公布的药食同源目录已纳入110余种物质，为功能性食品和天然药物开发提供了'
    '丰富资源。在睡眠健康领域，酸枣仁、百合、莲子、龙眼肉等药食同源食材的安神助眠功效已被'
    '《中国药典》收录，酸枣仁百合睡眠肽复合物产品已上市[1]。然而，这些食材中的特定活性肽段'
    '是什么？它们通过哪些分子靶点起作用？能否通过计算方法高通量预测？这些问题尚未被系统回答。',
    indent=True
)

add_heading_styled('2.2 食源性生物活性肽的计算发现：从传统到AI', 2)
add_para(
    '传统的生物活性肽（BAPs）发现遵循"蛋白酶解→分离纯化→活性检测→鉴定"的路径，单个活性肽'
    '的发现通常需要数月时间和数万元成本。近年来，虚拟筛选+分子对接的计算策略已显著加速这一过程'
    '[2]。Mo等（2025）通过虚拟酶切结合肽组学，从山羊乳酪蛋白中鉴定出118个肽段，并利用分子'
    '对接筛选到助眠肽FAWPQY（GABA_AR结合能-10.2 kcal/mol），体内验证延长小鼠睡眠时间至37'
    '分钟[3]。Qian等（2024）利用随机森林模型建立了酪蛋白助眠肽的结构-活性关系，发现了YQKFPQY'
    '等3个新肽段[4]。2024年底，Chang等在Trends in Food Science & Technology（IF 15.3）系统'
    '总结了AI在食源性肽段筛选中的应用框架，指出深度学习（CNN/LSTM/Transformer）在抗炎、抗氧化、'
    '降血压肽筛选中显著优于传统机器学习，且多尺度化学空间特征的通用模型是未来方向[5]。',
    indent=True
)

add_heading_styled('2.3 福建省药食同源资源的开发潜力', 2)
add_para(
    '福建省拥有丰富的药食同源资源，其中建瓯核桃（Juglans regia L.）、武夷山桑葚（Morus alba L.）'
    '和闽北黑芝麻（Sesamum indicum L.）种植面积广、产量大。2024-2025年，核桃肽的神经保护作用'
    '已被多篇高质量论文验证：TWLPLPR肽可穿透血脑屏障保护突触功能[6]、SSA-LSTM深度学习模型'
    '以98%准确率识别核桃神经保护肽[7]、WSPSGR肽通过MLCK通路改善PC12细胞氧化应激（ABTS清除率'
    '94.83%）[8]。黑芝麻ACE抑制肽ITAPHW和LLLPYY已被鉴定[9]。桑葚肽的研究相对薄弱，仅有一篇'
    '2024年Food Chemistry Advances报道了白桑蛋白水解物的ACE抑制和抗氧化活性[10]，但任何药食'
    '同源食材的助眠肽AI筛选，目前均未见报道。',
    indent=True
)

add_heading_styled('2.4 前期研究基础：已具备单细胞虚拟细胞模型', 2)
add_para(
    '申请人前期已对睡眠剥夺小鼠的单细胞转录组数据（GSE137665[11]，24,778个细胞，涵盖下丘脑/脑干/皮层）'
    '完成了系统分析，成果发表于bioRxiv预印本[12]。关键前期成果包括：（1）发现Pomc为下丘脑睡眠剥夺的'
    '最显著差异基因（log2FC=-4.08，padj=1.3×10⁻¹⁵⁷）；（2）构建了Pomc为中心的基因调控网络'
    '（35节点）；（3）建立了网络传播虚拟敲除算法，预测Pomc缺失的下游效应；（4）用梯度提升回归'
    '（GRNBoost2[16]）验证了Pomc-Pcsk2正反馈环路。这些前期工作为本项目构建"虚拟细胞"模型提供了'
    '直接的计算基础设施。',
    indent=True
)
add_para(
    '此外，申请人近期完成了三种闽产药食同源食材的AI虚拟筛选预实验[17]，验证了"虚拟酶切→AI筛选'
    '→分子对接"管道的技术可行性。具体成果如下：从UniProt获取核桃（500条）、黑芝麻（800条）、'
    '桑葚（776条）三种食材蛋白序列，经虚拟胃肠酶切（pepsin/trypsin/chymotrypsin/papain）构建了'
    '共451,785条理论肽段库（核桃135,577 + 黑芝麻224,149 + 桑葚92,059），利用随机森林模型'
    '（AUC=0.874, F1=0.618）预测神经保护活性。三种食材中均鉴定到富含芳香族氨基酸（Y/W/F/P）'
    '的高评分候选肽段（如PFYY、FYPY、YPW等），为后续7靶点分子对接精筛（420次）和HT22细胞'
    '验证提供了明确的候选池和完整的可复用代码框架。',
    indent=True
)

add_heading_styled('2.5 科学假说', 2)
add_para(
    '闽产药食同源食材（核桃、桑葚、黑芝麻）的蛋白经胃肠酶切后释放的生物活性肽，可通过多靶点'
    '（GABA_AR、OX_2R、MC4R、5-HT1A、Keap1/Nrf2）协同调控下丘脑Pomc神经元的基因网络，'
    '部分逆转睡眠剥夺造成的转录组紊乱。',
    indent=True, bold=True
)

add_heading_styled('参考文献', 2)
refs = [
    '[1] 酸枣仁百合睡眠肽复合物改善失眠模型大鼠的睡眠. 现代食品科技, 2021.',
    '[2] Singh BP, et al. In silico and molecular docking approaches in food-derived bioactive peptide discovery. Food Res Int, 2025.',
    '[3] Mo L, et al. Identification of sleep-enhancing peptides from goat milk casein. Food Res Int, 2025, 221:117333.',
    '[4] Qian J, et al. Exploring structural features of sleep-enhancing peptides from casein hydrolysates. Food Chem, 2024, 461:140838.',
    '[5] Chang J, et al. AI in food bioactive peptides screening. Trends Food Sci Technol, 2024, 156:104845.',
    '[6] Min W, et al. Walnut-derived peptides cross the blood-brain barrier and ameliorate Abeta-induced hypersynchronous neural network activity. Food Res Int, 2024, 197(2):115302.',
    '[7] Lin L, et al. SSA-LSTM-VD: An in silico scheme for optimizing the enzymatic acquisition of natural biologically active peptides based on machine learning and virtual digestion. Anal Chim Acta, 2024, 1298:342419.',
    '[8] Li J, et al. Structure-activity relationship and antioxidant properties of walnut peptides with high intestinal absorption capacity. Food Funct, 2025, 16:9063-9079.',
    '[9] Du T, et al. ACE inhibitory peptides from enzymatic hydrolysate of fermented black sesame seed. Food Chem, 2024, 437:137921.',
    '[10] Kirac FT, et al. Isolation and characterization of a new potential source of bioactive peptides: White mulberry (Morus alba) fruits and its leaves. Food Chem Adv, 2024, 4:100597.',
    '[11] Jha PK, et al. Single-cell transcriptomics and cell-specific proteomics reveals molecular signatures of sleep. Commun Biol, 2022, 5:846.',
    '[12] Yang YX. Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub. bioRxiv, 2026 (manuscript in preparation).',
    '[13] Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. J Comput Chem, 2010, 31:455-461.',
    '[14] Abraham MJ, et al. GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX, 2015, 1-2:19-25.',
    '[15] Morimoto BH, Koshland DE Jr. Excitatory amino acid uptake and N-methyl-D-aspartate-mediated secretion in a neural cell line. Proc Natl Acad Sci USA, 1990, 87:3518-3521.',
    '[16] Aibar S, et al. SCENIC: single-cell regulatory network inference and clustering. Nat Methods, 2017, 14:1083-1086.',
    '[17] Yang YX. AI-driven virtual screening of sleep-enhancing bioactive peptides from three Fujian medicinal-food homologous ingredients (Juglans regia, Sesamum indicum, Morus alba). bioRxiv, 2026 (manuscript in preparation).',
]
for r in refs:
    add_ref(r)

doc.add_page_break()

# ================================================================
# 三、研究目标
# ================================================================
add_heading_styled('三、研究目标', 1)
add_para('总体目标：建立一套"AI驱动的药食同源助眠活性肽高效发现管道"，并以三种闽产食材为实例完成筛选和初步验证。', indent=True, bold=True)
add_blank()
add_para('具体目标：', bold=True)
add_para('1. 构建三种食材的虚拟酶切肽段库（预计5,000+条），利用深度学习模型预测助眠/神经保护活性', indent=True)
add_para('2. 完成候选肽段与5个睡眠相关受体的分子对接筛选', indent=True)
add_para('3. 利用虚拟细胞模型预测候选肽段对睡眠剥夺转录组的逆转效应', indent=True)
add_para('4. 委托第三方合成TOP 3肽段，在HT22细胞中完成初步活性验证', indent=True)
add_blank()

# ================================================================
# 四、研究内容
# ================================================================
add_heading_styled('四、研究内容', 1)

# Content 1
add_heading_styled('内容一：蛋白组获取与虚拟酶切（第1-3月）', 2)
add_rich_para([('(1) UniProt蛋白序列获取', True)])
add_para('检索关键词：Juglans regia、Morus alba、Sesamum indicum', indent=True)
add_para('获取三种食材的全部已注释蛋白序列，筛选条件：蛋白长度 > 50 AA', indent=True)
add_rich_para([('(2) 虚拟胃肠酶切', True)])
add_para('模拟消化酶：胃蛋白酶（pepsin）、胰蛋白酶（trypsin）、糜蛋白酶（chymotrypsin）、木瓜蛋白酶（papain）', indent=True)
add_para('酶切参数：允许最多2个漏切位点；肽段长度2-20 AA；去冗余', indent=True)
add_para('预计产出：每种食材1,000-2,000条理论肽段，合计5,000+条', indent=True)

# Content 2
add_heading_styled('内容二：AI驱动的生物活性预测（第3-6月）', 2)
add_rich_para([('(1) 特征工程', True)])
add_para('氨基酸组成（AAC）+ 伪氨基酸组成（PseAAC）+ 二肽组成 → ~440维特征向量', indent=True)
add_rich_para([('(2) 传统机器学习模型', True)])
add_para('随机森林（RF）：复现Qian 2024的酪蛋白肽段QSAR框架', indent=True)
add_rich_para([('(3) 深度学习模型', True)])
add_para('BiLSTM + Attention机制；参考SSA-LSTM-VD架构（98%准确率）', indent=True)
add_rich_para([('(4) 生物活性多标签预测', True)])
add_para('预测活性：神经保护、抗氧化、ACE抑制、抗炎 → 每种食材TOP 20（合计60条候选）', indent=True)

# Content 3
add_heading_styled('内容三：多靶点分子对接筛选（第6-9月）', 2)
add_rich_para([('(1) 受体结构准备', True)])
target_rows = [
    ['GABA_AR', '6D6T', '抑制性神经递质受体'],
    ['OX_2R', '7TAC', '食欲素受体2'],
    ['MC4R', '6W25', '黑皮质素4受体'],
    ['5-HT1A', '7E2Z', '血清素受体'],
    ['Keap1', '2FLU', 'Nrf2抑制蛋白'],
]
add_table_with_data(['受体', 'PDB ID', '功能'], target_rows)
add_rich_para([('(2) 对接流程', True)])
add_para('软件：AutoDock Vina[13]；结合能 < -8 kcal/mol为强结合候选', indent=True)
add_rich_para([('(3) 分子动力学模拟验证（TOP 5肽段）', True)])
add_para('GROMACS[14] 50ns模拟；分析RMSD/RMSF/MM-PBSA结合自由能', indent=True)

# Content 4
add_heading_styled('内容四：虚拟细胞模型与网络扰动预测（第9-12月）', 2)
add_rich_para([('(1) Pomc神经元虚拟细胞构建', True)])
add_para('基于GSE137665下丘脑scRNA-seq数据；提取Pomc⁺神经元→三种状态GRN', indent=True)
add_rich_para([('(2) 虚拟肽段给药模拟', True)])
add_para('肽段-受体结合强度→TF活性变化→GRN边权重传播→全基因表达变化向量', indent=True)
add_rich_para([('(3) 逆转分数计算', True)])
add_para('计算虚拟肽段转录组 vs 睡眠剥夺转录组的余弦距离；逆转分数 > 0.5为有效', indent=True)

# Content 5
add_heading_styled('内容五：体外实验初步验证（第12-18月）', 2)
add_rich_para([('(1) 委托第三方合成TOP 3肽段', True)])
add_para('综合AI评分+分子对接+逆转分数排名前三的肽段；纯度>95%（HPLC）；每条10mg', indent=True)
add_rich_para([('(2) HT22细胞活力检测', True)])
add_para('HT22小鼠海马神经元[15]；H₂O₂诱导氧化应激（200μM, 24h）', indent=True)
add_para('肽段梯度浓度（10/50/100 μM）；CCK-8细胞活力+ROS荧光探针', indent=True)

doc.add_page_break()

# ================================================================
# 机制示意图
# ================================================================
add_heading_styled('研究假说与机制示意图', 1)
add_para('下图概括了本课题的整体研究设计：以三种闽产药食同源食材（建瓯核桃、武夷桑葚、闽北黑芝麻）为起点，通过"虚拟酶切->AI筛选->分子对接->虚拟细胞验证->体外实验"五步管道，系统发现并验证具有多靶点协同作用的助眠活性肽。', indent=True)
add_blank()

# Insert mechanism diagram image
mech_img_path = 'D:/药膳/mechanism_diagram_tcm.png'
if os.path.exists(mech_img_path):
    p_img = doc.add_paragraph()
    p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run_img = p_img.add_run()
    run_img.add_picture(mech_img_path, width=Cm(14))
    add_blank()
    add_para('图1. AI驱动的闽产药食同源食材助眠活性肽高效发现管道示意图', indent=True,
             alignment=WD_ALIGN_PARAGRAPH.CENTER, size=9)
else:
    add_para(f'[机制图: {mech_img_path} — 请确保已运行 generate_mechanism_diagrams.py 生成]', indent=True)

doc.add_page_break()

# ================================================================
# 五、技术路线
# ================================================================
add_heading_styled('五、技术路线', 1)

# ---- Embed roadmap image (template style) ----
add_blank()
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
run = p.add_run()
run.add_picture('D:/药膳/roadmap_tcm.png', width=Inches(5.5))
add_blank()

doc.add_page_break()

# ================================================================
# 六、关键科学问题
# ================================================================
add_heading_styled('六、关键科学问题', 1)
add_para('1. 深度学习能否从虚拟酶切肽段库中准确识别具有助眠/神经保护活性的肽段？即AI模型对"未见食材"的泛化能力。', indent=True)
add_para('2. 药膳肽段是否能通过多靶点协同（GABA_AR镇静 + MC4R促Pomc恢复 + Keap1抗氧化）实现睡眠剥夺转录组的"多通路逆转"？即中药"多成分多靶点"理论在肽段层面的计算验证。', indent=True)
add_para('3. 虚拟细胞模型预测的转录组效应能否被体外实验（至少部分）验证？即计算预测与实验的一致性评估。', indent=True)
add_blank()

# ================================================================
# 七、可行性分析
# ================================================================
add_heading_styled('七、可行性分析', 1)

add_heading_styled('7.1 前期基础扎实', 2)
add_para('已完成GSE137665全部分析，拥有成熟的GRN构建和虚拟敲除算法；成果已发表于bioRxiv预印本；团队具备完整的Python/R生信分析能力', indent=True)

add_heading_styled('7.2 数据和方法全部公开免费', 2)
add_para('UniProt、GEO、PDB均为免费公开数据库；AutoDock Vina、GROMACS、Scanpy、PyTorch均为免费开源软件；Google Colab提供免费GPU资源', indent=True)

add_heading_styled('7.3 技术路线成熟', 2)
add_para('虚拟酶切方法是肽段研究的成熟范式；AI筛选方法有SSA-LSTM-VD（2024）可直接参考；分子对接是食源性BAPs发现的标配工具；HT22细胞模型操作简单，结果重复性好', indent=True)

add_heading_styled('7.4 风险与应对', 2)
risk_rows = [
    ['UniProt中三种食材蛋白序列较少', '使用同属近缘物种序列补充；考虑转录组组装辅助'],
    ['AI模型训练数据不足', '数据增强（保守氨基酸替换）；迁移学习（预训练蛋白语言模型）'],
    ['分子对接对线性肽准确度有限', '用MD模拟补充验证；与已知活性肽FAWPQY横向对比'],
    ['委托合成肽段溶解性差', '合成前做溶解度预测；必要时做N端/C端修饰'],
    ['HT22细胞中肽段活性不明显', '增加浓度范围；考虑不同应激模型（谷氨酸盐/缺氧）'],
]
add_table_with_data(['潜在风险', '应对措施'], risk_rows)

add_blank()

# ================================================================
# 八、创新点
# ================================================================
add_heading_styled('八、创新点', 1)
add_para('1. 首次将AI深度学习（BiLSTM+Attention）应用于药食同源食材的助眠活性肽筛选：现有AI筛选集中在抗炎/抗氧化/降压，睡眠方向完全空白', indent=True)
add_para('2. 首次构建"虚拟酶切→AI筛选→多靶点对接→虚拟细胞验证→体外初筛"的完整计算管道：将食品科学+计算系统生物学+睡眠医学交叉融合', indent=True)
add_para('3. 首次提出"逆转分数(Rescue Score)"概念：定量评估候选肽段对疾病转录组的网络层面逆转效应', indent=True)
add_para('4. 聚焦福建地域特色资源：建瓯核桃+武夷桑葚+闽北黑芝麻，为"闽产药膳"提供分子层面的现代科学依据', indent=True)

doc.add_page_break()

# ================================================================
# 九、预期成果
# ================================================================
add_heading_styled('九、预期成果', 1)
outcome_rows = [
    ['SCI论文', '食品科学/中医药领域期刊（J Funct Foods / J Ethnopharmacol / Food Funct等, IF 3-6）', '1-2篇'],
    ['候选肽段', '经计算+体外验证的助眠活性肽序列', '3条'],
    ['计算管道', '可复用的"药膳肽段AI筛选"全套代码', '1套'],
    ['数据库', '三种闽产食材的虚拟酶切肽段库', '5000+条'],
]
add_table_with_data(['成果类型', '具体内容', '数量/指标'], outcome_rows)
add_blank()

# ================================================================
# 十、研究基础
# ================================================================
add_heading_styled('十、研究基础', 1)
add_heading_styled('10.1 已发表/预发表的相关成果', 2)
add_para('1. Yang YX. Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub. bioRxiv, 2026 (manuscript in preparation).', indent=True)
add_para('2. Yang YX. AI-driven virtual screening of sleep-enhancing bioactive peptides from three Fujian medicinal-food homologous ingredients (Juglans regia, Sesamum indicum, Morus alba). bioRxiv, 2026 (manuscript in preparation).', indent=True)
add_para('3. （申请人其他已发表论文）', indent=True)
add_heading_styled('10.2 已有的计算分析条件', 2)
add_para('已建立的单细胞RNA-seq全流程分析管道；已建立的GRN推断和网络传播虚拟敲除算法；已建立的机器学习生物标志物筛选流程；已建立的药食同源肽段AI虚拟筛选管道（UniProt→虚拟酶切→RF预测→分子对接）', indent=True)
add_heading_styled('10.3 协作条件', 2)
add_para('本课题纯计算部分在现有设备即可完成；肽段合成委托第三方公司；细胞实验可与省内高校实验室协作', indent=True)
add_blank()

# ================================================================
# 十一、经费预算
# ================================================================
add_heading_styled('十一、经费预算', 1)
budget_rows = [
    ['一、设备费', '0', '使用现有计算机'],
    ['二、材料费', '1.8', ''],
    ['  多肽合成（3条）', '1.2', '3条×4000元/条（>95% HPLC，含溶解度测试和质谱鉴定）'],
    ['  细胞实验试剂', '0.6', 'HT22细胞系、CCK-8试剂盒、ROS探针、培养基等'],
    ['三、测试化验加工费', '0', '—'],
    ['四、差旅/会议', '0.3', '国内学术会议1次'],
    ['五、论文发表', '0.5', '期刊版面费（优先非OA或利用单位发表基金）'],
    ['六、劳务费', '0.4', '科研辅助'],
    ['合计', '3.0', ''],
]
add_table_with_data(['科目', '金额（万元）', '计算依据'], budget_rows)
add_blank()

# ================================================================
# 十二、年度计划
# ================================================================
add_heading_styled('十二、年度计划', 1)
plan_rows = [
    ['第一年\n(2027)', '1-3月', 'UniProt蛋白组获取 + 虚拟酶切', '三种食材肽段库'],
    ['', '3-6月', 'BiLSTM+RF模型构建与训练', 'AI预测评分表'],
    ['', '6-9月', '5靶点分子对接（60肽×5靶点）', '对接结合能矩阵'],
    ['', '9-12月', '虚拟细胞模型+网络扰动+逆转分数', '虚拟筛选完整结果'],
    ['第二年\n(2028)', '1-6月', 'TOP 3肽段合成 + HT22 CCK-8/ROS', '体外验证数据'],
    ['', '7-12月', '数据分析整合+论文撰写+投稿+结题', 'SCI论文1篇'],
]
add_table_with_data(['年度', '时间', '研究内容', '考核指标'], plan_rows)
add_blank()

# ================================================================
# 十三、附件清单
# ================================================================
add_heading_styled('十三、附件清单', 1)
add_para('1. 申请人简历（含近5年发表论文列表）', indent=True)
add_para('2. 前期研究代表性成果（bioRxiv预印本全文）', indent=True)
add_para('3. 依托单位实验平台条件证明', indent=True)
add_para('4. 多肽合成第三方公司报价单（可选）', indent=True)
add_para('5. 协作单位意向书（如有细胞实验协作）', indent=True)

# ================================================================
# 十四、数据与代码可获取性声明
# ================================================================
add_heading_styled('十四、数据与代码可获取性声明', 1)
add_para('本项目所有计算分析和软件工具均基于开源平台，确保研究可复现性：', indent=True)
add_para('1. 三种食材虚拟酶切肽段库及AI预测评分将在GitHub公开', indent=True)
add_para('2. 虚拟细胞模型代码及网络传播算法将以Python包形式发布', indent=True)
add_para('3. 后续产生的RNA-seq数据（如有）将上传至GEO公共数据库', indent=True)
add_para('4. 项目预注册于OSF（Open Science Framework）', indent=True)

add_blank()
p = doc.add_paragraph()
p.paragraph_format.line_spacing = 1.3
run = p.add_run('* 本申请书基于2项bioRxiv预印本前期发现（睡眠剥夺单细胞分析 & 核桃肽AI虚拟筛选）及Chang 2024 (Trends Food Sci Technol) AI筛选框架设计。')
run.font.name = 'SimSun'
run.font.size = Pt(10)
run.italic = True

# ---- Save ----
os.makedirs(OUT_DIR, exist_ok=True)
output_path = f'{OUT_DIR}/福建省中医药管理局_课题申请书.docx'
doc.save(output_path)
print(f'Grant application DOCX saved to: {output_path}')
print(f'Grant application DOCX saved to: {output_path}')
