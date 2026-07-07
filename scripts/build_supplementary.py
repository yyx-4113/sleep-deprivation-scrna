"""
Build supplementary materials DOCX with all supplementary figures and tables.
"""
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
import os

FIG_DIR = 'results/figures'
TAB_DIR = 'results/tables'
SUPP_DIR = 'bioRxiv_submission'

doc = Document()

for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(11)
style.paragraph_format.line_spacing = 1.5

def add_heading_styled(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(text, bold=False, italic=False, alignment=None, size=11):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 1.5
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic
    if alignment:
        p.alignment = alignment
    return p

def add_figure(img_path, width_cm=14, caption=''):
    if not os.path.exists(img_path):
        add_para(f'[Figure missing: {os.path.basename(img_path)}]', italic=True)
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.space_after = Pt(4)
    run = p.add_run()
    run.add_picture(img_path, width=Cm(width_cm))
    if caption:
        cp = doc.add_paragraph()
        cp.paragraph_format.line_spacing = 1.3
        cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run_c = cp.add_run(caption)
        run_c.font.name = 'Times New Roman'
        run_c.font.size = Pt(9)
        run_c.bold = True

def add_table_from_csv(csv_path, caption='', max_rows=None):
    import pandas as pd
    if not os.path.exists(csv_path):
        add_para(f'[Table missing: {os.path.basename(csv_path)}]', italic=True)
        return
    df = pd.read_csv(csv_path)
    if max_rows and len(df) > max_rows:
        df = df.head(max_rows)
    if caption:
        cp = doc.add_paragraph()
        cp.paragraph_format.space_before = Pt(8)
        cp.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = cp.add_run(caption)
        run.font.name = 'Times New Roman'
        run.font.size = Pt(9)
        run.bold = True
    rows, cols = len(df) + 1, len(df.columns)
    table = doc.add_table(rows=rows, cols=cols, style='Table Grid')
    from docx.enum.table import WD_TABLE_ALIGNMENT
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for j, col_name in enumerate(df.columns):
        cell = table.rows[0].cells[j]
        cell.text = ''
        p = cell.paragraphs[0]
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = p.add_run(str(col_name))
        run.font.name = 'Times New Roman'
        run.font.size = Pt(7)
        run.bold = True
        shading = cell._element.get_or_add_tcPr()
        shading_elm = shading.makeelement(qn('w:shd'), {qn('w:fill'): 'D9E2F3', qn('w:val'): 'clear'})
        shading.append(shading_elm)
    for i, (_, row) in enumerate(df.iterrows()):
        for j, val in enumerate(row):
            cell = table.rows[i+1].cells[j]
            cell.text = ''
            p = cell.paragraphs[0]
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run = p.add_run(f'{val}' if not isinstance(val, float) else f'{val:.4g}')
            run.font.name = 'Times New Roman'
            run.font.size = Pt(7)
    doc.add_paragraph()

# ================================================================
add_heading_styled('Supplementary Information', 1)
add_para(
    'This document contains supplementary figures and tables accompanying the manuscript:\n'
    '"Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central '
    'Regulatory Hub and Predicts Downstream Transcriptional Consequences of Its Loss"'
)
add_para('')

# -- Supplementary Figures --
add_heading_styled('Supplementary Figures', 1)

add_figure(f'{FIG_DIR}/CellType_Composition.png', 14,
           'Figure S1. Cell type composition across brain regions and conditions.')
add_figure(f'{FIG_DIR}/Marker_Dotplot.png', 14,
           'Figure S2. Marker gene dotplot for cell type annotation.')
add_figure(f'{FIG_DIR}/Hypothalamus_UMAP.png', 12,
           'Figure S3. Hypothalamus UMAP with cluster annotations.')
add_figure(f'{FIG_DIR}/Top_DEGs_Heatmap.png', 14,
           'Figure S4. Heatmap of top DEGs across all brain regions and conditions.')
add_figure(f'{FIG_DIR}/Response_by_Region.png', 12,
           'Figure S5. Differential expression response comparison by brain region.')
add_figure(f'{FIG_DIR}/TF_Activity_Network.png', 14,
           'Figure S6. TF activity regulatory network around Pomc.')
add_figure(f'{FIG_DIR}/Validation_GSE211088_Barplot.png', 12,
           'Figure S7. Key gene expression in GSE211088 validation dataset.')
add_figure(f'{FIG_DIR}/Validation_GSE237419_Barplot.png', 12,
           'Figure S8. Key gene expression in GSE237419 validation dataset.')
add_figure(f'{FIG_DIR}/Validation_ourStudy_Barplot.png', 12,
           'Figure S9. Key gene expression in our study (GSE137665 hypothalamus).')

doc.add_page_break()

# -- Supplementary Tables --
add_heading_styled('Supplementary Tables', 1)
add_para('Tables S1-S13 are also provided as machine-readable CSV files in the supplementary_tables/ folder.')

add_table_from_csv(f'{TAB_DIR}/Top10_DEGs_Hypothalamus.csv',
                   'Table S1. Top 10 differentially expressed genes: Hypothalamus (A2 vs A1).')
add_table_from_csv(f'{TAB_DIR}/Top10_DEGs_Brainstem.csv',
                   'Table S2. Top 10 differentially expressed genes: Brainstem (A2 vs A1).')
add_table_from_csv(f'{TAB_DIR}/Top10_DEGs_Cortex.csv',
                   'Table S3. Top 10 differentially expressed genes: Cortex (A2 vs A1).')
add_table_from_csv(f'{TAB_DIR}/Pomc_coexpression.csv',
                   'Table S4. Pomc co-expression results (top 20 genes by Spearman |r|).', max_rows=20)
add_table_from_csv(f'{TAB_DIR}/Pomc_network_TFs.csv',
                   'Table S5. Transcription factors in Pomc co-expression network (|r| > 0.3).')
add_table_from_csv(f'{TAB_DIR}/Pomc_knockout_prediction.csv',
                   'Table S6. In silico Pomc knockout prediction results.', max_rows=20)
add_table_from_csv(f'{TAB_DIR}/ML_consensus_biomarkers.csv',
                   'Table S7. Machine learning consensus biomarkers.')
add_table_from_csv(f'{TAB_DIR}/External_Validation_Results.csv',
                   'Table S8. External validation: cross-dataset log2FC comparison.', max_rows=20)
add_table_from_csv(f'{TAB_DIR}/Differential_TF_Activity.csv',
                   'Table S9. Differential transcription factor activity (all TFs).', max_rows=15)
add_table_from_csv(f'{TAB_DIR}/TF_regulons.csv',
                   'Table S10. TF regulon definitions (|r| > 0.25, p < 0.01).', max_rows=15)
add_table_from_csv(f'{TAB_DIR}/GRNBoost2_links.csv',
                   'Table S11. Gradient boosting regression (GRNBoost2-equivalent) regulatory links.', max_rows=15)
add_table_from_csv(f'{TAB_DIR}/Pomc_KO_GRNBoost2.csv',
                   'Table S12. GBR-informed in silico Pomc knockout prediction.', max_rows=15)

# Save
output_path = f'{SUPP_DIR}/supplementary_materials.docx'
doc.save(output_path)
print(f'Supplementary materials saved to: {output_path}')
