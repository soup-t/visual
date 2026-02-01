import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import json
import random
from typing import Optional
import requests

# 设置页面配置
st.set_page_config(page_title="Rice GWAS Interactive Visualization", layout="wide")

# --- AI分析函数 ---
def generate_ai_analysis(filtered_df, total_genes, favorable, unfavorable, neutral_unknown, 
                        selected_category, segments_df, lang_code):
    """调用阿里云百炼API生成分析报告"""
    API_KEY = "sk-1308247f1cdc494daf93c4b999f36cf1"
    API_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1/chat/completions"
    
    try:
        # 准备数据摘要
        trait_counts = filtered_df['Trait_Category'].value_counts().head(10).to_dict()
        
        # 获取各类基因详情
        favorable_by_trait = {}
        unfavorable_by_trait = {}
        
        for trait in filtered_df['Trait_Category'].unique()[:8]:
            trait_df = filtered_df[filtered_df['Trait_Category'] == trait]
            fav_genes = trait_df[trait_df['Evaluation'] == 'Favorable']['GeneName'].tolist()[:5]
            unfav_genes = trait_df[trait_df['Evaluation'] == 'Unfavorable']['GeneName'].tolist()[:5]
            
            if fav_genes:
                favorable_by_trait[trait] = fav_genes
            if unfav_genes:
                unfavorable_by_trait[trait] = unfav_genes
        
        # 祖源信息
        ancestry_summary = ""
        if segments_df is not None and not segments_df.empty:
            for hap in ['hap1', 'hap2']:
                hap_data = segments_df[segments_df['Haplotype'] == hap]
                if not hap_data.empty:
                    group_counts = hap_data['Group'].value_counts()
                    ancestry_summary += f"\n{hap}: "
                    for group_code, count in group_counts.head(3).items():
                        if pd.notna(group_code):
                            group_name = GROUP_NAME_MAP.get(int(group_code), str(group_code))
                            ancestry_summary += f"{group_name}, "
        
        # 构建提示词
        if lang_code == 'zh':
            prompt = f"""你是一位资深的水稻育种专家。请基于以下基因数据提供简洁专业的育种分析报告。

**样本基因组概况：**
- 总基因数：{total_genes}
- 有利基因：{favorable}个 ({favorable/total_genes*100:.1f}%)
- 不利基因：{unfavorable}个 ({unfavorable/total_genes*100:.1f}%)
- 中性/未知：{neutral_unknown}个

**性状分类分布：**
{json.dumps(trait_counts, ensure_ascii=False, indent=2)}

**有利基因分布（按性状）：**
{json.dumps(favorable_by_trait, ensure_ascii=False, indent=2)}

**不利基因分布（按性状）：**
{json.dumps(unfavorable_by_trait, ensure_ascii=False, indent=2)}

**祖源组成：**{ancestry_summary}

请严格按照以下格式输出分析报告：

## ✅ 核心优势 (3项)

### 1. [第一个优势标题]
• [具体基因名称]：[功能描述]
• [具体基因名称]：[功能描述]

### 2. [第二个优势标题]
• [具体基因名称]：[功能描述]

### 3. [第三个优势标题]
• [具体基因名称]：[功能描述]

## ⚠️ 待改良缺陷 (4项)

### 1. [第一个缺陷标题] ！
• [具体问题]：[具体基因和影响]
• [改良建议]

### 2. [第二个缺陷标题] ！！
• [具体问题]
• [风险说明]

### 3. [第三个缺陷标题] ！
• [具体问题]
• [改良建议]

### 4. [第四个缺陷标题]
• [具体问题]
• [建议]

要求：
1. 必须基于提供的实际基因数据，列出具体基因名称
2. 每个优势/缺陷要具体到基因功能和育种意义
3. 使用专业术语但保持简洁
4. 重要程度用！数量表示
5. 格式严格按照示例，使用Markdown语法
"""
        else:
            prompt = f"""You are a senior rice breeding expert. Please provide a concise professional breeding analysis based on the following gene data.

**Genome Overview:**
- Total genes: {total_genes}
- Favorable genes: {favorable} ({favorable/total_genes*100:.1f}%)
- Unfavorable genes: {unfavorable} ({unfavorable/total_genes*100:.1f}%)
- Neutral/Unknown: {neutral_unknown}

**Trait Distribution:**
{json.dumps(trait_counts, indent=2)}

**Favorable Genes by Trait:**
{json.dumps(favorable_by_trait, indent=2)}

**Unfavorable Genes by Trait:**
{json.dumps(unfavorable_by_trait, indent=2)}

**Ancestry:**{ancestry_summary}

Please provide analysis in this format:

## ✅ Core Advantages (3 items)

### 1. [First advantage title]
• [Specific genes]: [Function description]

### 2. [Second advantage]
• [Specific genes]: [Function]

### 3. [Third advantage]
• [Specific genes]: [Function]

## ⚠️ Areas for Improvement (4 items)

### 1. [First deficiency] !
• [Specific issue]: [Genes and impact]
• [Improvement suggestion]

### 2. [Second deficiency] !!
• [Issue]
• [Risk note]

### 3. [Third deficiency] !
• [Issue]
• [Suggestion]

### 4. [Fourth deficiency]
• [Issue]
• [Recommendation]

Requirements:
1. Must be based on actual gene data provided
2. List specific gene names and breeding significance
3. Use professional but concise language
4. Use ! to indicate importance level
5. Follow format strictly using Markdown
"""
        
        # 调用API
        headers = {
            "Authorization": f"Bearer {API_KEY}",
            "Content-Type": "application/json"
        }
        
        data = {
            "model": "qwen-plus",
            "messages": [
                {"role": "system", "content": "你是一位专业的水稻育种专家，擅长基因组数据分析和育种决策。" if lang_code == 'zh' else "You are a professional rice breeding expert specializing in genomic data analysis."},
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.7,
            "max_tokens": 2500
        }
        
        response = requests.post(API_URL, headers=headers, json=data, timeout=60)
        
        if response.status_code == 200:
            result = response.json()
            return result['choices'][0]['message']['content']
        else:
            return f"❌ API调用失败: {response.status_code}\n{response.text}"
            
    except Exception as e:
        return f"❌ 分析出错: {str(e)}"

# --- 加载翻译文件 ---
@st.cache_data
def load_translations():
    json_path = os.path.join(os.path.dirname(__file__), 'translations.json')
    if os.path.exists(json_path):
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}

translations = load_translations()

# --- 加载基因数据库 ---
@st.cache_data
def load_gene_db():
    json_path = os.path.join(os.path.dirname(__file__), 'gene_info.json')
    if os.path.exists(json_path):
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {}

gene_db = load_gene_db()

# --- 加载 QTN 效应值 ---
@st.cache_data
def load_qtn_effect():
    """加载 QTN_effect 文件，返回 {Chr_Pos: effect} 的字典"""
    file_path = os.path.join(os.path.dirname(__file__), 'QTN_effect')
    qtn_dict = {}
    if os.path.exists(file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines[1:]:  # 跳过表头
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        chr_val = parts[0].strip()
                        pos_val = parts[1].strip()
                        effect_val = parts[2].strip().lower()
                        # 构造 key: "Chr6_2940130"
                        key = f"{chr_val}_{pos_val}"
                        qtn_dict[key] = effect_val
        except Exception as e:
            st.error(f"读取 QTN_effect 文件时出错: {e}")
    return qtn_dict

qtn_effects = load_qtn_effect()

# --- 加载 QTN 基因位置数据 ---
@st.cache_data
def load_qtn_genes():
    """加载 QTN 文件，返回 {Chr: [(GeneName, Position_Mb)]} 的字典"""
    file_path = os.path.join(os.path.dirname(__file__), 'QTN')
    qtn_genes = {}
    if os.path.exists(file_path):
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                for line in lines[1:]:  # 跳过表头
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        gene_name = parts[0].strip()
                        chr_val = parts[1].strip().replace('chr', '')  # chr1 -> 1
                        pos_bp = parts[2].strip()
                        try:
                            pos_mb = float(pos_bp) / 1_000_000
                            if chr_val not in qtn_genes:
                                qtn_genes[chr_val] = []
                            qtn_genes[chr_val].append((gene_name, pos_mb))
                        except ValueError:
                            continue
        except Exception as e:
            st.error(f"读取 QTN 文件时出错: {e}")
    return qtn_genes

qtn_genes = load_qtn_genes()

# --- 全局配置 ---
GROUP_NAME_MAP = {
    0: 'Indica',
    1: 'Temperate japonica',
    2: 'Tropical japonica',
    3: 'Basmati',
    4: 'Aus',
    5: 'O.rufipogon'
}

GROUP_COLOR_MAP = {
    0: '#1f77b4',
    1: '#ff7f0e',
    2: '#2ca02c',
    3: '#d62728',
    4: '#9467bd',
    5: '#8c564b'
}

TRAIT_COLOR_MAP = {
    '产量组成相关': "#d62728",
    '植株形态': '#ff7f0e',
    '抽穗期': "#2ca02c",
    '生物胁迫': '#1f77b4',
    '非生物胁迫': '#9467bd',
    '口感品质': "#8c564b",
    '种子形态': '#e377c2',
    '次生代谢相关': "#17becf",
    '其他': '#bcbd22'
}

EVALUATION_COLOR_MAP = {
    'Favorable': 'red',
    'Unfavorable': 'blue',
    'Neutral': 'black',
    'Unknown': 'gray'
}

# --- 辅助函数 ---

def get_text(key, lang, section='UI', format_args=None):
    """获取翻译文本"""
    lang_code = 'zh' if lang == '中文' else 'en'
    text = translations.get(section, {}).get(key, {}).get(lang_code, key)
    return text.format(*format_args) if format_args else text

def load_data(uploaded_file):
    """读取并清洗数据"""
    try:
        df = pd.read_csv(
            uploaded_file, 
            sep='\t', 
            skiprows=0, 
            comment='[', 
            header=0
        )
        
        # 清洗列名
        df.columns = df.columns.str.strip()

        # 获取原始文件的最后一列作为样本列，并重命名为 'genetype'
        original_sample_col = df.columns[-1]
        df.rename(columns={original_sample_col: 'genetype'}, inplace=True)
        sample_col = 'genetype'
        
        # 必须包含的列
        required_cols = ['Chr', 'Pos_7.0', 'Alt_Allele_Func', 'GeneName']
        if not all(col in df.columns for col in required_cols):
            st.error(f"上传的文件缺少必要的列: {required_cols}")
            return None, None, None

        # --- 基础数据清洗 ---
        df['Chr'] = df['Chr'].astype(str).str.strip()
        df['Chr_Clean'] = df['Chr'].str.replace('Chr', '', case=False).str.strip()
        df['Pos_7.0'] = df['Pos_7.0'].astype(str).str.strip()
        df['Pos_7.0_Clean'] = df['Pos_7.0'].apply(lambda x: x.split('-')[0])
        df['Pos_7.0_Num'] = pd.to_numeric(df['Pos_7.0_Clean'], errors='coerce').astype('Int64')
        df['Alt_Allele_Func'] = df['Alt_Allele_Func'].fillna('Unknown Function')
        df['GeneName'] = df['GeneName'].astype(str).str.strip().fillna('N/A')
        df['genetype'] = df['genetype'].astype(str).str.strip()
        df['Pos_Mb'] = df['Pos_7.0_Num'] / 1_000_000

        # --- 向量化数据匹配 ---
        key_simple = df['Chr_Clean'] + '_' + df['Pos_7.0']
        key_chr = 'Chr' + df['Chr_Clean'] + '_' + df['Pos_7.0']
        key_name = df['GeneName']
        
        def lookup_attribute(attr_name, default_val=None):
            res1 = key_simple.map(lambda k: gene_db.get(k))
            res2 = key_chr.map(lambda k: gene_db.get(k))
            res3 = key_name.map(lambda k: gene_db.get(k))
            final_obj = res1.combine_first(res2).combine_first(res3)
            
            if attr_name == 'whole_object':
                return final_obj
            return final_obj.map(lambda x: x.get(attr_name) if isinstance(x, dict) else default_val)

        info_series = lookup_attribute('whole_object')
        df['Trait_Category'] = info_series.map(lambda x: x.get('Trait') if isinstance(x, dict) else None).fillna('其他')
        df['RAP_Locus'] = info_series.map(lambda x: x.get('RAP_Locus', 'N/A') if isinstance(x, dict) else 'N/A')
        df['MSU_Locus'] = info_series.map(lambda x: x.get('MSU_Locus', 'N/A') if isinstance(x, dict) else 'N/A')
        
        def evaluate_genotype(row_tuple):
            info, g_type = row_tuple
            if not isinstance(info, dict) or not info.get('Evaluation'):
                return "Unknown"
            if g_type in ["|", ".|."] or g_type.startswith("DEL|"):
                return "Unknown"
            
            ref_status = info['Evaluation']
            if ref_status == "Neutral" or g_type == '0|0':
                return ref_status
            return "Unfavorable" if ref_status == "Favorable" else "Favorable"

        df['Evaluation'] = list(map(evaluate_genotype, zip(info_series, df['genetype'])))
        
        # --- 染色体与位置处理 ---

        df['Is_Valid_Chr'] = df['Chr_Clean'].apply(lambda c: c.isdigit() and 1 <= int(c) <= 12)
        df['Plot_Chr'] = df['Chr_Clean']
        
        mask_need_random = ~df['Is_Valid_Chr'] | (df['Is_Valid_Chr'] & df['Pos_Mb'].isna())
        df.loc[mask_need_random, 'Plot_Chr'] = 'Other Gene'
        
        if mask_need_random.any():
            def get_random_pos(name):
                random.seed(str(name))
                return random.uniform(0, 50)
            df.loc[mask_need_random, 'Pos_Mb'] = df.loc[mask_need_random, 'GeneName'].apply(get_random_pos)
            
        df['Is_Plot_Valid'] = df['Pos_Mb'].notna()
        unique_chrs = df[df['Is_Plot_Valid']]['Plot_Chr'].unique()
        numeric_chrs = sorted([int(c) for c in unique_chrs if str(c).isdigit()])
        other_chrs = sorted([c for c in unique_chrs if not str(c).isdigit()])
        valid_chromosomes = [str(c) for c in numeric_chrs] + other_chrs
        
        chr_y_map = {str(chr_name): i + 1 for i, chr_name in enumerate(valid_chromosomes)}
        df['Y_pos'] = df['Plot_Chr'].astype(str).map(chr_y_map)
        df['Chr_Label'] = df['Chr_Clean'].apply(lambda x: f"Chr{x}" if str(x).isdigit() else str(x))

        return df, chr_y_map, sample_col

    except Exception as e:
        st.error(f"读取文件时出错: {e}")
        return None, None, None


@st.cache_data
def load_segments_data(uploaded_file):
    """读取祖源分段文件（tsv），格式示例：Chr1 1131 19328410 hap2 0"""
    try:
        seg_df = pd.read_csv(
            uploaded_file,
            sep='\t',
            header=None,
            comment='[',
            dtype=str
        )
        if seg_df.shape[1] < 5:
            st.error("Segments TSV 列数不足，期望至少 5 列：Chr, Start, End, Hap, Group")
            return None

        seg_df = seg_df.iloc[:, :5]
        seg_df.columns = ['Chr', 'Start', 'End', 'Haplotype', 'Group']

        seg_df['Chr'] = seg_df['Chr'].astype(str).str.strip()
        seg_df['Chr_Clean'] = seg_df['Chr'].str.replace('Chr', '', case=False).str.strip()
        seg_df['Start'] = pd.to_numeric(seg_df['Start'], errors='coerce')
        seg_df['End'] = pd.to_numeric(seg_df['End'], errors='coerce')
        seg_df['Group'] = pd.to_numeric(seg_df['Group'], errors='coerce').astype('Int64')
        seg_df['Haplotype'] = seg_df['Haplotype'].astype(str).str.strip()

        seg_df = seg_df.dropna(subset=['Start', 'End', 'Group', 'Chr_Clean', 'Haplotype']).copy()
        seg_df = seg_df[seg_df['Start'] <= seg_df['End']]
        seg_df['Start_Mb'] = seg_df['Start'] / 1_000_000
        seg_df['End_Mb'] = seg_df['End'] / 1_000_000

        return seg_df
    except Exception as e:
        st.error(f"读取 segments 文件时出错: {e}")
        return None

def display_ancestry_stats(stats: dict, lang_code: str):
    """展示祖源统计数据"""
    col_hap1, col_hap2 = st.columns(2)
    
    for col, hap_name in [(col_hap1, 'hap1'), (col_hap2, 'hap2')]:
        with col:
            st.markdown(f"**{hap_name.capitalize()}**")
            if hap_name in stats:
                hap_data = stats[hap_name]
                for group_name in sorted(hap_data.keys()):
                    if group_name == 'total_length':
                        continue
                    group_info = hap_data[group_name]
                    # 使用斜体显示祖源名称
                    st.markdown(
                        f"**_{group_name}_**: {group_info['percentage']:.1f}% "
                        f"({group_info['length']/1_000_000:.2f} Mb)"
                    )

def plot_segments_for_chr_both_haps(seg_df: pd.DataFrame, chr_name: str, lang_code: str):
    """对单条染色体绘制 hap1 + hap2，并返回统计数据"""
    if seg_df is None or seg_df.empty:
        return None, None

    df_chr = seg_df[seg_df['Chr_Clean'].astype(str) == str(chr_name)].copy()
    if df_chr.empty:
        return None, None

    hap_y = {'hap1': 1, 'hap2': 0}
    fig = go.Figure()
    present_haps = [h for h in ['hap1', 'hap2'] if (df_chr['Haplotype'] == h).any()]
    if not present_haps:
        return None, None

    group_codes = sorted(df_chr['Group'].dropna().unique().tolist())
    shown_legend_groups = set()
    stats_dict = {}
    
    for hap in present_haps:
        df_h = df_chr[df_chr['Haplotype'] == hap]
        y0 = hap_y.get(hap, 0)
        hap_stats = {}
        total_length = 0
        
        for group_code in group_codes:
            if pd.isna(group_code):
                continue
            group_code_int = int(group_code)
            sub = df_h[df_h['Group'] == group_code]
            if sub.empty:
                continue
            
            group_length = (sub['End'] - sub['Start']).sum()
            total_length += group_length
            hap_stats[GROUP_NAME_MAP.get(group_code_int, str(group_code_int))] = group_length

            show_in_legend = group_code_int not in shown_legend_groups
            if show_in_legend:
                shown_legend_groups.add(group_code_int)

            x_vals, y_vals = [], []
            for _, r in sub.iterrows():
                x_vals.extend([r['Start_Mb'], r['End_Mb'], None])
                y_vals.extend([y0, y0, None])

            # 使用斜体格式的名称
            group_name_italic = f"<i>{GROUP_NAME_MAP.get(group_code_int, str(group_code_int))}</i>"
            
            fig.add_trace(
                go.Scatter(
                    x=x_vals, y=y_vals, mode='lines',
                    name=group_name_italic,
                    legendgroup=str(group_code_int), showlegend=show_in_legend,
                    line=dict(color=GROUP_COLOR_MAP.get(group_code_int, '#7f7f7f'), width=12),
                    hovertemplate=f"Group=%{{fullData.name}}<br>Haplotype={hap}<br>Position=%{{x:.3f}} Mb<extra></extra>"
                )
            )
        
        hap_stats_with_pct = {}
        for group_name, length in hap_stats.items():
            percentage = (length / total_length * 100) if total_length > 0 else 0
            hap_stats_with_pct[group_name] = {'length': length, 'percentage': percentage}
        hap_stats_with_pct['total_length'] = total_length
        stats_dict[hap] = hap_stats_with_pct

    # 先计算染色体范围（用于后续布局计算）
    title = f"Chr{chr_name}" if str(chr_name).isdigit() else str(chr_name)
    x_max, x_min = float(df_chr['End_Mb'].max()), float(df_chr['Start_Mb'].min())
    pad = max(1.0, (x_max - x_min) * 0.02)
    
    # 添加 hap1 的 QTN 基因标记
    if str(chr_name) in qtn_genes:
        genes_on_chr = sorted(qtn_genes[str(chr_name)], key=lambda x: x[1])  # 按位置排序
        
        line_color = '#c0392b'  # 统一红色
        label_y = 1.28  # 统一标签高度
        min_label_spacing = 1.0  # 标签最小间距(Mb)，字体加大后增加间距
        
        # 计算每个基因标签的x位置，确保不重叠
        label_positions = []
        for idx, (gene_name, pos_mb) in enumerate(genes_on_chr):
            if idx == 0:
                label_x = pos_mb
            else:
                prev_label_x = label_positions[-1]
                # 如果实际位置间距够大就用实际位置，否则往右推
                if pos_mb >= prev_label_x + min_label_spacing:
                    label_x = pos_mb
                else:
                    label_x = prev_label_x + min_label_spacing
            label_positions.append(label_x)
        
        # 计算标签需要的额外空间
        if label_positions:
            max_label_x = max(label_positions)
            extra_pad = max(0, max_label_x - x_max + 2)  # 额外padding
        else:
            extra_pad = 0
        
        # 绘制所有基因标记
        for idx, (gene_name, pos_mb) in enumerate(genes_on_chr):
            label_x = label_positions[idx]
            
            # 判断是否需要斜线（标签位置偏离实际位置）
            if abs(label_x - pos_mb) < 0.1:
                # 几乎垂直，画垂直线
                fig.add_trace(
                    go.Scatter(
                        x=[pos_mb, pos_mb],
                        y=[1.08, label_y],
                        mode='lines',
                        line=dict(color=line_color, width=1),
                        showlegend=False,
                        hovertemplate=f"<b>{gene_name}</b><br>Position: {pos_mb:.3f} Mb<extra></extra>"
                    )
                )
            else:
                # 斜线连接
                fig.add_trace(
                    go.Scatter(
                        x=[pos_mb, label_x],
                        y=[1.08, label_y],
                        mode='lines',
                        line=dict(color=line_color, width=1),
                        showlegend=False,
                        hovertemplate=f"<b>{gene_name}</b><br>Position: {pos_mb:.3f} Mb<extra></extra>"
                    )
                )
            
            # 添加基因名（斜体）
            fig.add_annotation(
                x=label_x, y=label_y,
                text=f"<i>{gene_name}</i>",
                showarrow=False,
                font=dict(size=10, color='#2c3e50'),
                textangle=-45,
                xanchor='left',
                yanchor='bottom'
            )
    else:
        extra_pad = 0

    fig.update_layout(
        title=title, height=260, margin=dict(l=10, r=10, t=70, b=30), plot_bgcolor='white',
        xaxis=dict(
            title='Position (Mb)' if lang_code != 'zh' else '位置 (Mb)',
            range=[max(0, x_min - pad), x_max + pad + extra_pad],
            showgrid=True, gridcolor='rgba(0,0,0,0.06)', zeroline=False
        ),
        yaxis=dict(
            title='', tickmode='array', tickvals=[0, 1], ticktext=['hap2', 'hap1'],
            range=[-0.6, 1.65], showgrid=False, zeroline=False
        ),
        legend=dict(orientation='h', yanchor='bottom', y=1.28, xanchor='left', x=0)    # 图标注高度
    )
    
    return fig, stats_dict

# --- 侧边栏: 语言选择 ---
st.sidebar.title("🌐 Language")

if 'lang_code' not in st.session_state:
    st.session_state.lang_code = st.session_state.get('language', '中文') == '中文' and 'zh' or 'en'

lang_label = "选择语言" if st.session_state.lang_code == 'zh' else "Select Language"
format_func = lambda x: {"zh": "中文", "en": "英语" if st.session_state.lang_code == 'zh' else "English"}[x]

selected_code = st.sidebar.radio(
    lang_label, ['zh', 'en'],
    index=0 if st.session_state.lang_code == 'zh' else 1,
    format_func=format_func, key="lang_radio_code"
)

if selected_code != st.session_state.lang_code:
    st.session_state.lang_code = selected_code
    st.rerun()

lang = '中文' if st.session_state.lang_code == 'zh' else 'English'

# --- 侧边栏: 数据上传与分类选择 ---

st.sidebar.title(get_text("sidebar_title", lang))

# 1. 文件上传部分
st.sidebar.header(get_text("data_source_header", lang))
data_source_options = [get_text("use_demo", lang), get_text("upload_file_option", lang)]
# 为了逻辑判断方便，我们需要知道用户选的是第几个选项，而不是依赖文本
data_source_idx = st.sidebar.radio(
    get_text("select_data_label", lang), 
    range(len(data_source_options)), 
    format_func=lambda x: data_source_options[x]
)

df = None
chr_y_map = None
sample_col = None

segments_df = None

if data_source_idx == 0: # 使用 Demo 数据
    demo_path = os.path.join(os.path.dirname(__file__), 'demo data', 'HHZ.geno')
    if os.path.exists(demo_path):
        with open(demo_path, 'r') as f:
            df, chr_y_map, sample_col = load_data(f)
        st.sidebar.success(get_text("load_success_demo", lang))
    else:
        st.sidebar.error(get_text("file_not_found", lang))

    # demo 模式下也默认加载 demo segments
    segments_demo_path = os.path.join(os.path.dirname(__file__), 'demo data', 'HHZ.segments.tsv')
    if os.path.exists(segments_demo_path):
        with open(segments_demo_path, 'r', encoding='utf-8') as f:
            segments_df = load_segments_data(f)
    else:
        segments_df = None
else: # 上传文件
    uploaded_file = st.sidebar.file_uploader(get_text("upload_label", lang), type=['geno', 'txt', 'csv', 'tsv'])
    if uploaded_file is not None:
        df, chr_y_map, sample_col = load_data(uploaded_file)
        st.sidebar.success(get_text("load_success_file", lang, format_args=[uploaded_file.name]))

    # segments 上传入口合并到 geno 上传下面
    segments_uploaded = st.sidebar.file_uploader(
        "上传 segments TSV（祖源分析）" if st.session_state.lang_code == 'zh' else "Upload segments TSV (ancestry)",
        type=['tsv', 'txt'],
        key="segments_uploader"
    )
    if segments_uploaded is not None:
        segments_df = load_segments_data(segments_uploaded)

# --- 搜索功能 ---
st.sidebar.markdown("---")
st.sidebar.header(get_text("search_header", lang, section="UI") if "search_header" in translations.get("UI", {}) else ("Search Gene" if st.session_state.lang_code != 'zh' else "搜索基因"))
search_gene_query = st.sidebar.text_input(
    "",
    placeholder="Gene ID / MSU ID / RAP ID" if st.session_state.lang_code != 'zh' else "基因ID / MSU ID / RAP ID",
    label_visibility="collapsed"
).strip()

# --- 共享染色体筛选：对 geno 与 segments 同时生效 ---
selected_chr = 'All'
shared_chr_options = []
if df is not None:
    shared_chr_options.extend(df['Plot_Chr'].dropna().astype(str).unique().tolist())
if segments_df is not None and not segments_df.empty:
    shared_chr_options.extend(segments_df['Chr_Clean'].dropna().astype(str).unique().tolist())

shared_chr_options = sorted(set(shared_chr_options), key=lambda x: int(x) if str(x).isdigit() else 999)
if shared_chr_options:
    st.sidebar.header(get_text("chr_filter_header", lang))
    all_chr_label = get_text("all_chr_option", lang)
    chr_options = [all_chr_label] + shared_chr_options
    selected_chr_display = st.sidebar.selectbox(get_text("select_chr_label", lang), chr_options, key="shared_chr_select")
    if selected_chr_display != all_chr_label:
        selected_chr = selected_chr_display

# --- 主界面逻辑 ---

if df is not None:
    # 3. 可视化模式选择
    st.sidebar.header(get_text("view_mode_header", lang))
    view_mode_options = [get_text("mode_trait", lang), get_text("mode_eval", lang)]
    view_mode_idx = st.sidebar.radio(
        get_text("view_mode_label", lang), 
        range(len(view_mode_options)), 
        format_func=lambda x: view_mode_options[x]
    )

    # 4. 分类筛选 (模拟左侧大分类)
    st.sidebar.header(get_text("filter_header", lang))
    
    # 获取所有可用分类 (原始中文 Key)
    raw_categories = sorted(df['Trait_Category'].unique().tolist())
    
    # 构建显示用的分类列表 (翻译后)
    # 使用 'All' 的翻译作为第一个选项
    all_option_label = get_text("all_option", lang)
    
    # 创建一个映射：显示文本 -> 原始中文 Key
    display_to_raw_cat = {all_option_label: 'All'}
    display_categories = [all_option_label]
    
    for cat in raw_categories:
        # 尝试从 Categories 部分获取翻译，如果没有则显示原文
        trans_cat = get_text(cat, lang, section='Categories')
        display_categories.append(trans_cat)
        display_to_raw_cat[trans_cat] = cat
    
    selected_display_category = st.sidebar.radio(get_text("filter_label", lang), display_categories)
    selected_raw_category = display_to_raw_cat[selected_display_category]

    # 根据选择筛选数据
    if selected_raw_category != 'All':
        filtered_df = df[df['Trait_Category'] == selected_raw_category].copy()
        # 标题也要翻译分类名称
        cat_trans = get_text(selected_raw_category, lang, section='Categories')
        main_title = get_text("main_title_filtered", lang, format_args=[cat_trans])
    else:
        filtered_df = df.copy()
        main_title = get_text("main_title_all", lang)

    # 再按染色体筛选
    if selected_chr != 'All':
        filtered_df = filtered_df[filtered_df['Plot_Chr'].astype(str) == selected_chr]
        main_title += f" - {selected_chr}"

    # --- 主区域 ---
    
    st.title(main_title)
    
    # 统计仪表盘
    total_genes = len(filtered_df)
    favorable = len(filtered_df[filtered_df['Evaluation'] == 'Favorable'])
    unfavorable = len(filtered_df[filtered_df['Evaluation'] == 'Unfavorable'])
    neutral_unknown = total_genes - favorable - unfavorable
    
    # 图片路径
    img_dir = os.path.join(os.path.dirname(__file__), 'picture')
    
    # 创建4列布局，增加间距
    col1, col2, col3, col4 = st.columns([1.2, 1.2, 1.2, 1.2], gap="medium")
    
    # 第1列：基因总数
    with col1:
        subcol1, subcol2 = st.columns([1.2, 1.5])
        with subcol1:
            total_img_path = os.path.join(img_dir, 'total.png')
            if os.path.exists(total_img_path):
                st.image(total_img_path, use_container_width=True)
        with subcol2:
            st.markdown(f"<p style='font-size:14px; color:#666; margin-bottom:0;'> {'Total Genes' if st.session_state.lang_code != 'zh' else '基因总数'}</p>", unsafe_allow_html=True)
            st.markdown(f"<h1 style='margin-top:0; margin-bottom:10px;'>{total_genes}</h1>", unsafe_allow_html=True)
    
    # 第2列：有利基因
    with col2:
        subcol1, subcol2 = st.columns([1.2, 1.5])
        with subcol1:
            favorable_img_path = os.path.join(img_dir, 'Favorable.png')
            if os.path.exists(favorable_img_path):
                st.image(favorable_img_path, use_container_width=True)
        with subcol2:
            st.markdown(f"<p style='font-size:14px; color:#666; margin-bottom:0;'> {'Favorable' if st.session_state.lang_code != 'zh' else '有利基因'}</p>", unsafe_allow_html=True)
            st.markdown(f"<h1 style='margin-top:0; margin-bottom:5px;'>{favorable}</h1>", unsafe_allow_html=True)
            if total_genes > 0:
                st.markdown(f"<div style='background-color:#d4edda; color:#28a745; padding:4px 12px; border-radius:15px; display:inline-block; font-size:14px; font-weight:600;'>↑ {favorable/total_genes*100:.1f}%</div>", unsafe_allow_html=True)
    
    # 第3列：不利基因
    with col3:
        subcol1, subcol2 = st.columns([1.2, 1.5])
        with subcol1:
            unfavorable_img_path = os.path.join(img_dir, 'UnFavorable.png')
            if os.path.exists(unfavorable_img_path):
                st.image(unfavorable_img_path, use_container_width=True)
        with subcol2:
            st.markdown(f"<p style='font-size:14px; color:#666; margin-bottom:0;'> {'Unfavorable' if st.session_state.lang_code != 'zh' else '不利基因'}</p>", unsafe_allow_html=True)
            st.markdown(f"<h1 style='margin-top:0; margin-bottom:5px;'>{unfavorable}</h1>", unsafe_allow_html=True)
            if total_genes > 0:
                st.markdown(f"<div style='background-color:#f8d7da; color:#dc3545; padding:4px 12px; border-radius:15px; display:inline-block; font-size:14px; font-weight:600;'>↓ {unfavorable/total_genes*100:.1f}%</div>", unsafe_allow_html=True)
    
    # 第4列：中性/未知
    with col4:
        subcol1, subcol2 = st.columns([1.2, 1.5])
        with subcol1:
            natural_img_path = os.path.join(img_dir, 'Natural.png')
            if os.path.exists(natural_img_path):
                st.image(natural_img_path, use_container_width=True)
        with subcol2:
            st.markdown(f"<p style='font-size:14px; color:#666; margin-bottom:0;'> {'Neutral/Unknown' if st.session_state.lang_code != 'zh' else '中性/未知'}</p>", unsafe_allow_html=True)
            st.markdown(f"<h1 style='margin-top:0; margin-bottom:5px;'>{neutral_unknown}</h1>", unsafe_allow_html=True)
            if total_genes > 0:
                st.markdown(f"<div style='background-color:#e7e7f5; color:#6c757d; padding:4px 12px; border-radius:15px; display:inline-block; font-size:14px; font-weight:600;'>{neutral_unknown/total_genes*100:.1f}%</div>", unsafe_allow_html=True)
    
    # AI智能分析按钮
    st.markdown("---")
    col_ai_btn = st.columns([2, 1, 2])
    with col_ai_btn[1]:
        if st.button("🤖 AI 智能分析" if st.session_state.lang_code == 'zh' else "🤖 AI Smart Analysis", 
                     use_container_width=True, type="primary"):
            with st.spinner("🧠 AI分析中..." if st.session_state.lang_code == 'zh' else "🧠 AI Analyzing..."):
                analysis_result = generate_ai_analysis(
                    filtered_df, total_genes, favorable, unfavorable, neutral_unknown,
                    selected_raw_category, segments_df, st.session_state.lang_code
                )
                st.session_state.ai_analysis_result = analysis_result
    
    # 显示AI分析结果
    if 'ai_analysis_result' in st.session_state and st.session_state.ai_analysis_result:
        st.markdown("---")
        
        # 创建一个带边框的容器显示结果
        result_container = st.container()
        with result_container:
            col_title, col_download = st.columns([3, 1])
            with col_title:
                st.markdown("### 📊 " + ("AI 分析报告" if st.session_state.lang_code == 'zh' else "AI Analysis Report"))
            with col_download:
                st.download_button(
                    label="📥 " + ("下载报告" if st.session_state.lang_code == 'zh' else "Download"),
                    data=st.session_state.ai_analysis_result,
                    file_name=f"breeding_analysis_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.md",
                    mime="text/markdown",
                    use_container_width=True
                )
            
            # 显示分析内容
            st.markdown(st.session_state.ai_analysis_result)
            
            # 清除按钮
            if st.button("🗑️ " + ("清除报告" if st.session_state.lang_code == 'zh' else "Clear Report")):
                del st.session_state.ai_analysis_result
                st.rerun()
    
    # 快速筛选提示
    if total_genes > 0:
        st.caption(
            f"💡 Tip: Use the sidebar filters to focus on specific categories or chromosomes" 
            if st.session_state.lang_code != 'zh' 
            else f"💡 提示：使用侧边栏筛选器可聚焦于特定分类或染色体"
        )
    
    st.markdown("---")
    
    # --- 交互式绘图模块 ---
    
    # 准备绘图数据 (仅包含有效坐标的行)
    plot_df = filtered_df[filtered_df['Is_Plot_Valid']].copy()
    
    # 使用原始 df 计算染色体长度和背景，确保所有染色体都显示
    all_plot_df = df[df['Is_Plot_Valid']].copy()
    if not all_plot_df.empty:
        chr_lengths = all_plot_df.groupby('Plot_Chr')['Pos_Mb'].max().reset_index()
        chr_lengths['Max_Pos'] = chr_lengths['Pos_Mb'] + 2
        chr_lengths['Sort_Key'] = chr_lengths['Plot_Chr'].apply(lambda v: int(v) if str(v).isdigit() else 999)
        chr_lengths = chr_lengths.sort_values('Sort_Key')
        chr_lengths['Y_pos'] = chr_lengths['Plot_Chr'].map(lambda x: chr_y_map.get(str(x)))
        X_MAX = chr_lengths['Max_Pos'].max()
    else:
        X_MAX = 50
        chr_lengths = pd.DataFrame()

    plot_df['row_id'] = plot_df.index

    plot_df['Alt_Allele_Func_Display'] = plot_df['Alt_Allele_Func'].map(lambda x: get_text(x, lang, section='Functions'))
    plot_df['Evaluation_Display'] = plot_df['Evaluation'].map(lambda x: get_text(x, lang, section='Evaluation'))
    plot_df['Trait_Category_Display'] = plot_df['Trait_Category'].map(lambda x: get_text(x, lang, section='Categories'))

    if view_mode_idx == 0:
        color_map = TRAIT_COLOR_MAP
        plot_color_map = {get_text(k, lang, section='Categories'): v for k, v in color_map.items()}
        plot_color_col = 'Trait_Category_Display'
    else:
        color_map = EVALUATION_COLOR_MAP
        plot_color_map = {get_text(k, lang, section='Evaluation'): v for k, v in color_map.items()}
        plot_color_col = 'Evaluation_Display'

    # 绘制散点图
    if not plot_df.empty:
        # 大小说明提示
        st.caption(
            "💡 Marker size indicates QTN importance: ● Large (Major effect) ● Medium (Minor effect) ● Small (Other loci)"
            if st.session_state.lang_code != 'zh'
            else "💡 标记大小表示位点重要性：● 大（主效位点）● 中（次要位点）● 小（其他位点）"
        )
        fig = px.scatter(
            plot_df, 
            x='Pos_Mb', 
            y='Y_pos', 
            color=plot_color_col, # 使用翻译过的列
            color_discrete_map=plot_color_map, # 使用翻译过的颜色映射
            hover_name="GeneName", 
            custom_data=['row_id'], # 将 row_id 传入 custom_data
            hover_data={
                'Chr_Label': True,          
                'Pos_Mb': False, # 不显示计算出的 Mb 位置           
                'Pos_7.0': True, # 显示原始位置信息 (可能是数字也可能是字符串)          
                'GeneName': False,          
                'Alt_Allele_Func': False, # 隐藏原始列
                'Alt_Allele_Func_Display': True, # 显示翻译后的列
                sample_col: True, # 显示样本基因型
                'Evaluation': False, # 隐藏原始评估
                'Evaluation_Display': True, # 始终显示翻译后的评估结果
                'Y_pos': False,             
                'Chr': False, 
                'Chr_Clean': False,
                'row_id': False,
                'Trait_Category': False, # 隐藏原始中文列
                'Trait_Category_Display': True # 始终显示翻译后的性状分类
            },
            labels={
                'Pos_Mb': get_text("col_position", lang) + ' (Mb)', 
                'Pos_7.0': get_text("col_position", lang), # 使用翻译后的 "位置"
                'Trait_Category': get_text("Trait_Category", lang, section='Columns'),
                'Trait_Category_Display': get_text("Trait_Category", lang, section='Columns'),
                'Alt_Allele_Func': get_text("Alt_Allele_Func", lang, section='Columns'),
                'Alt_Allele_Func_Display': get_text("Alt_Allele_Func", lang, section='Columns'),
                sample_col: get_text("genetype", lang, section='Columns'),
                'Evaluation': get_text("Evaluation", lang, section='Columns'),
                'Evaluation_Display': get_text("Evaluation", lang, section='Columns'),
                'Chr_Label': get_text("Chr", lang, section='Columns') # 将 Chr_Label 映射为 "染色体" 或 "Chr"
            },
            height=600
        )

        # 添加染色体背景线
        shapes = []
        if not chr_lengths.empty:
            for index, row in chr_lengths.iterrows():
                if pd.notna(row['Y_pos']):
                    shapes.append(
                        go.layout.Shape(
                            type="line",
                            x0=0, y0=row['Y_pos'], 
                            x1=row['Max_Pos'], # 使用动态计算的长度
                            y1=row['Y_pos'],
                            line=dict(color="lightgray", width=10),
                            layer="below"
                        )
                    )

        fig.update_traces(marker=dict(opacity=0.8, line=dict(width=1, color='white')))
        
        for trace in fig.data:
            if trace.customdata is not None:
                sizes = []
                for pt_data in trace.customdata:
                    row_id = pt_data[0]
                    if row_id in plot_df.index:
                        key = f"{plot_df.loc[row_id, 'Chr']}_{str(plot_df.loc[row_id, 'Pos_7.0']).strip()}"
                        effect = qtn_effects.get(key, None)
                        sizes.append(18 if effect == 'high' else 13 if effect == 'low' else 8)
                    else:
                        sizes.append(6)
                trace.marker.size = sizes
        if search_gene_query:
            # 支持GeneName、MSU_Locus、RAP_Locus三种方式搜索
            search_mask = (
                plot_df['GeneName'].astype(str).str.contains(search_gene_query, case=False, na=False) |
                plot_df['MSU_Locus'].astype(str).str.contains(search_gene_query, case=False, na=False) |
                plot_df['RAP_Locus'].astype(str).str.contains(search_gene_query, case=False, na=False)
            )
            highlight_df = plot_df[search_mask]
            
            if not highlight_df.empty:
                fig.add_trace(
                    go.Scatter(
                        x=highlight_df['Pos_Mb'], y=highlight_df['Y_pos'], mode='markers',
                        marker=dict(size=20, color='red', symbol='circle-open', line=dict(width=3, color='red')),
                        name='Search Result', hoverinfo='skip', showlegend=True
                    )
                )
                msg = f"Found {len(highlight_df)} genes matching '{search_gene_query}'"
                st.toast(msg if st.session_state.lang_code != 'zh' else f"找到 {len(highlight_df)} 个匹配 '{search_gene_query}' 的基因")
            else:
                msg = f"No genes found matching '{search_gene_query}'"
                st.toast(msg if st.session_state.lang_code != 'zh' else f"未找到匹配 '{search_gene_query}' 的基因")
        # 使用 chr_lengths 中的所有染色体来确定Y轴刻度
        present_y_pos = set()
        if not chr_lengths.empty:
            present_y_pos.update(chr_lengths['Y_pos'].dropna().unique())
        # 如果有筛选后的数据点，也包含它们
        if not plot_df.empty:
            present_y_pos.update(plot_df['Y_pos'].dropna().unique())
             
        y_ticks = sorted([item for item in chr_y_map.items() if item[1] in present_y_pos], key=lambda x: x[1])
        tick_vals = [item[1] for item in y_ticks]
        tick_text = [f'Chr{name}' if name.isdigit() else str(name) for name, _ in y_ticks]
        
        fig.update_layout(
            shapes=shapes,
            yaxis=dict(
                title="",
                tickvals=tick_vals,
                ticktext=tick_text,
                autorange="reversed", 
                showgrid=False,
                zeroline=False
            ),
            xaxis=dict(title=get_text("col_position", lang) + " (Mb)", range=[0, X_MAX]),
            plot_bgcolor='white',
            legend=dict(orientation="h", yanchor="bottom", y=0.98, xanchor="right", x=1),
            clickmode='event+select', # 允许点击选择
            hoverlabel=dict(
                font=dict(size=16)
            )
        )

        # 显示图表
        selection = st.plotly_chart(fig, use_container_width=True, on_select="rerun")
    else:
        st.warning(get_text("no_plot_data", lang))
        selection = None
    
    # --- 4. 详情与链接模块 ---
    
    st.markdown("---")
    # 列表标题翻译
    cat_display = get_text(selected_raw_category, lang, section='Categories') if selected_raw_category != 'All' else get_text("all_option", lang)
    st.subheader(get_text("list_header", lang, format_args=[cat_display]))
    st.caption(get_text("list_caption", lang))

    # 准备显示的数据表
    display_cols = ['Chr', 'Pos_7.0', 'GeneName', 'Alt_Allele_Func', 'Trait_Category', sample_col, 'Evaluation', 'RAP_Locus', 'MSU_Locus']
    
    try:
        if selection and 'selection' in selection and 'points' in selection['selection']:
            selected_row_ids = [p['customdata'][0] for p in selection['selection']['points'] if 'customdata' in p]
            if selected_row_ids:
                selected_data = filtered_df.loc[selected_row_ids]
                st.info(get_text("selected_info", lang, format_args=[len(selected_data)]))
                display_df = selected_data[display_cols].copy()
            else:
                display_df = filtered_df[display_cols].copy()
        else:
            display_df = filtered_df[display_cols].copy()
    except:
        display_df = filtered_df[display_cols].copy()

    for col, section in [('Trait_Category', 'Categories'), ('Evaluation', 'Evaluation'), ('Alt_Allele_Func', 'Functions')]:
        display_df[col] = display_df[col].map(lambda x: get_text(x, lang, section=section))
    
    if search_gene_query:
        # 支持GeneName、MSU_Locus、RAP_Locus三种方式搜索
        search_mask = (
            display_df['GeneName'].astype(str).str.contains(search_gene_query, case=False, na=False) |
            display_df['MSU_Locus'].astype(str).str.contains(search_gene_query, case=False, na=False) |
            display_df['RAP_Locus'].astype(str).str.contains(search_gene_query, case=False, na=False)
        )
        if search_mask.any():
            display_df = display_df[search_mask].copy()
            msg = f"🔍 Table filtered to show {len(display_df)} genes matching '{search_gene_query}'"
            st.info(msg if st.session_state.lang_code != 'zh' else f"🔍 表格已筛选，显示 {len(display_df)} 个匹配 '{search_gene_query}' 的基因")
    
    def to_url(val, base_url, param):
        s = str(val).strip()
        return f"{base_url}?{param}={s}" if s and s.upper() not in ['NA', 'N/A', ''] else None

    display_df['RAP_Locus'] = display_df['RAP_Locus'].apply(lambda v: to_url(v, 'https://rapdb.dna.affrc.go.jp/locus/', 'name'))
    display_df['MSU_Locus'] = display_df['MSU_Locus'].apply(lambda v: to_url(v, 'https://rice.uga.edu/cgi-bin/ORF_infopage.cgi', 'orf'))

    # 重命名列
    col_map = {
        'Chr': get_text("Chr", lang, section='Columns'),
        'Pos_7.0': get_text("Pos_7.0", lang, section='Columns'),
        'GeneName': get_text("GeneName", lang, section='Columns'),
        'Alt_Allele_Func': get_text("Alt_Allele_Func", lang, section='Columns'),
        'Trait_Category': get_text("Trait_Category", lang, section='Columns'),
        sample_col: get_text("genetype", lang, section='Columns'),
        'Evaluation': get_text("Evaluation", lang, section='Columns'),
        'RAP_Locus': 'RAP Locus',
        'MSU_Locus': 'MSU Locus'
    }
    display_df = display_df.rename(columns=col_map)
    
    # 获取翻译后的列名用于配置
    pos_col = get_text("Pos_7.0", lang, section='Columns')

    # --- 下载按钮和表格信息 ---
    col_download, col_info = st.columns([1, 3])
    
    with col_download:
        csv_data = display_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 CSV",
            data=csv_data,
            file_name=f'gene_data_{selected_raw_category if selected_raw_category else "all"}.csv',
            mime='text/csv',
            use_container_width=True
        )
    
    with col_info:
        st.caption(
            f"📋 Showing {len(display_df)} rows • Click column headers to sort" 
            if st.session_state.lang_code != 'zh' 
            else f"📋 显示 {len(display_df)} 行数据 • 点击列标题可排序"
        )

    st.dataframe(
        display_df,
        column_config={
            pos_col: st.column_config.TextColumn(
                get_text("col_position", lang)
            ),
            'RAP Locus': st.column_config.LinkColumn(
                "RAP Locus",
                display_text=r"name=(.*)"
            ),
            'MSU Locus': st.column_config.LinkColumn(
                "MSU Locus",
                display_text=r"orf=(.*)"
            )
        },
        hide_index=True,
        use_container_width=True
    )

else:
    # 更友好的空状态显示
    st.markdown("### " + ("Welcome to Rice GWAS Interactive Visualization" if st.session_state.lang_code != 'zh' else "欢迎使用水稻 GWAS 交互式可视化系统"))
    
    col1, col2 = st.columns([2, 1])
    with col1:
        st.markdown(
            """
            **Getting Started:**
            
            1. 📁 Choose a data source from the sidebar
            2. 📊 Explore genome-wide visualization
            3. 🔍 Search for specific genes
            4. 📥 Download filtered results
            
            **Features:**
            - Interactive chromosome viewer
            - Gene annotation and evaluation
            - Ancestry analysis visualization
            """
            if st.session_state.lang_code != 'zh'
            else """
            **快速开始：**
            
            1. 📁 从侧边栏选择数据源
            2. 📊 探索全基因组可视化
            3. 🔍 搜索特定基因
            4. 📥 下载筛选结果
            
            **功能特点：**
            - 交互式染色体浏览器
            - 基因注释与评估
            - 祖源分析可视化
            """
        )
    with col2:
        st.info(
            "💡 **Pro Tip:** Start with the demo data to explore all features!"
            if st.session_state.lang_code != 'zh'
            else "💡 **专业提示：** 从演示数据开始探索所有功能！"
        )

# --- 祖源分析（segments TSV）图：与 geno 图/列表解耦，但共享染色体筛选 ---
st.markdown("---")
st.subheader("Ancestral analysis" if st.session_state.lang_code != 'zh' else "祖源分析")
st.caption("The uncolored regions represent segments with " \
            "insufficient confidence in ancestry determination and have not been assigned to any reference ancestry type."
            if st.session_state.lang_code != 'zh' else"未着色区域为祖源判定置信度不足的区段，未被分配至任何参考祖源类型。"
        )

if segments_df is None or segments_df.empty:
    st.info(
        "Upload a segments TSV (or use Demo) to visualize."
        if st.session_state.lang_code != 'zh'
        else "请在左侧上传 segments TSV（或使用 Demo）以进行可视化。"
    )
else:
    available_chr = sorted(
        segments_df['Chr_Clean'].dropna().astype(str).unique().tolist(),
        key=lambda x: int(x) if str(x).isdigit() else 999
    )

    chrs_to_show = available_chr if selected_chr == 'All' else [selected_chr]
    
    if len(chrs_to_show) > 1:
        # 使用 Tabs 展示多条染色体
        tabs = st.tabs([f"Chr {c}" for c in chrs_to_show])
        for i, chr_name in enumerate(chrs_to_show):
            with tabs[i]:
                fig, stats = plot_segments_for_chr_both_haps(segments_df, chr_name, st.session_state.lang_code)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
                    if stats:
                        st.markdown("**" + ("Ancestry Composition" if st.session_state.lang_code != 'zh' else "祖源组成比例") + "**")
                        display_ancestry_stats(stats, st.session_state.lang_code)
                else:
                    st.info(f"No segment data for Chr {chr_name}" if st.session_state.lang_code != 'zh' else f"Chr {chr_name} 无祖源分段数据")
    else:
        # 单条染色体直接展示
        chr_name = chrs_to_show[0]
        fig, stats = plot_segments_for_chr_both_haps(segments_df, chr_name, st.session_state.lang_code)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
            if stats:
                st.markdown("**" + ("Ancestry Composition" if st.session_state.lang_code != 'zh' else "祖源组成比例") + "**")
                display_ancestry_stats(stats, st.session_state.lang_code)
        else:
            st.info(f"No segment data for Chr {chr_name}" if st.session_state.lang_code != 'zh' else f"Chr {chr_name} 无祖源分段数据")

