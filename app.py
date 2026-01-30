import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import json
import random
from typing import Optional

# 设置页面配置
st.set_page_config(page_title="Rice GWAS Interactive Visualization", layout="wide")

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

# --- 辅助函数 ---

def get_text(key, lang, section='UI', format_args=None):
    """获取翻译文本"""
    # 映射语言代码
    lang_code = 'zh' if lang == '中文' else 'en'
    try:
        text = translations.get(section, {}).get(key, {}).get(lang_code, key)
        if format_args:
            return text.format(*format_args)
        return text
    except:
        return key

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

        # --- 1. 基础数据清洗 ---
        
        # 清理 Chr，保留原始列但转字符串
        df['Chr'] = df['Chr'].astype(str).str.strip()
        df['Chr_Clean'] = df['Chr'].str.replace('Chr', '', case=False).str.strip()
        
        # 处理 Pos_7.0
        df['Pos_7.0'] = df['Pos_7.0'].astype(str).str.strip()
        df['Pos_7.0_Clean'] = df['Pos_7.0'].apply(lambda x: x.split('-')[0]) # 这个操作比较轻量，可以保留 apply
        df['Pos_7.0_Num'] = pd.to_numeric(df['Pos_7.0_Clean'], errors='coerce').astype('Int64')
        
        df['Alt_Allele_Func'] = df['Alt_Allele_Func'].fillna('Unknown Function')
        df['GeneName'] = df['GeneName'].astype(str).str.strip().fillna('N/A')
        df['genetype'] = df['genetype'].astype(str).str.strip()
        
        # 计算 Pos_Mb
        df['Pos_Mb'] = df['Pos_7.0_Num'] / 1_000_000

        # --- 2. 向量化数据匹配 (替代 apply) ---

        # 构造用于查找的 Key 列表
        # 原始逻辑是尝试: "Chr_Pos", "ChrX_Pos", "GeneName"
        # 我们可以先构造几个 Series 作为 Key
        
        # Key 1: "1_1000"
        key_simple = df['Chr_Clean'] + '_' + df['Pos_7.0']
        # Key 2: "Chr1_1000" (如果 Chr_Clean 是数字)
        key_chr = 'Chr' + df['Chr_Clean'] + '_' + df['Pos_7.0']
        # Key 3: GeneName
        key_name = df['GeneName']
        
        # 为了加速查找，我们将 gene_db 的某个字段提取出来做成查找表 (Map)
        # 例如要查 'Trait'
        # 但 gene_db 混杂了各种 key。直接用 map(gene_db.get) 是可行的，因为 Python 字典查找很快。
        # 相比于 df.apply(axis=1) 的主要开销在于 Series 的构建和上下文切换。
        # 直接对 Series 使用 map 会更快。

        def lookup_attribute(attr_name, default_val=None):
            # 将 gene_db 里的特定属性提取出来做成小字典，可能更快，但直接查 gene_db 也行
            # 这里我们直接映射整个对象，然后提取属性
            
            # 使用列表推导式可能比 map 更快，因为 map 会处理 NaN
            # 但为了代码简洁，我们用 map。注意 gene_db.get 返回的是 dict 或 None
            
            # 优先顺序：Key1 -> Key2 -> Key3
            
            # 1. 尝试 Key 1
            res1 = key_simple.map(lambda k: gene_db.get(k))
            
            # 2. 尝试 Key 2 (用于填补 res1 为 None 的)
            # 只有当 res1 为空时才需要查，但向量化通常全部查完再 combine
            res2 = key_chr.map(lambda k: gene_db.get(k))
            
            # 3. 尝试 Key 3
            res3 = key_name.map(lambda k: gene_db.get(k))
            
            # 合并结果：res1 优先，res2 次之，res3 最后
            # combine_first: Update null elements with value in the same location in 'other'.
            final_obj = res1.combine_first(res2).combine_first(res3)
            
            # 提取属性
            if attr_name == 'whole_object':
                 return final_obj

            return final_obj.map(lambda x: x.get(attr_name) if isinstance(x, dict) else default_val)

        # 获取 Info 对象 (只做一次大合并)
        info_series = lookup_attribute('whole_object')
        
        # 1. Trait_Category
        df['Trait_Category'] = info_series.map(lambda x: x.get('Trait') if isinstance(x, dict) else None).fillna('其他')
        
        # 2. RAP_Locus & MSU_Locus
        df['RAP_Locus'] = info_series.map(lambda x: x.get('RAP_Locus', 'N/A') if isinstance(x, dict) else 'N/A')
        df['MSU_Locus'] = info_series.map(lambda x: x.get('MSU_Locus', 'N/A') if isinstance(x, dict) else 'N/A')

        # 3. Evaluation
        # 逻辑较复杂，需要 info 和 genetype 配合
        # evaluation logic:
        # if not info or not info['Evaluation']: "Unknown"
        # if genetype in ["|", ".|.", "DEL|..."]: "Unknown"
        # if eval == "Neutral": "Neutral"
        # if genetype == "0|0": eval
        # else: flip(eval)
        
        def vectorized_evaluate(row_tuple):
            # row_tuple = (info_dict, genetype_str)
            info, g_type = row_tuple
            if not isinstance(info, dict):
                return "Unknown"
            
            ref_status = info.get('Evaluation')
            if not ref_status:
                return "Unknown"
            
            if g_type in ["|", ".|."] or g_type.startswith("DEL|"):
                return "Unknown"
            
            if ref_status == "Neutral":
                return "Neutral"
            
            if g_type == '0|0':
                return ref_status
            else:
                if ref_status == "Favorable": return "Unfavorable"
                if ref_status == "Unfavorable": return "Favorable"
            return "Unknown"

        # 构造临时 tuple series 用于 map (虽然不是纯向量化，但避免了 apply(axis=1) 的全行扫描)
        df['Evaluation'] = list(map(vectorized_evaluate, zip(info_series, df['genetype'])))
        
        # --- 3. 染色体与位置处理 (绘图准备) ---

        def is_valid_chr(c):
            return c.isdigit() and 1 <= int(c) <= 12

        df['Is_Valid_Chr'] = df['Chr_Clean'].apply(is_valid_chr) # 只有12个值，apply没问题

        # 分配 Other Gene 的位置
        # 逻辑：有效染色体但 Pos 为空 -> 标记为 -1 (虚拟)
        # 逻辑：无效染色体 -> 标记为 -1 (虚拟) + Chr 改为 Other Gene
        
        # 1. 确定 Pos_Mb
        # 已经计算了 Pos_Mb，如果有 NaN，需要处理
        
        # 2. 确定 Plot_Chr
        # 默认 Plot_Chr = Chr_Clean
        df['Plot_Chr'] = df['Chr_Clean']
        
        # 再次明确哪些是 "Other"
        # 条件 A: Is_Valid_Chr is False
        # 条件 B: Is_Valid_Chr is True BUT Pos_Mb is NaN
        
        mask_invalid_chr = ~df['Is_Valid_Chr']
        mask_valid_no_pos = df['Is_Valid_Chr'] & df['Pos_Mb'].isna()
        
        mask_need_random = mask_invalid_chr | mask_valid_no_pos
        
        # 对这些行，设置 Plot_Chr 为 'Other Gene'
        df.loc[mask_need_random, 'Plot_Chr'] = 'Other Gene'
        
        # 对这些行，生成随机 Pos
        if mask_need_random.any():
            # 生成随机数
            # 为了确定性，使用 GeneName 的 hash 
            # hash 可能会随 python 会话改变，使用 zlib.adler32 或者简单的 seeded random
            # apply 在这里只针对需要处理的行，量级通常较小
            def get_random_pos(name):
                random.seed(str(name))
                return random.uniform(0, 50)
            
            df.loc[mask_need_random, 'Pos_Mb'] = df.loc[mask_need_random, 'GeneName'].apply(get_random_pos)
            
        df['Is_Plot_Valid'] = df['Pos_Mb'].notna()

        # Generate Y axis mapping
        unique_chrs = df[df['Is_Plot_Valid']]['Plot_Chr'].unique()
        numeric_chrs = sorted([int(c) for c in unique_chrs if str(c).isdigit()])
        other_chrs = sorted([c for c in unique_chrs if not str(c).isdigit()])
        
        valid_chromosomes = [str(c) for c in numeric_chrs] + other_chrs
        chr_y_map = {str(chr_name): i + 1 for i, chr_name in enumerate(valid_chromosomes)}
        
        # Map Y pos
        df['Y_pos'] = df['Plot_Chr'].astype(str).map(chr_y_map)
        
        # Chr Label
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

        # 丢弃无效行（例如带表头的情况）
        seg_df = seg_df.dropna(subset=['Start', 'End', 'Group', 'Chr_Clean', 'Haplotype']).copy()
        seg_df = seg_df[seg_df['Start'] <= seg_df['End']]

        # 转 Mb 便于可视化
        seg_df['Start_Mb'] = seg_df['Start'] / 1_000_000
        seg_df['End_Mb'] = seg_df['End'] / 1_000_000

        return seg_df
    except Exception as e:
        st.error(f"读取 segments 文件时出错: {e}")
        return None

def plot_segments_for_chr_both_haps(seg_df: pd.DataFrame, chr_name: str, lang_code: str):
    """对单条染色体绘制 hap1 + hap2（同一坐标轴、上下两条轨道），并返回统计数据。
    
    Returns:
        tuple: (fig, stats_dict) 或 (None, None)
    """
    if seg_df is None or seg_df.empty:
        return None, None

    df_chr = seg_df[seg_df['Chr_Clean'].astype(str) == str(chr_name)].copy()
    if df_chr.empty:
        return None, None

    group_name_map = {
        0: 'indica',
        1: 'temperate_japonica',
        2: 'tropical_japonica',
        3: 'basmati',
        4: 'aus',
        5: 'O._rufipogon'
    }

    group_color_map = {
        0: '#1f77b4',
        1: '#ff7f0e',
        2: '#2ca02c',
        3: '#d62728',
        4: '#9467bd',
        5: '#8c564b'
    }

    hap_y = {'hap1': 1, 'hap2': 0}
    fig = go.Figure()

    present_haps = [h for h in ['hap1', 'hap2'] if (df_chr['Haplotype'] == h).any()]
    if not present_haps:
        return None, None

    group_codes = sorted(df_chr['Group'].dropna().unique().tolist())
    shown_legend_groups = set()

    # 计算统计数据
    stats_dict = {}
    
    for hap in present_haps:
        df_h = df_chr[df_chr['Haplotype'] == hap]
        y0 = hap_y.get(hap, 0)
        
        # 计算该单倍体的总长度和各祖源长度
        hap_stats = {}
        total_length = 0
        
        for group_code in group_codes:
            if pd.isna(group_code):
                continue
            group_code_int = int(group_code)
            sub = df_h[df_h['Group'] == group_code]
            if sub.empty:
                continue
            
            # 计算该祖源类型的总长度
            group_length = (sub['End'] - sub['Start']).sum()
            total_length += group_length
            hap_stats[group_name_map.get(group_code_int, str(group_code_int))] = group_length

            show_in_legend = False
            if group_code_int not in shown_legend_groups:
                show_in_legend = True
                shown_legend_groups.add(group_code_int)

            x_vals = []
            y_vals = []
            for _, r in sub.iterrows():
                x_vals.extend([r['Start_Mb'], r['End_Mb'], None])
                y_vals.extend([y0, y0, None])

            fig.add_trace(
                go.Scatter(
                    x=x_vals,
                    y=y_vals,
                    mode='lines',
                    name=group_name_map.get(group_code_int, str(group_code_int)),
                    legendgroup=str(group_code_int),
                    showlegend=show_in_legend,
                    line=dict(color=group_color_map.get(group_code_int, '#7f7f7f'), width=12),
                    hovertemplate=(
                        "Group=%{fullData.name}<br>" +
                        f"Haplotype={hap}<br>" +
                        "Position=%{x:.3f} Mb" +
                        "<extra></extra>"
                    )
                )
            )
        
        # 计算百分比
        hap_stats_with_pct = {}
        for group_name, length in hap_stats.items():
            percentage = (length / total_length * 100) if total_length > 0 else 0
            hap_stats_with_pct[group_name] = {
                'length': length,
                'percentage': percentage
            }
        hap_stats_with_pct['total_length'] = total_length
        stats_dict[hap] = hap_stats_with_pct

    title = f"Chr{chr_name}" if str(chr_name).isdigit() else str(chr_name)
    if lang_code == 'zh':
        title = f"Chr{chr_name}" if str(chr_name).isdigit() else str(chr_name)

    x_max = float(df_chr['End_Mb'].max())
    x_min = float(df_chr['Start_Mb'].min())
    pad = max(1.0, (x_max - x_min) * 0.02)

    fig.update_layout(
        title=title,
        height=260,
        margin=dict(l=10, r=10, t=40, b=30),
        plot_bgcolor='white',
        xaxis=dict(
            title='Position (Mb)' if lang_code != 'zh' else '位置 (Mb)',
            range=[max(0, x_min - pad), x_max + pad],
            showgrid=True,
            gridcolor='rgba(0,0,0,0.06)',
            zeroline=False
        ),
        yaxis=dict(
            title='',
            tickmode='array',
            tickvals=[0, 1],
            ticktext=['hap2', 'hap1'],
            range=[-0.6, 1.6],
            showgrid=False,
            zeroline=False
        ),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0)
    )

    return fig, stats_dict

# --- 侧边栏: 语言选择 ---
st.sidebar.title("🌐 Language")

# 初始化语言状态 (使用代码 'zh'/'en' 以便更灵活控制显示)
if 'lang_code' not in st.session_state:
    # 尝试迁移旧状态
    if 'language' in st.session_state:
        st.session_state.lang_code = 'zh' if st.session_state.language == '中文' else 'en'
    else:
        st.session_state.lang_code = 'zh'

# 根据当前状态决定 Label 和 选项显示
if st.session_state.lang_code == 'zh':
    lang_label = "选择语言"
    # 选项值保持 'zh', 'en'，但显示为 "中文", "英语"
    format_func = lambda x: "中文" if x == 'zh' else "英语"
else:
    lang_label = "Select Language"
    # 选项值保持 'zh', 'en'，但显示为 "Chinese", "English"
    format_func = lambda x: "Chinese" if x == 'zh' else "English"

# 语言选择组件
selected_code = st.sidebar.radio(
    lang_label, 
    ['zh', 'en'],
    index=0 if st.session_state.lang_code == 'zh' else 1,
    format_func=format_func,
    key="lang_radio_code"
)

# 如果选择改变，更新状态并重新运行
if selected_code != st.session_state.lang_code:
    st.session_state.lang_code = selected_code
    st.rerun()

# 保持向下兼容，定义 lang 变量
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
    "Gene Name / ID" if st.session_state.lang_code != 'zh' else "基因名称 / ID", 
    placeholder="e.g. Os01g01010"
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
    stat_col1, stat_col2, stat_col3, stat_col4 = st.columns(4)
    
    total_genes = len(filtered_df)
    favorable = len(filtered_df[filtered_df['Evaluation'] == 'Favorable'])
    unfavorable = len(filtered_df[filtered_df['Evaluation'] == 'Unfavorable'])
    neutral_unknown = total_genes - favorable - unfavorable
    
    with stat_col1:
        st.metric(
            label="📊 " + ("Total Genes" if st.session_state.lang_code != 'zh' else "基因总数"),
            value=f"{total_genes:,}"
        )
    with stat_col2:
        st.metric(
            label="✅ " + ("Favorable" if st.session_state.lang_code != 'zh' else "有利基因"),
            value=f"{favorable:,}",
            delta=f"{favorable/total_genes*100:.1f}%" if total_genes > 0 else "0%"
        )
    with stat_col3:
        st.metric(
            label="❌ " + ("Unfavorable" if st.session_state.lang_code != 'zh' else "不利基因"),
            value=f"{unfavorable:,}",
            delta=f"-{unfavorable/total_genes*100:.1f}%" if total_genes > 0 else "0%"
        )
    with stat_col4:
        st.metric(
            label="⚪ " + ("Neutral/Unknown" if st.session_state.lang_code != 'zh' else "中性/未知"),
            value=f"{neutral_unknown:,}"
        )
    
    # 快速筛选提示
    if total_genes > 0:
        st.caption(
            f"💡 Tip: Use the sidebar filters to focus on specific categories or chromosomes" 
            if st.session_state.lang_code != 'zh' 
            else f"💡 提示：使用侧边栏筛选器可聚焦于特定分类或染色体"
        )
    
    st.markdown("---")
    
    # 布局：左侧/上方是图表，右侧/下方是详情列表
    # 为了更好的交互，我们将图表放在顶部，详情放在下面或侧边
    
    # --- 3. 交互式绘图模块 ---
    
    # 颜色映射
    # 使用用户指定的9种颜色，映射到 gene_info.json 中的中文分类
    trait_color_map = {
        '产量组成相关': '#fb6a4b',
        '植株形态': '#fcbca1',
        '抽穗期': '#fdbb84',
        '生物胁迫': '#9fcbe2',
        '非生物胁迫': '#bbbddc',
        '口感品质': '#9e9bc7',
        '种子形态': '#4292c7',
        '次生代谢相关': '#acde8b',
        '其他': '#015932'
    }
    
    evaluation_color_map = {
        'Favorable': 'red',
        'Unfavorable': 'blue',
        'Neutral': 'black',
        'Unknown': 'gray'
    }
    
    # 准备绘图数据 (仅包含有效坐标的行)
    plot_df = filtered_df[filtered_df['Is_Plot_Valid']].copy()
    
    # 计算染色体长度用于背景 (基于 plot_df)
    if not plot_df.empty:
        # 计算每条染色体的最大位置，并稍微增加一点长度 (例如 + 2Mb)
        chr_lengths = plot_df.groupby('Plot_Chr')['Pos_Mb'].max().reset_index()
        chr_lengths['Max_Pos'] = chr_lengths['Pos_Mb'] + 2
        
        # 排序逻辑：数字染色体按数字排，Other Gene 排最后
        def get_sort_key(val):
            if str(val).isdigit():
                return int(val)
            return 999 # 放在最后

        chr_lengths['Sort_Key'] = chr_lengths['Plot_Chr'].apply(get_sort_key)
        chr_lengths = chr_lengths.sort_values('Sort_Key')
        
        # 映射 Y 轴位置
        chr_lengths['Y_pos'] = chr_lengths['Plot_Chr'].apply(lambda x: chr_y_map.get(str(x)))
        
        # 计算 X 轴最大范围，用于设置图表宽度
        X_MAX = chr_lengths['Max_Pos'].max()
    else:
        X_MAX = 50
        chr_lengths = pd.DataFrame()

    # 为 plot_df 添加一个临时列用于索引追踪，确保点击事件能获取正确的数据行
    plot_df['row_id'] = plot_df.index

    # 翻译 Alt_Allele_Func 用于 hover (使用 map 优化性能)
    plot_df['Alt_Allele_Func_Display'] = plot_df['Alt_Allele_Func'].map(
        lambda x: get_text(x, lang, section='Functions')
    )
    
    # 翻译 Evaluation 用于 hover (始终生成)
    plot_df['Evaluation_Display'] = plot_df['Evaluation'].map(
        lambda x: get_text(x, lang, section='Evaluation')
    )

    # 翻译 Trait_Category 用于 hover (始终生成，避免在 Evaluation 模式下报错)
    plot_df['Trait_Category_Display'] = plot_df['Trait_Category'].map(
        lambda x: get_text(x, lang, section='Categories')
    )

    # 根据模式选择颜色配置
    if view_mode_idx == 0: # 按性状分类
        color_col = 'Trait_Category'
        color_map = trait_color_map
        
        # 更新颜色映射的 Key 为当前语言
        plot_color_map = {
            get_text(k, lang, section='Categories'): v 
            for k, v in trait_color_map.items()
        }
        plot_color_col = 'Trait_Category_Display'
            
    else: # 按评估分类
        color_col = 'Evaluation'
        color_map = evaluation_color_map
        
        # 同样处理评估的翻译 (Evaluation_Display 已经在上面生成了)
        plot_color_map = {
            get_text(k, lang, section='Evaluation'): v 
            for k, v in evaluation_color_map.items()
        }
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

        # 更新点的大小，突出显示重点基因
        # 首先设置默认样式 (移除固定 size)
        fig.update_traces(marker=dict(opacity=0.8, line=dict(width=1, color='white')))
        
        # 遍历所有 trace，根据 GeneName 更新 size 和 color
        # 注意：Plotly Express 生成的 trace 可能已经被分组，所以我们必须小心处理颜色
        # 如果要高亮某个点，最简单的方法是添加一个新的 Trace
        
        # 1. 设置基础大小（根据 QTN 效应值）
        for trace in fig.data:
            if trace.customdata is not None:
                sizes = []
                # customdata 是一个 numpy 数组或列表的列表
                for pt_data in trace.customdata:
                    row_id = pt_data[0]
                    if row_id in plot_df.index:
                        # 构造 key: "Chr" + Chr_Clean + "_" + Pos_7.0
                        chr_val = plot_df.loc[row_id, 'Chr']
                        pos_val = str(plot_df.loc[row_id, 'Pos_7.0']).strip()
                        key = f"{chr_val}_{pos_val}"
                        
                        effect = qtn_effects.get(key, None)
                        
                        if effect == 'high':
                            sizes.append(18)  # 主效位点（最大）
                        elif effect == 'low':
                            sizes.append(13)  # 次要位点（中等）
                        else:
                            sizes.append(8)   # 普通位点（最小）
                    else:
                        sizes.append(6)
                trace.marker.size = sizes

        # 2. 如果有搜索结果，添加高亮 Trace
        if search_gene_query:
            # 模糊匹配：只要包含搜索词就算
            # 优先精确匹配
            highlight_mask = plot_df['GeneName'].astype(str).str.contains(search_gene_query, case=False, na=False)
            highlight_df = plot_df[highlight_mask]
            
            if not highlight_df.empty:
                # 添加一个红色的高亮层
                fig.add_trace(
                    go.Scatter(
                        x=highlight_df['Pos_Mb'],
                        y=highlight_df['Y_pos'],
                        mode='markers',
                        marker=dict(
                            size=20,
                            color='red',
                            symbol='circle-open', # 空心圆圈
                            line=dict(width=3, color='red')
                        ),
                        name='Search Result',
                        hoverinfo='skip', # 悬停不显示额外信息，利用底层的点显示
                        showlegend=True
                    )
                )
                st.toast(
                    f"Found {len(highlight_df)} genes matching '{search_gene_query}'" 
                    if st.session_state.lang_code != 'zh' 
                    else f"找到 {len(highlight_df)} 个匹配 '{search_gene_query}' 的基因"
                )
            else:
                if search_gene_query:
                     st.toast(
                        f"No genes found matching '{search_gene_query}'"
                        if st.session_state.lang_code != 'zh' 
                        else f"未找到匹配 '{search_gene_query}' 的基因"
                     )
        
        # 确保 Y 轴标签正确
        # 根据 chr_y_map 反向生成标签列表
        # 仅显示当前 plot_df 中存在的染色体
        present_y_pos = set(plot_df['Y_pos'].dropna().unique())
        # 如果有背景线，也要包含背景线的 Y_pos
        if not chr_lengths.empty:
             present_y_pos.update(chr_lengths['Y_pos'].dropna().unique())
             
        y_ticks = sorted([item for item in chr_y_map.items() if item[1] in present_y_pos], key=lambda x: x[1])
        
        tick_vals = [item[1] for item in y_ticks]
        tick_text = []
        for name, _ in y_ticks:
            if name == 'Other Gene':
                tick_text.append('Other Gene')
            elif name.isdigit():
                tick_text.append(f'Chr{name}')
            else:
                tick_text.append(str(name))
        
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
    
    # 如果用户在图表上选择了点 (Streamlit 1.35+ 支持 on_select)
    try:
        # 尝试解析 selection
        if selection and 'selection' in selection and 'points' in selection['selection']:
             points = selection['selection']['points']
             # 提取 custom_data 中的 row_id
             selected_row_ids = [p['customdata'][0] for p in points if 'customdata' in p]
             
             if selected_row_ids:
                 # 使用 row_id 筛选数据 (注意：这里要从 filtered_df 中筛选，因为 plot_df 是 filtered_df 的子集)
                 # 但 row_id 是 plot_df 的 index，而 plot_df 是 filtered_df 的切片，所以 index 是一致的
                 selected_data = filtered_df.loc[selected_row_ids]
                 st.info(get_text("selected_info", lang, format_args=[len(selected_data)]))
                 # 仅显示选中的数据
                 display_df = selected_data[display_cols].copy()
             else:
                 display_df = filtered_df[display_cols].copy()
        else:
            display_df = filtered_df[display_cols].copy()
    except Exception as e:
        # 降级处理或错误提示
        # st.error(f"解析选择时出错: {e}")
        display_df = filtered_df[display_cols].copy()

    # 翻译表格内容 (使用 map 优化性能)
    # 1. Trait_Category
    display_df['Trait_Category'] = display_df['Trait_Category'].map(lambda x: get_text(x, lang, section='Categories'))
    # 2. Evaluation
    display_df['Evaluation'] = display_df['Evaluation'].map(lambda x: get_text(x, lang, section='Evaluation'))
    # 3. Alt_Allele_Func
    display_df['Alt_Allele_Func'] = display_df['Alt_Allele_Func'].map(lambda x: get_text(x, lang, section='Functions'))
    
    # 如果有搜索查询，在表格中也筛选出匹配项
    if search_gene_query:
        search_mask = display_df['GeneName'].astype(str).str.contains(search_gene_query, case=False, na=False)
        if search_mask.any():
            display_df = display_df[search_mask].copy()
            st.info(
                f"🔍 Table filtered to show {len(display_df)} genes matching '{search_gene_query}'"
                if st.session_state.lang_code != 'zh'
                else f"🔍 表格已筛选，显示 {len(display_df)} 个匹配 '{search_gene_query}' 的基因"
            )
    
    # 处理 RAP 和 MSU 链接
    def transform_rap_url(val):
        s = str(val).strip()
        if s and s.upper() not in ['NA', 'N/A', '']:
            return f"https://rapdb.dna.affrc.go.jp/locus/?name={s}"
        return None

    def transform_msu_url(val):
        s = str(val).strip()
        if s and s.upper() not in ['NA', 'N/A', '']:
            return f"https://rice.uga.edu/cgi-bin/ORF_infopage.cgi?orf={s}"
        return None

    display_df['RAP_Locus'] = display_df['RAP_Locus'].apply(transform_rap_url)
    display_df['MSU_Locus'] = display_df['MSU_Locus'].apply(transform_msu_url)

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
            - Export-ready publication charts
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
            - 可导出的出版级图表
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
                    
                    # 展示统计数据
                    if stats:
                        st.markdown("**" + ("Ancestry Composition" if st.session_state.lang_code != 'zh' else "祖源组成比例") + "**")
                        
                        # 为 hap1 和 hap2 创建两列
                        col_hap1, col_hap2 = st.columns(2)
                        
                        with col_hap1:
                            st.markdown("**Hap1**")
                            if 'hap1' in stats:
                                hap1_data = stats['hap1']
                                total_len = hap1_data.get('total_length', 0)
                                
                                for group_name in sorted(hap1_data.keys()):
                                    if group_name == 'total_length':
                                        continue
                                    group_info = hap1_data[group_name]
                                    st.metric(
                                        label=group_name.replace('_', ' ').title(),
                                        value=f"{group_info['percentage']:.1f}%",
                                        delta=f"{group_info['length']/1_000_000:.2f} Mb"
                                    )
                        
                        with col_hap2:
                            st.markdown("**Hap2**")
                            if 'hap2' in stats:
                                hap2_data = stats['hap2']
                                total_len = hap2_data.get('total_length', 0)
                                
                                for group_name in sorted(hap2_data.keys()):
                                    if group_name == 'total_length':
                                        continue
                                    group_info = hap2_data[group_name]
                                    st.metric(
                                        label=group_name.replace('_', ' ').title(),
                                        value=f"{group_info['percentage']:.1f}%",
                                        delta=f"{group_info['length']/1_000_000:.2f} Mb"
                                    )
                else:
                    st.info(f"No segment data for Chr {chr_name}" if st.session_state.lang_code != 'zh' else f"Chr {chr_name} 无祖源分段数据")
    else:
        # 单条染色体直接展示
        chr_name = chrs_to_show[0]
        fig, stats = plot_segments_for_chr_both_haps(segments_df, chr_name, st.session_state.lang_code)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
            
            # 展示统计数据
            if stats:
                st.markdown("**" + ("Ancestry Composition" if st.session_state.lang_code != 'zh' else "祖源组成比例") + "**")
                
                # 为 hap1 和 hap2 创建两列
                col_hap1, col_hap2 = st.columns(2)
                
                with col_hap1:
                    st.markdown("**Hap1**")
                    if 'hap1' in stats:
                        hap1_data = stats['hap1']
                        total_len = hap1_data.get('total_length', 0)
                        
                        for group_name in sorted(hap1_data.keys()):
                            if group_name == 'total_length':
                                continue
                            group_info = hap1_data[group_name]
                            st.metric(
                                label=group_name.replace('_', ' ').title(),
                                value=f"{group_info['percentage']:.1f}%",
                                delta=f"{group_info['length']/1_000_000:.2f} Mb"
                            )
                
                with col_hap2:
                    st.markdown("**Hap2**")
                    if 'hap2' in stats:
                        hap2_data = stats['hap2']
                        total_len = hap2_data.get('total_length', 0)
                        
                        for group_name in sorted(hap2_data.keys()):
                            if group_name == 'total_length':
                                continue
                            group_info = hap2_data[group_name]
                            st.metric(
                                label=group_name.replace('_', ' ').title(),
                                value=f"{group_info['percentage']:.1f}%",
                                delta=f"{group_info['length']/1_000_000:.2f} Mb"
                            )
        else:
            st.info(f"No segment data for Chr {chr_name}" if st.session_state.lang_code != 'zh' else f"Chr {chr_name} 无祖源分段数据")

