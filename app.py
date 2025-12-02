import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import json
import random

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

def categorize_trait(row):
    """
    根据 GeneName 查找数据库分类，如果找不到则归类为 '其他'
    """
    gene_name = str(row['GeneName']).strip()
    
    # 优先使用数据库中的分类
    if gene_name in gene_db:
        db_trait = gene_db[gene_name].get('Trait')
        if db_trait:
            return db_trait
            
    return '其他'

def evaluate_gene(row):
    """
    评估基因型是有利还是不利
    """
    gene_name = str(row['GeneName']).strip()
    sample_geno = str(row['genetype']).strip()
    
    if gene_name in gene_db:
        info = gene_db[gene_name]
        ref_status = info.get('RefStatus')
        
        if not ref_status: return "Unknown"
        
        if ref_status == "Neutral":
            return "Neutral"
            
        if sample_geno == '0|0':
            return ref_status
        else:
            # 反转状态
            if ref_status == "Favorable": return "Unfavorable"
            if ref_status == "Unfavorable": return "Favorable"
            
    return "Unknown"

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

        # 数据清洗
        # 注意：不再直接删除 Pos_7.0 为空的行，以便在列表中展示
        # 但为了绘图，我们需要处理 Pos_7.0
        
        # 尝试清理 Chr
        df['Chr_Clean'] = df['Chr'].astype(str).str.replace('Chr', '', case=False).str.strip()
        
        # 处理 Pos_7.0
        df['Pos_7.0_Clean'] = df['Pos_7.0'].apply(lambda x: str(x).split('-')[0])
        df['Pos_7.0_Num'] = pd.to_numeric(df['Pos_7.0_Clean'], errors='coerce').astype('Int64')
        
        # 计算 Pos_Mb (仅针对有效数值)
        df['Pos_Mb'] = df['Pos_7.0_Num'] / 1_000_000

        # 仅保留主要的染色体 (1到12) 用于绘图逻辑的 Y 轴映射
        # 但原始数据保留在 df 中
        
        def is_valid_chr(c):
            return str(c).isdigit() and 1 <= int(c) <= 12

        # 标记有效的染色体
        df['Is_Valid_Chr'] = df['Chr_Clean'].apply(is_valid_chr)

        # 处理未定位但有染色体信息的基因 (Pos_7.0_Num 为空 但 Is_Valid_Chr 为真)
        # 将它们随机分布在 46-49Mb 之间 (用户要求 45-50Mb)
        def assign_fallback_pos(row):
            if row['Is_Valid_Chr'] and pd.isna(row['Pos_Mb']):
                return random.uniform(46, 49)
            return row['Pos_Mb']
            
        df['Pos_Mb'] = df.apply(assign_fallback_pos, axis=1)

        df['Is_Plot_Valid'] = df['Is_Valid_Chr'] & df['Pos_Mb'].notna()
        
        # 仅对可绘图的数据进行排序优化，但我们需要保持原始 df 完整
        # 这里我们只对 Chr 进行标准化以便后续处理，不删除行
        
        # 填充缺失值
        df['Alt_Allele_Func'] = df['Alt_Allele_Func'].fillna('Unknown Function')
        df['GeneName'] = df['GeneName'].fillna('N/A')
        
        # 应用分类 (现在传入整行)
        df['Trait_Category'] = df.apply(categorize_trait, axis=1)
        
        # 应用评估 (Favorable/Unfavorable)
        df['Evaluation'] = df.apply(evaluate_gene, axis=1)
        
        # 生成 Y 轴位置 (仅针对有效染色体)
        valid_chromosomes = sorted([int(c) for c in df[df['Is_Plot_Valid']]['Chr_Clean'].unique()])
        chr_y_map = {chr_num: i + 1 for i, chr_num in enumerate(valid_chromosomes)}
        
        def get_y_pos(row):
            if row['Is_Plot_Valid']:
                return chr_y_map.get(int(row['Chr_Clean']))
            return None

        df['Y_pos'] = df.apply(get_y_pos, axis=1)
        df['Chr_Label'] = df['Chr_Clean'].apply(lambda x: f"Chr{x}" if str(x).isdigit() else str(x))
        
        # 生成占位符链接 (后续可替换为真实链接逻辑)
        # 假设链接格式为: http://www.ricedata.cn/gene/{GeneName}
        df['Link'] = df['GeneName'].apply(lambda x: f"http://www.ricedata.cn/gene/search?key={x}" if x != 'N/A' else "#")

        return df, chr_y_map, sample_col

    except Exception as e:
        st.error(f"读取文件时出错: {e}")
        return None, None, None

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

if data_source_idx == 0: # 使用 Demo 数据
    demo_path = os.path.join(os.path.dirname(__file__), 'demo data', 'HHZ.geno')
    if os.path.exists(demo_path):
        with open(demo_path, 'r') as f:
            df, chr_y_map, sample_col = load_data(f)
        st.sidebar.success(get_text("load_success_demo", lang))
    else:
        st.sidebar.error(get_text("file_not_found", lang))
else: # 上传文件
    uploaded_file = st.sidebar.file_uploader(get_text("upload_label", lang), type=['geno', 'txt', 'csv', 'tsv'])
    if uploaded_file is not None:
        df, chr_y_map, sample_col = load_data(uploaded_file)
        st.sidebar.success(get_text("load_success_file", lang, format_args=[uploaded_file.name]))

# --- 主界面逻辑 ---

if df is not None:
    # 2. 可视化模式选择
    st.sidebar.header(get_text("view_mode_header", lang))
    view_mode_options = [get_text("mode_trait", lang), get_text("mode_eval", lang)]
    view_mode_idx = st.sidebar.radio(
        get_text("view_mode_label", lang), 
        range(len(view_mode_options)), 
        format_func=lambda x: view_mode_options[x]
    )

    # 3. 分类筛选 (模拟左侧大分类)
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

    # --- 主区域 ---
    
    st.title(main_title)
    
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
        chr_lengths = plot_df.groupby('Chr_Clean')['Pos_Mb'].max().reset_index()
        # 确保 Chr_Clean 是数字以便排序
        chr_lengths['Chr_Num'] = pd.to_numeric(chr_lengths['Chr_Clean'])
        chr_lengths = chr_lengths.sort_values('Chr_Num')
        chr_lengths['Y_pos'] = chr_lengths['Chr_Num'].map(chr_y_map)
        # 将最大长度限制为 50MB
        X_MAX = 50
    else:
        X_MAX = 50
        chr_lengths = pd.DataFrame()

    # 为 plot_df 添加一个临时列用于索引追踪，确保点击事件能获取正确的数据行
    plot_df['row_id'] = plot_df.index

    # 翻译 Alt_Allele_Func 用于 hover
    plot_df['Alt_Allele_Func_Display'] = plot_df['Alt_Allele_Func'].apply(
        lambda x: get_text(x, lang, section='Functions')
    )
    
    # 翻译 Evaluation 用于 hover (始终生成)
    plot_df['Evaluation_Display'] = plot_df['Evaluation'].apply(
        lambda x: get_text(x, lang, section='Evaluation')
    )

    # 翻译 Trait_Category 用于 hover (始终生成，避免在 Evaluation 模式下报错)
    plot_df['Trait_Category_Display'] = plot_df['Trait_Category'].apply(
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
                'Link': False, 
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
                            x1=X_MAX, y1=row['Y_pos'],
                            line=dict(color="lightgray", width=10),
                            layer="below"
                        )
                    )

        fig.update_traces(marker=dict(size=10, opacity=0.8, line=dict(width=1, color='white')))
        
        # 确保 Y 轴标签正确
        tick_vals = list(chr_y_map.values())
        tick_text = [f'Chr{c}' for c in sorted(chr_y_map.keys())]
        
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
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            clickmode='event+select' # 允许点击选择
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
    display_cols = ['Chr', 'Pos_7.0', 'GeneName', 'Alt_Allele_Func', 'Trait_Category', sample_col, 'Evaluation']
    
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
                 # 同时也需要 Link 列来生成 View details
                 display_df['Link'] = selected_data['Link']
             else:
                 display_df = filtered_df[display_cols].copy()
                 display_df['Link'] = filtered_df['Link']
        else:
            display_df = filtered_df[display_cols].copy()
            display_df['Link'] = filtered_df['Link']
    except Exception as e:
        # 降级处理或错误提示
        # st.error(f"解析选择时出错: {e}")
        display_df = filtered_df[display_cols].copy()
        display_df['Link'] = filtered_df['Link']

    # 翻译表格内容
    # 1. Trait_Category
    display_df['Trait_Category'] = display_df['Trait_Category'].apply(lambda x: get_text(x, lang, section='Categories'))
    # 2. Evaluation
    display_df['Evaluation'] = display_df['Evaluation'].apply(lambda x: get_text(x, lang, section='Evaluation'))
    # 3. Alt_Allele_Func
    display_df['Alt_Allele_Func'] = display_df['Alt_Allele_Func'].apply(lambda x: get_text(x, lang, section='Functions'))
    
    # 重命名列
    col_map = {
        'Chr': get_text("Chr", lang, section='Columns'),
        'Pos_7.0': get_text("Pos_7.0", lang, section='Columns'),
        'GeneName': get_text("GeneName", lang, section='Columns'),
        'Alt_Allele_Func': get_text("Alt_Allele_Func", lang, section='Columns'),
        'Trait_Category': get_text("Trait_Category", lang, section='Columns'),
        sample_col: get_text("genetype", lang, section='Columns'),
        'Evaluation': get_text("Evaluation", lang, section='Columns'),
        'Link': get_text("More Info", lang, section='Columns')
    }
    display_df = display_df.rename(columns=col_map)
    
    # 获取翻译后的列名用于配置
    more_info_col = get_text("More Info", lang, section='Columns')
    pos_col = get_text("Pos_7.0", lang, section='Columns')

    st.dataframe(
        display_df,
        column_config={
            more_info_col: st.column_config.LinkColumn(
                more_info_col,
                help=get_text("col_more_info_help", lang),
                display_text="View details"
            ),
            pos_col: st.column_config.TextColumn(
                get_text("col_position", lang)
            )
        },
        hide_index=True,
        use_container_width=True
    )

else:
    st.info(get_text("start_prompt", lang))

