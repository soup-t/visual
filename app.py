import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
import json

# 设置页面配置
st.set_page_config(page_title="Rice GWAS Interactive Visualization", layout="wide")

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
        # 为了方便后续处理，我们标记是否可绘图
        
        def is_valid_chr(c):
            return str(c).isdigit() and 1 <= int(c) <= 12

        df['Is_Plot_Valid'] = df['Chr_Clean'].apply(is_valid_chr) & df['Pos_7.0_Num'].notna()
        
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

# --- 侧边栏: 数据上传与分类选择 ---

st.sidebar.title("🌾 水稻基因组可视化")

# 1. 文件上传部分
st.sidebar.header("1. 数据来源")
data_source = st.sidebar.radio("选择数据:", ["使用 Demo 数据", "上传 .geno 文件"])

df = None
chr_y_map = None
sample_col = None

if data_source == "使用 Demo 数据":
    demo_path = os.path.join(os.path.dirname(__file__), 'demo data', 'HHZ (3).geno')
    if os.path.exists(demo_path):
        with open(demo_path, 'r') as f:
            df, chr_y_map, sample_col = load_data(f)
        st.sidebar.success("已加载 Demo 数据: HHZ (3).geno")
    else:
        st.sidebar.error("Demo 文件未找到，请检查路径。")
else:
    uploaded_file = st.sidebar.file_uploader("上传你的 .geno 文件", type=['geno', 'txt', 'csv', 'tsv'])
    if uploaded_file is not None:
        df, chr_y_map, sample_col = load_data(uploaded_file)
        st.sidebar.success(f"已加载文件: {uploaded_file.name}")

# --- 主界面逻辑 ---

if df is not None:
    # 2. 可视化模式选择
    st.sidebar.header("2. 可视化模式")
    view_mode = st.sidebar.radio("选择着色模式:", ["按性状分类 (Trait Category)", "按有利/不利评估 (Evaluation)"])

    # 3. 分类筛选 (模拟左侧大分类)
    st.sidebar.header("3. 基因功能分类筛选")
    
    # 获取所有可用分类
    all_categories = ['All'] + sorted(df['Trait_Category'].unique().tolist())
    
    # 使用 radio 模拟左侧菜单栏的效果
    selected_category = st.sidebar.radio("选择查看的性状类别:", all_categories)

    # 根据选择筛选数据
    if selected_category != 'All':
        filtered_df = df[df['Trait_Category'] == selected_category].copy()
        main_title = f"{selected_category} 相关基因位点"
    else:
        filtered_df = df.copy()
        main_title = "全基因组重要位点概览"

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

    # 根据模式选择颜色配置
    if view_mode == "按性状分类 (Trait Category)":
        color_col = 'Trait_Category'
        color_map = trait_color_map
    else:
        color_col = 'Evaluation'
        color_map = evaluation_color_map

    # 绘制散点图
    if not plot_df.empty:
        fig = px.scatter(
            plot_df, 
            x='Pos_Mb', 
            y='Y_pos', 
            color=color_col,
            color_discrete_map=color_map,
            hover_name="GeneName", 
            custom_data=['row_id'], # 将 row_id 传入 custom_data
            hover_data={
                'Chr_Label': True,          
                'Pos_Mb': ':.2f',           
                'GeneName': False,          
                'Alt_Allele_Func': True,
                sample_col: True, # 显示样本基因型
                'Evaluation': True, # 显示评估结果
                'Y_pos': False,             
                'Chr': False, # 原始 Chr 列可能包含非数字，这里已被过滤，但 hover data 引用的是列名
                'Chr_Clean': False,
                'Link': False, # 不在 hover 中直接显示长链接
                'row_id': False # 隐藏 row_id
            },
            labels={
                'Pos_Mb': 'Position (Mb)', 
                'Trait_Category': 'Trait Category',
                'Alt_Allele_Func': 'Function',
                sample_col: 'genetype',
                'Evaluation': 'Evaluation'
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
            xaxis=dict(title="Position (Mb)", range=[0, X_MAX]),
            plot_bgcolor='white',
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            clickmode='event+select' # 允许点击选择
        )

        # 显示图表
        selection = st.plotly_chart(fig, use_container_width=True, on_select="rerun")
    else:
        st.warning("当前筛选条件下没有可绘图的数据 (缺少有效的染色体或位置信息)。")
        selection = None
    
    # --- 4. 详情与链接模块 ---
    
    st.markdown("---")
    st.subheader(f"🧬 {selected_category if selected_category != 'All' else '所有'} 基因详情列表")
    st.caption("点击上方图表中的点，或在下方列表中查找。")

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
                 st.info(f"你选择了 {len(selected_data)} 个位点")
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

    # 将 Link 列转换为可点击的 HTML 链接
    # Streamlit 的 dataframe 组件支持 LinkColumn
    
    # 重命名 Link 列为 "More Info" 或其他，并在 column_config 中设置 display_text
    display_df = display_df.rename(columns={'Link': 'More Info'})
    
    st.dataframe(
        display_df,
        column_config={
            "More Info": st.column_config.LinkColumn(
                "More Info",
                help="点击查看详情",
                display_text="View details"
            ),
            "Pos_7.0": st.column_config.TextColumn(
                "Position"
            )
        },
        hide_index=True,
        use_container_width=True
    )

else:
    st.info("👈 请在左侧侧边栏上传文件或选择 Demo 数据以开始。")

