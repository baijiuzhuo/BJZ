# GFIP 详细技术路线图 (Gene Family Identification Pipeline)

> **版本**: v3.1 (Quad-Core Edition)  
> **最后更新**: 2026-01-24  
> **文档目的**: 详细记录 Pipeline 每一步的技术细节、数据流、处理逻辑和关键参数

---

## 目录

1. [总体架构](#总体架构)
2. [Phase 1: 种子序列获取](#phase-1-种子序列获取)
3. [Phase 2: HMM模型构建](#phase-2-hmm模型构建)
4. [Phase 3: 四核心搜索](#phase-3-四核心搜索)
5. [Phase 4: 域验证](#phase-4-域验证)
6. [Phase 5: 序列提取](#phase-5-序列提取)
7. [Phase 6: 多序列比对与Motif分析](#phase-6-多序列比对与motif分析)
8. [Phase 7: 系统发育树构建](#phase-7-系统发育树构建)
9. [Phase 8: Ka/Ks选择压力分析](#phase-8-kaks选择压力分析)
10. [Phase 9: 共线性分析](#phase-9-共线性分析)
11. [Phase 10: 启动子分析](#phase-10-启动子分析)
12. [关键技术细节备忘](#关键技术细节备忘)

---

## 总体架构

```
用户输入 (query, domains, genome, proteome, cds, gff)
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 1: 种子获取 (retrieve_seeds.py)                                │
│   三轨并行: NCBI + UniProt + InterPro                                │
│   分类输出: seeds_gold.fasta + seeds_broad.fasta                     │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 2: HMM构建 (build_hmm.py)                                      │
│   去重 → 长度过滤 → MAFFT比对 → Gap修剪 → hmmbuild                    │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 3: 四核心搜索 (search_extract.py x 2 + BLAST x 2)              │
│   HMM-Gold + HMM-Broad + BLAST-Gold + BLAST-Broad                    │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 4: 域验证 (scan_cdd_ncbi.py + interproscan_runner.py)          │
│   CDD + InterPro → Union/Intersection                                │
└─────────────────────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────────────────────┐
│ Phase 5: 序列提取 (universal_family_extractor.py)                    │
│   最长异构体过滤 + ID映射 → PEP/CDS/Gene/Promoter                     │
└─────────────────────────────────────────────────────────────────────┘
    │
    ├────────────────┬─────────────────┬─────────────────┐
    ▼                ▼                 ▼                 ▼
┌─────────┐   ┌───────────┐    ┌───────────┐    ┌────────────┐
│ MSA     │   │ Ka/Ks     │    │ 共线性     │    │ 启动子     │
│ 系统发育 │   │ 选择压力   │    │ 分析      │    │ 分析       │
└─────────┘   └───────────┘    └───────────┘    └────────────┘
    │
    ▼
HTML报告 + 可视化
```

---

## Phase 1: 种子序列获取

### 脚本: `retrieve_seeds.py` (463行)

#### 1.1 数据源与API

| 数据源 | 函数 | API端点 | 认证 |
|--------|------|---------|------|
| **NCBI Protein** | `get_ncbi_many()` | Entrez ESearch + EFetch | API Key (可选) |
| **UniProt** | `get_uniprot_seeds_by_name()` | `rest.uniprot.org/uniprotkb/search` | 无 |
| **InterPro** | `get_interpro_seeds()` | `ebi.ac.uk/interpro/api` | 无 |

#### 1.2 请求重试机制

```python
# get_session() - 创建带重试逻辑的Session
session = requests.Session()
retries = Retry(
    total=5,                # 最多重试5次
    backoff_factor=1,       # 1s, 2s, 4s, 8s, 16s 指数退避
    status_forcelist=[429, 500, 502, 503, 504],  # 这些状态码触发重试
    allowed_methods=["HEAD", "GET", "POST"]
)
```

#### 1.3 NCBI速率限制处理

```python
# 关键参数
batch_size = 400                    # 每批下载的序列数
max_workers = 3                     # 并发线程数 (减少以避免429)
wait_time = 2 + random.uniform(0, 1)  # 查询间隔 (2-3秒)

# 429错误处理
if resp.status_code == 429:
    time.sleep(10 + retry_count * 5 + random.uniform(0, 2))  # 10-30秒+
```

#### 1.4 种子分类逻辑 (Gold vs Silver)

```python
def classify_sequence(rec, source):
    """
    GOLD (高可信度):
    ├── NCBI: NP_* 或 YP_* (已验证的RefSeq)
    ├── UniProt: Swiss-Prot (reviewed:true)
    └── InterPro: /protein/reviewed/ 端点

    SILVER (预测/未审核):
    ├── NCBI: XP_* (预测蛋白), WP_* (非冗余)
    ├── UniProt: TrEMBL (reviewed:false)
    └── InterPro: /protein/unreviewed/ 端点
    """
```

#### 1.5 InterPro分批下载

```python
# fetch_interpro_endpoint() 详细流程:
# 1. 收集Accession IDs (分页请求)
while next_url and len(accessions) < limit:
    resp = session.get(next_url, headers={"Accept": "application/json"})
    for res in data["results"]:
        accessions.append(res["metadata"]["accession"])
    next_url = data.get("next")

# 2. 批量下载序列
batch_size = 100
max_workers = 20  # 高并发下载
for batch in batches:
    url = "https://rest.uniprot.org/uniprotkb/accessions"
    params = {"accessions": ",".join(batch), "format": "fasta"}
```

#### 1.6 去重与输出

```python
# 基于序列内容去重 (非ID)
seen = set()
for rec in all_records:
    seq_str = str(rec.seq).upper()
    if seq_str in seen: continue
    seen.add(seq_str)
    
    # 根据分类标签分配
    if rec._classification == "GOLD":
        gold_recs.append(rec)
    else:
        silver_recs.append(rec)

# 输出文件:
# - {family}_seeds_gold.fasta: 仅Gold序列
# - {family}_seeds_broad.fasta: Gold + Silver (受max_seeds限制)
```

---

## Phase 2: HMM模型构建

### 脚本: `build_hmm.py` (212行)

#### 2.1 序列合并与去重

```python
def merge_and_deduplicate(input_files, output_file):
    seen_sequences = set()
    unique_records = []
    
    for record in SeqIO.parse(file_path, "fasta"):
        seq_str = str(record.seq).strip().upper()
        if seq_str not in seen_sequences:
            seen_sequences.add(seq_str)
            unique_records.append(record)
```

#### 2.2 ⚠️ 长度过滤 (关键步骤)

```python
# 基于中位数长度过滤异常值
lengths = [len(r.seq) for r in unique_records]
lengths.sort()
median_len = lengths[len(lengths)//2]

# 保留范围: 中位数的 75% ~ 125%
min_len = int(median_len * 0.75)
max_len = int(median_len * 1.25)

filtered = [r for r in unique_records if min_len <= len(r.seq) <= max_len]

# 示例: 中位数=500aa → 保留375-625aa的序列
```

#### 2.3 MAFFT多序列比对

```bash
# 实际执行命令
mafft --auto --thread {cpu} seeds_merged.fasta > seeds_aligned.fasta

# --auto: 自动选择最佳算法
#   - <200条: L-INS-i (最准确)
#   - 200-2000条: FFT-NS-i
#   - >2000条: FFT-NS-2 (最快)
```

#### 2.4 ⚠️ MSA Gap修剪 (HMM构建前)

```python
def clean_msa(input_file, output_file, max_gap_fraction=0.9):
    """
    移除Gap比例>90%的列
    (注意: 这里用90%是为HMM构建优化, 系统发育树用50%)
    """
    for i in range(seq_len):
        col = [rec.seq[i] for rec in alignment]
        gap_count = col.count("-")
        if gap_count / len(alignment) <= max_gap_fraction:
            keep_cols.append(i)
```

#### 2.5 hmmbuild构建HMM

```bash
# 实际执行命令
hmmbuild -n {family_name} family.hmm seeds_aligned.fasta

# 输出: family_gold.hmm, family_broad.hmm
```

---

## Phase 3: 四核心搜索

### 脚本: `search_extract.py` (140行)

#### 3.1 四条搜索流配置

| 流 | 模型/种子 | 工具 | 输出文件 |
|----|-----------|------|----------|
| **Stream 1** | `family_gold.hmm` | hmmsearch | `hits_hmm_gold.fasta` |
| **Stream 2** | `family_broad.hmm` | hmmsearch | `hits_hmm_broad.fasta` |
| **Stream 3** | `seeds_gold.fasta` | BLASTp | `hits_blast_gold.fasta` |
| **Stream 4** | `seeds_broad.fasta` | BLASTp | `hits_blast_broad.fasta` |

#### 3.2 hmmsearch执行细节

```python
def run_hmmsearch(hmm_file, target_proteome, output_tbl, threads=4, evalue=1e-5):
    cmd = [
        "hmmsearch",
        "--tblout", output_tbl,     # 表格输出
        "--cpu", str(threads),
        "-E", str(evalue),          # E-value阈值
        hmm_file,
        target_proteome
    ]
    # stdout重定向到DEVNULL以避免大量比对输出
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
```

#### 3.3 结果解析与序列提取

```python
def parse_and_extract(tbl_file, target_proteome, output_fasta, evalue_cutoff):
    valid_ids = set()
    
    with open(tbl_file, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            
            target_name = parts[0]  # 列0: 目标序列名
            evalue = float(parts[4])  # 列4: Full sequence E-value
            
            if evalue <= evalue_cutoff:
                valid_ids.add(target_name)
    
    # 使用索引提取序列 (内存高效)
    proteome_dict = SeqIO.index(target_proteome, "fasta")
    for pid in valid_ids:
        if pid in proteome_dict:
            extracted_records.append(proteome_dict[pid])
```

#### 3.4 候选合并

```python
# run_pipeline_v3.py 中的合并逻辑
hits_files = [hmm_gold, hmm_broad, blast_gold, blast_broad]
seen_ids = set()

for hit_file in hits_files:
    if os.path.exists(hit_file):
        for rec in SeqIO.parse(hit_file, "fasta"):
            if rec.id not in seen_ids:
                seen_ids.add(rec.id)
                all_candidates.append(rec)
```

---

## Phase 4: 域验证

### 4.1 CDD扫描

#### 脚本: `scan_cdd_ncbi.py` (191行)

##### 4.1.1 提交流程

```python
def submit_search(fasta_file, db="cdd", evalue=0.01):
    # 序列清洗: 移除终止密码子符号
    clean_seq = str(seq.seq).replace("*", "")
    
    payload = {
        "db": db,           # cdd, pfam, smart 等
        "smode": "auto",    # 自动模式选择
        "useid1": "true",   # 使用用户提供的ID
        "filter": "true",   # 低复杂度过滤
        "evalue": str(evalue),
        "tdata": "hits",    # 命中表格输出
        "queries": query_str
    }
    
    resp = requests.post(CDSEARCH_URL, data=payload, timeout=30)
    
    # 从响应中提取Job ID (多种格式)
    # 格式1: Search-ID: QM3-xxx
    # 格式2: value="QM3-xxx" name="cdsid"
    # 格式3: #cdsid QM3-xxx
```

##### 4.1.2 ⚠️ 指数退避重试

```python
# 提交失败时的重试策略
max_retries = 5
for attempt in range(max_retries):
    cdsid = submit_search(args.input, args.db, args.evalue)
    
    if cdsid: break
    
    # 指数退避: 30s, 60s, 120s, 240s, 480s + 抖动
    wait_time = 30 * (2 ** attempt) + random.uniform(-5, 5)
    wait_time = max(10, wait_time)  # 最少10秒
    time.sleep(wait_time)
```

##### 4.1.3 结果轮询

```python
def retrieve_results(cdsid, output_file):
    while True:
        resp = requests.post(CDSEARCH_URL, data={"cdsid": cdsid, "tdata": "hits"})
        
        # 检查是否完成 (包含表头)
        if "Q#" in content and "Hit type" in content:
            # 结果就绪
            with open(output_file, "w") as f:
                f.write(content)
            return True
            
        time.sleep(10)  # 每10秒轮询一次
```

### 4.2 InterPro扫描

#### 脚本: `interproscan_runner.py` (345行)

##### 4.2.1 两种运行模式

| 模式 | 参数 | 说明 |
|------|------|------|
| **API模式** | 默认 | 使用EBI REST API |
| **本地模式** | `--local_path` | 使用本地interproscan.sh |

##### 4.2.2 ⚠️ 序列清洗 (关键步骤)

```python
# InterProScan对输入序列非常敏感
# 问题字符会导致提交失败

# 全局清洗 (run_robust_pipeline)
for s in sequences:
    # 1. 移除终止密码子 (*)
    # 2. 移除Gap字符 (.)
    # 3. 去除空白
    cleaned_seq = str(s.seq).upper().replace('*', '').replace('.', '').strip()
    
    # 4. 过滤短序列
    if len(cleaned_seq) > 10:
        s.seq = Seq(cleaned_seq)  # 重要: 必须用Seq对象!
        clean_sequences.append(s)
```

##### 4.2.3 分批提交配置

```python
MAX_CONCURRENT_JOBS = 5   # 同时运行的Job数
batch_size = 20           # 每批序列数

# 申请的分析类型
appl = "PfamA,SuperFamily,CDD,Phobius,TMHMM,SignalP_EUK"
```

##### 4.2.4 Job状态轮询

```python
def process_batch(batch_data):
    job_id = submit_batch(batch_seqs, idx, email)
    
    # 等待完成 (支持所有等待状态)
    while status in ["RUNNING", "QUEUED", "STARTED", "PENDING"]:
        time.sleep(10)
        status = check_status(job_id)
        
        if wait_count > 180:  # 最长等待30分钟
            break
```

##### 4.2.5 断点续传 (Checkpoint)

```python
# 每批完成后立即保存
batch_file = os.path.join(temp_dir, f"batch_{idx}.tsv")

if os.path.exists(batch_file) and os.path.getsize(batch_file) > 0:
    # 跳过已完成的批次
    print(f"Batch {idx} already exists. Skipping.")
    return cached_content
```

---

## Phase 5: 序列提取

### 脚本: `universal_family_extractor.py` (802行)

#### 5.1 ID归一化

```python
def normalize_id(pid):
    """
    将各平台ID转换为统一格式:
    - XP_028218932.1 → XP_028218932_1
    - cds-CAA33989.1 → cds-CAA33989_1
    - EVM0001234.1 → EVM0001234_1
    
    规则: 将 . 和 _ 统一转为 _
    """
    return re.sub(r'[._]', '_', pid.split()[0])
```

#### 5.2 ⚠️ CDS映射策略 (7种策略)

```python
def build_cds_map(cds_path, proteome_path=None):
    """
    构建 Protein ID → CDS Record ID 的映射
    """
    
    # === Step 1: Ensembl特殊处理 (从蛋白文件提取transcript链接) ===
    if proteome_path:
        for record in SeqIO.parse(proteome_path, "fasta"):
            # Ensembl蛋白: >cds-CAA33989.1 pep ... transcript:transcript-rps2 ...
            m = re.search(r'transcript:(\S+)', desc)
            if m:
                prot_to_transcript[pid] = m.group(1)
    
    # === Step 2: 从CDS文件构建映射 ===
    for record in SeqIO.parse(cds_path, "fasta"):
        candidates = []
        
        # Strategy 1: NCBI [protein_id=XP_123.1]
        m = re.search(r'\[protein_id=([^\]]+)\]', desc)
        if m: candidates.append(m.group(1))
        
        # Strategy 2: NCBI lcl|..._prot_XP_123.1
        if '_prot_' in rid:
            m2 = re.search(r'([XY]P_\d+\.\d+)', sub)
            
        # Strategy 3: NCBI lcl|..._cds_XP_123.1
        if '_cds_' in rid:
            m3 = re.search(r'([XY]P_\d+[._]\d+)', sub)
        
        # Strategy 4: Ensembl gene:xxx
        m4 = re.search(r'gene:(\S+)', desc)
        
        # Strategy 5: EVM gene=xxx
        m5 = re.search(r'gene=(\S+)', desc)
        
        # Strategy 6: 直接使用record ID (简单格式如连翘)
        candidates.append(rid)
        
        # Strategy 7: 第一个token
        first_token = desc.split()[0]
        candidates.append(first_token)
    
    # === Step 3: 通过transcript链接Ensembl ===
    for prot_id, transcript_id in prot_to_transcript.items():
        if transcript_id in cds_records:
            mapping[prot_id] = cds_records[transcript_id]
```

#### 5.3 最长异构体过滤

```python
def filter_longest_isoforms(candidate_ids, proteome_path, gff_path):
    """
    按Gene分组, 每组保留最长转录本
    
    步骤:
    1. 解析GFF获取 mRNA→Gene 关系
    2. 获取每个候选的蛋白长度
    3. 按Gene ID分组
    4. 每组选择最长的
    """
    # 构建 mRNA → Gene 映射
    mrna_to_gene = {}
    for line in gff:
        if ftype == "mRNA":
            parent = attrs.get("Parent", "")
            mrna_to_gene[mrna_id] = parent
    
    # 分组
    gene_to_candidates = defaultdict(list)
    for cid in candidate_ids:
        gene = resolve_gene(cid)
        gene_to_candidates[gene].append((cid, length))
    
    # 选最长
    final_ids = []
    for gene, candidates in gene_to_candidates.items():
        candidates.sort(key=lambda x: x[1], reverse=True)
        final_ids.append(candidates[0][0])
```

#### 5.4 输出文件

| 文件 | 内容 | 提取逻辑 |
|------|------|----------|
| `family_members.pep.fasta` | 蛋白序列 | 直接从proteome提取 |
| `family_members.cds.fasta` | CDS序列 | 通过ID映射提取 |
| `family_members.gene.fasta` | 基因组序列 | 从GFF获取坐标 + 从genome提取 |
| `family_members.promoter.fasta` | 启动子序列 | TSS上游2000bp |
| `family_members.gff3` | 基因结构 | 从原始GFF过滤 |

---

## Phase 6: 多序列比对与Motif分析

### 6.1 MAFFT比对

```python
# pipeline_utils.run_mafft_alignment()
cmd = ["mafft", "--auto", "--amino", "--thread", str(threads)]
# 输入: family_members.pep.fasta
# 输出: family_members.aln.fasta
```

### 6.2 ⚠️ MSA Gap修剪 (系统发育树前)

```python
def trim_msa_by_gap(input_aln, output_aln, max_gap_ratio=0.5):
    """
    移除Gap比例>50%的列 (比HMM构建更严格)
    
    目的:
    - 提高系统发育树拓扑质量
    - 减少比对噪音
    
    输出: family_members.trimmed.aln
    """
    for i in range(length):
        col = alignment[:, i]
        gap_count = col.count("-")
        gap_ratio = gap_count / len(col)
        
        if gap_ratio <= max_gap_ratio:  # <=50%
            keep_cols.append(i)
```

### 6.3 Motif分析 (两种策略)

#### 策略1: MEME (Gold Standard)

```python
# 优先尝试本地MEME
cmd = [
    "meme", input_fasta,
    "-o", out_dir,
    "-protein",
    "-nmotifs", str(n_motifs),  # 默认15
    "-minw", "6",               # 最小motif宽度
    "-maxw", "50",              # 最大motif宽度
    "-mod", "zoops"             # Zero Or One Per Sequence
]

# Docker备选
if not local_meme:
    docker run memesuite/memesuite meme ...
```

#### 策略2: MSA-Based Fallback

```python
def extract_motifs_from_msa(msa_file, min_len=6, conservation_threshold=0.7):
    """
    从MSA直接提取保守区块
    
    步骤:
    1. 计算每列保守性 (最高频率氨基酸占比)
    2. 保守性 >= 70% 的列标记为保守
    3. 连续保守列组成区块
    4. 合并距离<=10的区块
    5. 过滤长度<6的区块
    """
    # 保守性计算
    for i in range(aln_len):
        col_nogap = [c for c in col if c != '-']
        top_char, count = Counter(col_nogap).most_common(1)[0]
        freq = count / num_seqs
        conserved_cols.append(freq >= conservation_threshold)
    
    # 区块合并
    for next_m in motifs_ranges[1:]:
        dist = next_m[0] - curr_m[1]
        if dist <= 10:  # 合并距离阈值
            curr_m = (curr_m[0], next_m[1])
```

---

## Phase 7: 系统发育树构建

### 7.1 工具优先级

```python
# run_pipeline_v3.py
iqtree_bin = shutil.which("iqtree2") or shutil.which("iqtree")

if iqtree_bin:
    run_iqtree(tree_input_aln, prefix, threads=cpu)
elif shutil.which("fasttree"):
    run_fasttree(tree_input_aln, tree_file, threads=cpu)
```

### 7.2 IQ-TREE参数

```python
def run_iqtree(aln_file, out_prefix, threads=4):
    cmd = [
        iqtree_bin,
        "-s", aln_file,
        "-m", "MFP",       # ModelFinder Plus: 自动模型选择
        "-bb", "1000",     # 1000次UFBoot
        "-alrt", "1000",   # 1000次SH-aLRT
        "-nt", str(threads),
        "-redo"            # 覆盖已有结果
    ]
```

### 7.3 FastTree参数

```bash
fasttree -lg < trimmed.aln > tree.nwk
# -lg: LG蛋白模型
```

---

## Phase 8: Ka/Ks选择压力分析

### 脚本: `run_kaks_analysis.py` (386行)

#### 8.1 完整流程

```
蛋白MSA (family_members.aln.fasta)
    │
    ▼ + CDS (family_members.cds.fasta)
┌─────────────────────────────────────────────────────┐
│ Step 1: 密码子比对 (protein_to_codon_alignment)     │
│   蛋白Gap → '---'                                   │
│   蛋白残基 → 对应密码子                              │
│   移除: 终止密码子, N碱基, Gap列                     │
└─────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────┐
│ Step 2: AXT格式转换 (write_axt_format)              │
│   生成所有配对组合 (All-vs-All)                      │
└─────────────────────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────────────────────┐
│ Step 3: KaKs_Calculator                             │
│   方法: MA (Model Averaging) 默认                   │
│   支持: YN (Yang-Nielsen), ML (Maximum Likelihood)  │
└─────────────────────────────────────────────────────┘
    │
    ▼
可视化 (直方图 + 散点图)
```

#### 8.2 ⚠️ 密码子比对函数

```python
def protein_to_codon_alignment(protein_aln_file, cds_file, output_file, remove_gaps=True):
    """
    将蛋白MSA"翻译"回密码子比对
    
    关键检查:
    1. 长度一致性: |CDS_len - Prot_len*3| <= 3
    2. 移除: 终止密码子 (TAA, TAG, TGA)
    3. 移除: 包含N的密码子
    4. 移除: Gap列 (如果remove_gaps=True)
    """
    # 映射逻辑
    for aa in pseq:
        if aa == "-":
            mapped_codons.append("---")
        else:
            codon = cseq[cds_idx : cds_idx+3]
            mapped_codons.append(codon)
            cds_idx += 3
    
    # 过滤列
    stops = {"TAA", "TAG", "TGA"}
    for i in range(aln_len):
        for rec in codon_aln_records:
            codon = rec['codons'][i]
            
            # 检查Gap
            if remove_gaps and codon == '---':
                keep = False
            
            # 检查N碱基
            if 'N' in codon:
                keep = False
            
            # 检查终止密码子
            if codon.upper() in stops:
                keep = False
```

#### 8.3 AXT格式转换

```python
def write_axt_format(codon_aln_file, output_axt):
    """
    生成All-vs-All配对的AXT格式
    
    格式:
    Seq1&Seq2
    ATGCATGC...
    ATGCATGC...
    
    (空行分隔)
    """
    from itertools import combinations
    
    for rec1, rec2 in combinations(records, 2):
        f.write(f"{rec1.id}&{rec2.id}\n")
        f.write(f"{str(rec1.seq)}\n")
        f.write(f"{str(rec2.seq)}\n")
        f.write("\n")
```

#### 8.4 并行KaKs计算

```python
# 分块策略
if n_threads > 1:
    chunk_size = math.ceil(total_pairs / n_threads)
    
    for i in range(n_threads):
        chunk_pairs = pairs[start:end]
        
        # 写入临时AXT
        with open(chunk_axt, 'w') as f:
            for block in chunk_pairs:
                f.writelines(block)
        
        # 提交计算
        cmd = ["KaKs_Calculator", "-i", chunk_axt, "-o", chunk_out, "-m", "MA"]
        executor.submit(run_cmd, cmd)

# 合并结果 (保留一个header)
with open(kaks_out_file, 'w') as outfile:
    header_written = False
    for chunk_out in chunk_files:
        lines = f.readlines()
        if not header_written:
            outfile.write(lines[0])
            header_written = True
        outfile.writelines(lines[1:])
```

---

## Phase 9: 共线性分析

### 脚本: `run_synteny_analysis.py` (819行)

#### 9.1 SyntenyAnalyzer类概览

```python
class SyntenyAnalyzer:
    def __init__(self, genome, gff, out_dir, prefix="species", 
                 threads=4, evalue="1e-5", cds=None, pep=None, highlights=None):
        """
        参数:
        - genome: 基因组FASTA文件
        - gff: GFF3注释文件
        - cds/pep: 可选的CDS/蛋白序列 (否则用gffread提取)
        - highlights: 高亮基因列表 (家族成员ID)
        """
        self.work_dir = os.path.join(self.out_dir, "Synteny_Work")
        
        # 输出文件路径
        self.bed_file = f"{prefix}.bed"
        self.blast_file = f"{prefix}.{prefix}.last"
        self.anchors_file = f"{prefix}.{prefix}.anchors"
```

#### 9.2 Step 1: 数据准备 (prepare_data)

##### 9.2.1 GFF→BED转换 (双遍扫描策略)

```python
def _fallback_gff_to_bed(self):
    """
    Pass 1: 构建完整映射
    """
    parent_map = {}       # FeatureID → ParentID
    gene_attr_map = {}    # FeatureID → GeneName
    feat_to_prot = {}     # FeatureID → ProteinID
    feat_coords = {}      # FeatureID → (chrom, start, end, strand)
    
    # 需要提取的属性键
    target_keys = {'protein_id', 'id', 'name', 'alias', 'product', 'transcript_id'}
    gene_keys = {'gene', 'gene_id', 'locus_tag', 'name', 'gene_name'}
    
    for line in gff_file:
        # 解析GFF9列格式
        parts = line.split('\t')
        ftype = parts[2]  # gene, mRNA, CDS等
        
        for chunk in attr.split(';'):
            k, v = chunk.split('=', 1)
            
            if k == 'ID': feat_id = v
            elif k == 'Parent': parent = v
            elif k in gene_keys: found_gene_attrs[k] = v
            elif k in target_keys: potential_prot_ids.append(v)
            elif k == 'Dbxref':
                # 解析 Dbxref=GeneID:xxx,Genbank:xxx
                for ref in v.split(','):
                    potential_prot_ids.append(ref.split(':')[-1])
    
    """
    Pass 2: 输出BED
    """
    target_types = {'mrna', 'transcript', 'cds'}
    
    for feature in gff_lines:
        if ftype not in target_types: continue
        
        # 解析ProteinID (优先级)
        prot_id = feat_to_prot.get(feat_id)
        if not prot_id and parent:
            prot_id = feat_to_prot.get(parent)
        
        # BED格式 (0-based start)
        bed_start = start - 1
        fout.write(f"{chrom}\t{bed_start}\t{end}\t{prot_id}\t1000\t{strand}\n")
```

##### 9.2.2 FASTA清洗

```python
def _clean_fasta(self, fasta_path):
    """
    JCVI/Diamond兼容性处理:
    1. 仅保留ID (移除description)
    2. 移除终止密码子 (*)
    3. 移除Gap字符 (.)
    """
    for record in SeqIO.parse(fin, "fasta"):
        clean_id = record.id.split()[0]
        clean_seq = str(record.seq).replace('.', '').replace('*', '')
        fout.write(f">{clean_id}\n{clean_seq}\n")
```

#### 9.3 Step 2: 同源搜索 (run_homology_search)

```python
def run_homology_search(self):
    # 工具优先级: Diamond > BLAST+
    
    if check_dependency("diamond"):
        # Diamond (比BLAST快100x)
        db_file = f"{prefix}.dmnd"
        
        # 1. 构建索引
        cmd_makedb = f"diamond makedb --in {pep_file} --db {db_file}"
        
        # 2. All-vs-All搜索
        # 输出格式: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        cmd_blastp = f"""diamond blastp 
            -q {pep_file} 
            --db {db_file} 
            -o {blast_file} 
            --evalue {evalue} 
            --threads {threads}"""
    
    elif check_dependency("blastp"):
        # BLAST+ (较慢但更准确)
        cmd_makedb = f"makeblastdb -in {pep_file} -dbtype prot -out {db_base}"
        cmd_blastp = f"blastp -query {pep_file} -db {db_base} -out {blast_file} -evalue {evalue} -num_threads {threads} -outfmt 6"
```

#### 9.4 Step 3: MCScanX共线性检测 (run_synteny)

```python
def run_synteny(self):
    # 切换到工作目录 (JCVI要求)
    os.chdir(self.work_dir)
    
    # JCVI MCScanX命令
    # --cscore=.7: 共线性评分阈值
    # --no_strip_names: 保留原始ID (重要!)
    cmd = f"python -m jcvi.compara.catalog ortholog {prefix} {prefix} --cscore=.7 --no_strip_names"
    
    # 输出: {prefix}.{prefix}.anchors (共线性基因对)
```

#### 9.5 Step 4: Circos可视化输入生成 (generate_circos_conf)

```python
def generate_circos_conf(self):
    # 1. 解析BED获取染色体尺寸和基因坐标
    for line in bed_file:
        chrom, start, end, name = parts[0:4]
        gene_coords[name] = (chrom, start, end)
        chrom_sizes[chrom] = max(chrom_sizes[chrom], end)
    
    # 2. 加载高亮基因 (家族成员)
    if self.highlights_file:
        for line in f:
            raw_id = line.strip()
            # ID归一化处理
            if '_' in raw_id:
                parts = raw_id.rsplit('_', 1)
                if parts[1].isdigit():
                    normalized_id = f"{parts[0]}.{parts[1]}"  # XP_xxx_1 → XP_xxx.1
                    highlight_ids.add(normalized_id)
    
    # 3. 染色体选择算法
    MAX_DISPLAY = 12
    
    # 优先保留包含目标基因的染色体
    target_chroms = set()
    for gid in highlight_ids:
        if gid in gene_coords:
            target_chroms.add(gene_coords[gid][0])
    
    # 预扫描anchors确定有高亮连接的染色体
    for line in anchors_file:
        gene_a, gene_b = parts[0], parts[1]
        if gene_a in highlight_ids or gene_b in highlight_ids:
            chroms_with_highlight_links.add(gene_coords[gene_a][0])
            chroms_with_highlight_links.add(gene_coords[gene_b][0])
    
    # 如果目标染色体 < 12, 用最长染色体填充
    if len(target_chroms) < MAX_DISPLAY:
        sorted_by_len = sorted(all_chroms, key=lambda k: chrom_sizes[k], reverse=True)
        for chrom in sorted_by_len:
            if chrom not in final_chrom_set:
                final_chrom_set.add(chrom)
                if len(final_chrom_set) >= MAX_DISPLAY: break
    
    # 4. 生成Karyotype文件
    # 格式: chr - chr_name label start end color
    for idx, chrom in enumerate(final_sorted):
        color = f"chr{idx % 20 + 1}"
        f.write(f"chr - {chrom} {chrom} 0 {chrom_sizes[chrom]} {color}\n")
    
    # 5. 生成Links文件
    for line in anchors_file:
        gene_a, gene_b = parts[0], parts[1]
        chr_a, start_a, end_a = gene_coords[gene_a]
        chr_b, start_b, end_b = gene_coords[gene_b]
        
        # 仅绘制选定染色体间的连线
        if chr_a in top_chrom_set and chr_b in top_chrom_set:
            # 高亮连线 (红色加粗)
            if gene_a in highlight_ids or gene_b in highlight_ids:
                link_opts = "color=red,thickness=4,z=10"
            else:
                link_opts = "color=grey_a4"
            
            fout.write(f"{chr_a} {start_a} {end_a} {chr_b} {start_b} {end_b} {gene_a} {gene_b} {link_opts}\n")
```

#### 9.6 Step 5: Matplotlib Circos绘图 (render_circos_plot)

```python
def _render_cartesian_circos(self, karyotype_path, links_path, output_path):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    
    # 1. 布局计算
    # 将染色体排列成圆形
    total_genome_len = sum(chrom_sizes.values())
    gap_angle = 2  # 染色体间隔角度
    
    # 2. 绘制染色体弧段
    for chrom, length in chroms:
        arc = patches.Arc(center, radius, radius, 
                         theta1=start_angle, theta2=end_angle)
        ax.add_patch(arc)
    
    # 3. 绘制共线性连线 (贝塞尔曲线)
    for link in links:
        # 计算起点和终点的角度位置
        angle_a = get_angle(chr_a, pos_a)
        angle_b = get_angle(chr_b, pos_b)
        
        # 通过中心点绘制贝塞尔曲线
        Path = mpath.Path
        path_data = [
            (Path.MOVETO, point_a),
            (Path.CURVE3, control_point),
            (Path.CURVE3, point_b)
        ]
        
        # 高亮连线使用红色
        if is_highlight:
            color = 'red'
            linewidth = 1.5
        else:
            color = 'gray'
            alpha = 0.3
    
    # 4. 保存图片
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
```

---

## Phase 10: 启动子分析

### 脚本: `run_promoter_analysis.py` (147行)

#### 10.1 启动子文件定位

```python
def main():
    # 优先使用用户提供的启动子文件
    if args.promoter_file:
        promoter_file = Path(args.promoter_file)
    else:
        # 尝试从pipeline输出目录查找
        candidates = list(final_dir.glob("*.promoter.fasta"))
        if candidates:
            promoter_file = candidates[0]
    
    # 检查文件有效性
    if promoter_file.stat().st_size == 0:
        logger.warning("Promoter file is empty. Skipping analysis.")
        sys.exit(0)
```

#### 10.2 Cis-element数据库定位 (多路径搜索)

```python
# 搜索顺序:
# 1. config/cis_elements_default.txt
# 2. ./cis_elements_default.txt (当前目录)
# 3. {script_dir}/cis_elements_default.txt (脚本目录)

cis_db_path = Path("config/cis_elements_default.txt")
if not cis_db_path.exists():
    cis_db_path = Path("cis_elements_default.txt")
if not cis_db_path.exists():
    script_dir = Path(__file__).parent.resolve()
    cis_db_path = script_dir / "cis_elements_default.txt"
```

#### 10.3 Golden List扫描 (正则匹配)

```python
# pipeline_utils.scan_promoters()
def scan_promoters(promoter_file, cis_db_path):
    """
    使用正则表达式扫描已知的顺式作用元件
    
    cis_elements_default.txt 格式:
    Element_Name\tSequence_Pattern
    ABRE\tACGTG[TC]
    W-box\tTTGAC[TC]
    MYB\tC[ACGT]GTT[AG]
    ...
    """
    cis_counts = {}  # Element → Count
    cis_details = {} # Gene → [(Element, Position, Sequence)]
    
    for gene_id, sequence in promoter_seqs:
        for element_name, pattern in cis_patterns:
            for match in re.finditer(pattern, sequence, re.IGNORECASE):
                cis_counts[element_name] += 1
                cis_details[gene_id].append({
                    'element': element_name,
                    'position': match.start(),
                    'sequence': match.group()
                })
    
    return cis_counts, cis_details
```

#### 10.4 De Novo Motif发现 (MEME)

```python
# pipeline_utils.run_meme_promoter()
def run_meme_promoter(promoter_file, out_dir, threads=1):
    """
    启动子专用MEME参数
    """
    cmd = [
        "meme", str(promoter_file),
        "-dna",              # DNA模式 (非蛋白)
        "-oc", str(out_dir),
        "-nmotifs", "10",    # 发现10个motif
        "-minw", "6",        # 最小宽度6bp
        "-maxw", "20",       # 最大宽度20bp
        "-mod", "zoops",     # Zero Or One Per Sequence
        "-p", str(threads)   # 并行线程
    ]
    
    subprocess.run(cmd, check=True)
```

#### 10.5 TOMTOM Motif验证

```python
# pipeline_utils.run_tomtom()
def run_tomtom(meme_out_dir, jaspar_db_path):
    """
    将MEME发现的motif与JASPAR数据库比对
    验证发现的motif是否对应已知转录因子结合位点
    """
    meme_file = meme_out_dir / "meme.txt"
    tomtom_out = meme_out_dir / "tomtom_results"
    
    cmd = [
        "tomtom",
        "-oc", str(tomtom_out),
        "-evalue", "10.0",      # 宽松阈值便于探索
        "-min-overlap", "5",    # 最小重叠5bp
        str(meme_file),
        str(jaspar_db_path)
    ]
    
    subprocess.run(cmd, check=True)

# JASPAR数据库文件:
# cis_elements_jaspar_plants.meme (植物转录因子MEME格式)
```

#### 10.6 结果输出

| 文件 | 内容 | 格式 |
|------|------|------|
| `cis_counts.json` | 元件计数 | `{"ABRE": 45, "W-box": 32}` |
| `cis_details.json` | 详细命中 | `{"Gene1": [{"element": "ABRE", "pos": 123}]}` |
| `MEME_Promoter/meme.txt` | MEME输出 | MEME标准格式 |
| `MEME_Promoter/tomtom_results/` | TOMTOM比对 | HTML + TSV |



---

## 关键技术细节备忘

### ⚠️ 必须注意的处理步骤

| 阶段 | 处理 | 参数 | 原因 |
|------|------|------|------|
| **种子获取** | 序列去重 | 基于序列内容 | 避免冗余影响HMM |
| **HMM构建前** | 长度过滤 | 75%-125%中位数 | 移除异常值 |
| **HMM构建前** | MSA Gap修剪 | >90%列移除 | 优化HMM profile |
| **InterPro提交前** | 序列清洗 | 移除`*`和`.` | 避免API错误 |
| **系统发育树前** | MSA Gap修剪 | >50%列移除 | 提高树拓扑质量 |
| **密码子比对后** | 终止密码子移除 | TAA/TAG/TGA | 避免KaKs错误 |
| **密码子比对后** | N碱基列移除 | 全部 | 确保比对有效 |
| **Motif分析** | 区块合并 | 距离<=10 | 连接相近保守区 |

### 🔧 重要配置参数参考

| 参数 | 默认值 | 脚本 | 说明 |
|------|--------|------|------|
| `evalue` | 1e-5 | search_extract | HMM/BLAST阈值 |
| `max_gap_ratio` (HMM) | 0.9 | build_hmm | HMM构建MSA修剪 |
| `max_gap_ratio` (Tree) | 0.5 | pipeline_utils | 系统发育MSA修剪 |
| `conservation_threshold` | 0.7 | pipeline_utils | Motif保守性阈值 |
| `merge_distance` | 10 | pipeline_utils | Motif区块合并距离 |
| `batch_size` (CDD) | 400 | scan_cdd_ncbi | CDD分批大小 |
| `batch_size` (IPS) | 20 | interproscan_runner | InterPro分批大小 |
| `MAX_CONCURRENT_JOBS` | 5 | interproscan_runner | InterPro并发数 |
| `upstream_len` | 2000 | run_promoter_analysis | 启动子长度 |
| `n_motifs` | 15 | pipeline_utils | MEME发现数 |
| `bootstrap` | 1000 | pipeline_utils | IQ-TREE bootstrap |

### 🔑 ID映射优先级

```
1. NCBI [protein_id=XP_xxx]          → 最高优先级
2. NCBI lcl|..._prot_XP_xxx          → 次优先级
3. NCBI lcl|..._cds_XP_xxx           → 次优先级
4. Ensembl transcript:xxx             → 新增(v3.1)
5. Ensembl gene:xxx                   → 中优先级
6. EVM gene=xxx                       → 中优先级
7. 直接使用 record ID                 → 兜底策略
8. 第一个 token                       → 备用策略
```

---

*文档结束 - 版本 v3.1*
