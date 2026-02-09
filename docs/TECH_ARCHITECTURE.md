# ChemVCS 技术架构草案

> v1.0 | 2026-02-09  
> 状态：Draft — 供工程评审  
> 范围：MVP（Phase 1）核心 + Phase 2 HPC 集成的架构预留

---

## 1. 设计目标与非目标

### 设计目标

| # | 目标 | 量化约束 |
|---|------|----------|
| D1 | **单进程、零服务** — 所有操作在一次 CLI 调用内完成，不启动后台进程或监听端口 | `lsof -i` 在 chemvcs 运行期间无新增 socket |
| D2 | **共享文件系统安全** — 在 Lustre / GPFS / NFS 上正确运行，即使多个用户在同一项目目录操作 | `metadata.db` 写操作串行化；blob 写入幂等 |
| D3 | **可离线构建、可离线运行** — 安装和运行均不需要网络访问 | 提供 `sdist` + `wheel`；运行时无 HTTP 调用 |
| D4 | **commit 操作 O(n) 于已暂存文件** — 不做全库扫描或索引重建 | 5 个 VASP 输入文件的 commit <2s（Lustre 冷缓存） |
| D5 | **数据完整性可校验** — 任意时刻可验证所有 blob 是否与元数据一致 | `chemvcs fsck` 遍历所有 blob，对比 hash，报告损坏文件 |

### 非目标

| # | 非目标 | 说明 |
|---|--------|------|
| N1 | 不实现分布式对象传输协议 | MVP 无 push/pull；远程同步将在 Phase 3 通过 rsync + manifest 实现 |
| N2 | 不实现 DAG 历史（分支/合并） | MVP 仅线性历史（单链表），避免三路合并的语义复杂度 |
| N3 | 不做 blob 层加密 | 对象存储明文写入；加密需求由文件系统或用户自行处理 |
| N4 | 不做增量 diff 存储（delta encoding） | blob 按全量存储 + 去重；delta 编码收益与复杂度不成比例（VASP 输入文件通常 <100 KB） |
| N5 | 不内建任务队列或调度器 | HPC 集成仅生成 sbatch 脚本并调用系统命令；不管理任务生命周期 |

---

## 2. 高层架构

```
┌───────────────────────────────────────────────────────────────────┐
│                          CLI 层 (Typer)                           │
│  init │ add │ commit │ log │ diff │ reproduce │ status │ fsck    │
└───┬───┴──┬──┴───┬────┴──┬──┴──┬───┴─────┬─────┴───┬────┴───┬────┘
    │      │      │       │     │         │         │        │
    ▼      ▼      ▼       ▼     ▼         ▼         ▼        ▼
┌───────────────────────────────────────────────────────────────────┐
│                       核心引擎层 (Core)                            │
│                                                                   │
│  ┌──────────┐  ┌──────────┐  ┌───────────┐  ┌────────────────┐  │
│  │ Staging  │  │ Commit   │  │ History   │  │ Reproduce      │  │
│  │ Manager  │  │ Builder  │  │ Walker    │  │ Engine         │  │
│  └────┬─────┘  └────┬─────┘  └─────┬─────┘  └───────┬────────┘  │
│       │             │              │                 │            │
│       ▼             ▼              ▼                 ▼            │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │                 存储抽象层 (Storage)                          │ │
│  │  ObjectStore    │   MetadataDB    │   Config               │ │
│  │  (blob r/w)     │   (SQLite)      │   (.json)              │ │
│  └────┬────────────┴───────┬─────────┴──────────────────────┘  │
│       │                    │                                     │
└───────┼────────────────────┼─────────────────────────────────────┘
        │                    │
        ▼                    ▼
┌──────────────┐   ┌─────────────────┐
│  .chemvcs/   │   │  .chemvcs/      │
│  objects/    │   │  metadata.db    │
│  (文件系统)   │   │  (SQLite file)  │
└──────────────┘   └─────────────────┘

    ┌───────────────────────────────────────────────────────────┐
    │                    语义解析层 (Parsers)                     │
    │                                                           │
    │  ┌──────────┐ ┌───────────┐ ┌──────────┐ ┌────────────┐ │
    │  │ INCAR    │ │ POSCAR    │ │ KPOINTS  │ │ OUTCAR     │ │
    │  │ Parser   │ │ Parser    │ │ Parser   │ │ Extractor  │ │
    │  └──────────┘ └───────────┘ └──────────┘ └────────────┘ │
    │  ┌──────────┐ ┌───────────┐                              │
    │  │ POTCAR   │ │ Generic   │   ← 非 VASP 文件退回文本 diff │
    │  │ RefOnly  │ │ TextDiff  │                              │
    │  └──────────┘ └───────────┘                              │
    └───────────────────────────────────────────────────────────┘

    ┌───────────────────────────────────────────────────────────┐
    │               HPC 集成层 (Phase 2 预留)                    │
    │                                                           │
    │  ┌──────────────┐  ┌──────────────┐  ┌────────────────┐  │
    │  │ SLURM        │  │ PBS          │  │ Environment    │  │
    │  │ Adapter      │  │ Adapter      │  │ Snapshot       │  │
    │  │ (sbatch/     │  │ (qsub/       │  │ (module list,  │  │
    │  │  squeue/     │  │  qstat/      │  │  VASP version) │  │
    │  │  scancel)    │  │  qdel)       │  │                │  │
    │  └──────────────┘  └──────────────┘  └────────────────┘  │
    └───────────────────────────────────────────────────────────┘

    ┌───────────────────────────────────────────────────────────┐
    │              插件层 (Phase 4 预留，MVP 不实现)               │
    │                                                           │
    │  chemvcs.parsers  entry_point group                       │
    │  → 第三方可注册 QE / CP2K / ORCA / LAMMPS 解析器           │
    └───────────────────────────────────────────────────────────┘
```

---

## 3. 存储设计

### 3.1 目录结构

```
.chemvcs/
├── VERSION                     # 存储格式版本号，纯文本 "1"
├── config.json                 # 项目级配置
├── HEAD                        # 当前 commit id（纯文本，40 hex chars）
├── metadata.db                 # SQLite — 可重建的索引/缓存
├── metadata.db-wal             # SQLite WAL（若适用）
├── metadata.db-shm             # SQLite shared memory（若适用）
├── objects/                    # content-addressable blob 存储（真相源）
│   ├── ab/
│   │   └── cdef0123456789...   # SHA-256 前 2 字符为目录
│   ├── fc/
│   │   └── 1a2b3c4d5e6f...
│   └── ...
├── commits/                    # commit 对象（JSON 文件，同样按 hash 寻址）
│   ├── de/
│   │   └── ad0000beef...json
│   └── ...
├── staging/
│   └── manifest.json           # 暂存区清单
├── tmp/                        # 写入中转区（同文件系统，确保 rename 原子性）
└── telemetry.jsonl             # 本地遥测日志（可选，默认关闭）
```

**关键不变量：**
- `objects/` 和 `commits/` 是 **append-only**，任何写入完成后不再修改
- `metadata.db` 是 **可重建的缓存**——可从 `commits/` + `objects/` 完全重建
- `HEAD` 文件始终指向最新 commit hash

### 3.2 Commit 对象格式

每个 commit 以 JSON 文件存储在 `commits/<hash前2>/<hash剩余>.json`。commit hash = SHA-256(JSON 内容的 canonical 序列化)。

```jsonc
// .chemvcs/commits/de/ad0000beef....json
{
  "chemvcs_version": "0.1.0",
  "format_version": 1,
  "id": "dead0000beef1234567890abcdef1234567890abcdef1234567890abcdef1234",
  "parent_id": "aabb0000ccdd...",   // null for first commit
  "timestamp": "2026-02-09T14:32:01.123456Z",
  "author": "liming@hpc-login01",
  "message": "PBE+U, U=3.5eV, kpoints 4x4x4",

  "files": [
    {
      "path": "INCAR",
      "blob_hash": "fc1a2b3c...",
      "size_bytes": 842,
      "file_type": "INCAR",
      "is_reference": false
    },
    {
      "path": "POSCAR",
      "blob_hash": "ab12cd34...",
      "size_bytes": 1204,
      "file_type": "POSCAR",
      "is_reference": false
    },
    {
      "path": "KPOINTS",
      "blob_hash": "ee44ff55...",
      "size_bytes": 56,
      "file_type": "KPOINTS",
      "is_reference": false
    },
    {
      "path": "POTCAR",
      "blob_hash": "0000aaaa...",  // 指向引用 blob（JSON 元信息）
      "size_bytes": 0,             // 原始文件大小记录为 0（不存内容）
      "file_type": "POTCAR_REF",
      "is_reference": true
    }
  ],

  "semantic_summary": {
    "incar": {
      "ENCUT": 520,
      "ISMEAR": 0,
      "SIGMA": 0.05,
      "ISPIN": 2,
      "LDAU": true,
      "LDAUU": [3.5, 0, 0],
      "LDAUJ": [0, 0, 0],
      "EDIFF": 1e-6,
      "NSW": 100,
      "IBRION": 2
    },
    "poscar": {
      "formula": "Li4Co3FeO8",
      "spacegroup": "R-3m",
      "natoms": 16,
      "lattice": {
        "a": 2.830, "b": 2.830, "c": 14.050,
        "alpha": 90.0, "beta": 90.0, "gamma": 120.0
      },
      "volume_A3": 97.45
    },
    "kpoints": {
      "mode": "Gamma",
      "grid": [4, 4, 4],
      "shift": [0, 0, 0]
    },
    "potcar": [
      {"element": "Li", "functional": "PBE", "hash": "a1b2c3d4...", "path": "/opt/vasp/pp/PBE/Li_sv/POTCAR"},
      {"element": "Co", "functional": "PBE", "hash": "e5f6a7b8...", "path": "/opt/vasp/pp/PBE/Co/POTCAR"},
      {"element": "O",  "functional": "PBE", "hash": "c9d0e1f2...", "path": "/opt/vasp/pp/PBE/O/POTCAR"}
    ]
  },

  "output_summary": null,  // 仅当 commit 含 OUTCAR 时填充

  "environment": {
    "hostname": "hpc-login01",
    "python_version": "3.9.7",
    "chemvcs_version": "0.1.0",
    "vasp_version": null,
    "loaded_modules": null,
    "slurm_job_id": null
  }
}
```

当 commit 包含 OUTCAR/vasprun.xml 时，`output_summary` 填充：

```jsonc
"output_summary": {
  "total_energy_eV": -45.2314,
  "energy_per_atom_eV": -2.8270,
  "is_converged": true,
  "ionic_steps": 23,
  "max_force_eV_A": 0.0012,
  "bandgap_eV": 1.82,
  "magnetic_moment_uB": 3.01,
  "warnings": [],
  "vasp_version": "6.3.2",
  "elapsed_time_s": 3421.5
}
```

### 3.3 索引 / 缓存策略：SQLite 的角色

**核心原则：SQLite 是可重建的索引，不是真相源。**

```
真相源 (Source of Truth)          缓存 (Rebuildable Index)
──────────────────────          ─────────────────────────
.chemvcs/commits/*.json    →    metadata.db: commits 表
.chemvcs/objects/*         →    metadata.db: files 表
commit JSON 内嵌字段       →    metadata.db: semantic_summary 表
                                metadata.db: output_summary 表
                                metadata.db: environment 表
.chemvcs/HEAD              →    metadata.db: head 表
```

**重建流程（`chemvcs rebuild-index`）：**

```python
def rebuild_index(repo_path):
    """从 commits/ 目录重建 metadata.db，时间复杂度 O(commits)"""
    db = init_empty_db(repo_path / ".chemvcs" / "metadata.db")
    for commit_file in walk_commits_dir(repo_path / ".chemvcs" / "commits"):
        commit = json.loads(commit_file.read_text())
        db.insert_commit(commit)
        for f in commit["files"]:
            db.insert_file(commit["id"], f)
        if commit.get("semantic_summary"):
            db.insert_semantic(commit["id"], commit["semantic_summary"])
        if commit.get("output_summary"):
            db.insert_output(commit["id"], commit["output_summary"])
    db.set_head(read_file(repo_path / ".chemvcs" / "HEAD"))
```

**共享文件系统上的 SQLite 策略：**

```python
import sqlite3, os, fcntl

def open_db(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(str(db_path), timeout=10)
    
    # 1. 尝试 WAL 模式（性能最好，读写并发）
    try:
        conn.execute("PRAGMA journal_mode=WAL")
        result = conn.execute("PRAGMA journal_mode").fetchone()[0]
        if result.upper() != "WAL":
            raise RuntimeError("WAL not supported")
    except Exception:
        # 2. 回退到 DELETE 模式（NFS/部分 Lustre 兼容）
        conn.execute("PRAGMA journal_mode=DELETE")
    
    conn.execute("PRAGMA busy_timeout=5000")  # 等待锁最多 5s
    conn.execute("PRAGMA synchronous=NORMAL")  # 平衡性能与安全
    return conn

class DBWriteLock:
    """外层文件锁——防止多进程同时写 metadata.db"""
    def __init__(self, lock_path: Path):
        self.lock_path = lock_path
    
    def __enter__(self):
        self.fd = open(self.lock_path, "w")
        fcntl.flock(self.fd, fcntl.LOCK_EX)  # 阻塞等待
        return self
    
    def __exit__(self, *args):
        fcntl.flock(self.fd, fcntl.LOCK_UN)
        self.fd.close()
```

**SQLite 可选替换路径（当 SQLite 完全不可用时）：**

如果检测到文件系统不支持 SQLite 的任何日志模式（极端情况），退回到 **纯 JSON 索引**：

```
.chemvcs/index.json   # 所有 commit 的精简列表，append-only
```

此路径性能下降（log/grep 需要扫描全文件），但保证可用性。在 MVP 中作为 fallback 实现，不作为默认路径。

### 3.4 Append-only Log

`commits/` 目录本身就是 append-only log：

- 每个 commit 是一个独立的 JSON 文件
- commit 通过 `parent_id` 形成单向链表
- `HEAD` 文件指向链尾
- **从不修改已提交的 commit 文件**（不可变性保证完整性校验可行）
- 读取历史 = 从 HEAD 沿 parent_id 回溯

```
HEAD → commit_C.json → commit_B.json → commit_A.json → null
            ↑ parent_id   ↑ parent_id    ↑ parent_id
```

---

## 4. 语义 Diff 设计

### 4.1 总体架构

```python
# chemvcs/diff/engine.py

class DiffEngine:
    """根据 file_type 分派到对应的语义 diff 实现"""
    
    REGISTRY: dict[str, type["BaseDiffer"]] = {}
    
    @classmethod
    def register(cls, file_type: str):
        def decorator(differ_cls):
            cls.REGISTRY[file_type] = differ_cls
            return differ_cls
        return decorator
    
    @classmethod
    def diff(cls, file_type: str, old_content: bytes, new_content: bytes,
             old_parsed: dict | None, new_parsed: dict | None) -> "DiffResult":
        differ = cls.REGISTRY.get(file_type, GenericTextDiffer)
        return differ().compute(old_content, new_content, old_parsed, new_parsed)

class DiffResult:
    file_type: str
    changes: list[ChangeItem]       # 结构化变更列表
    summary_text: str               # 人类可读摘要
    raw_json: dict                  # 机器可读 JSON（--format json 输出）
```

### 4.2 INCAR Diff — key-value 语义比较

```python
@DiffEngine.register("INCAR")
class IncarDiffer(BaseDiffer):
    """
    策略：
    1. 用 pymatgen.io.vasp.Incar.from_string() 解析为 dict
    2. 失败时退回正则解析 key = value
    3. 逐 key 比较：新增 / 删除 / 修改
    4. 数值型参数计算差值和单位
    """
    
    # 已知参数的单位与类型映射
    PARAM_META = {
        "ENCUT":  {"type": "float", "unit": "eV"},
        "SIGMA":  {"type": "float", "unit": "eV"},
        "EDIFF":  {"type": "float", "unit": "eV"},
        "EDIFFG": {"type": "float", "unit": "eV/Å"},
        "NSW":    {"type": "int",   "unit": "steps"},
        "ISMEAR": {"type": "int",   "unit": ""},
        "ISPIN":  {"type": "int",   "unit": ""},
        "LDAUU":  {"type": "float_list", "unit": "eV"},
        "LDAUJ":  {"type": "float_list", "unit": "eV"},
        "MAGMOM": {"type": "float_list", "unit": "μB"},
        # ... 约 30 个常用参数
    }
```

**输出示例：**

```
── INCAR diff ──────────────────────────────────
  MODIFIED  LDAUU    : 3.5 0 0 → 4.0 0 0  (+0.5 eV on site 1)
  MODIFIED  EDIFF    : 1e-05 → 1e-06  (tighter by 10×)
  ADDED     LWAVE    : .TRUE.
  DELETED   LCHARG   : (was .FALSE.)
  UNCHANGED : ENCUT=520, ISMEAR=0, SIGMA=0.05, ... (18 params)
```

**类型归一化规则：**

| INCAR 写法 | 归一化后 | 比较方式 |
|------------|----------|----------|
| `.TRUE.` / `T` / `.T.` | `True` (bool) | 布尔相等 |
| `3*0.0` | `[0.0, 0.0, 0.0]` | 列表逐元素比较 |
| `520` / `520.0` / `520.` | `520.0` (float) | 浮点相等（tolerance=1e-10） |
| `# comment` / `!comment` | 忽略 | 不参与 diff |

### 4.3 POSCAR Diff — 结构差异

```python
@DiffEngine.register("POSCAR")
class PoscarDiffer(BaseDiffer):
    """
    策略：
    1. 用 pymatgen.core.Structure.from_str(content, fmt="poscar") 解析
    2. 比较晶格参数（a, b, c, α, β, γ）
    3. 比较化学式和原子数
    4. 比较空间群（SpacegroupAnalyzer, symprec=0.1）
    5. 若原子数相同，计算坐标 RMSD（StructureMatcher）
    6. 若原子数不同，报告增删的原子种类与数量
    """
```

**输出示例：**

```
── POSCAR diff ─────────────────────────────────
  Formula    : Li4Co4O8 → Li4Co3FeO8  (substitution: Co→Fe ×1)
  Spacegroup : R-3m → R-3m  (unchanged)
  Lattice a  : 2.830 → 2.845 Å  (+0.015)
  Lattice c  : 14.050 → 14.120 Å  (+0.070)
  Volume     : 97.45 → 98.12 Å³  (+0.69, +0.7%)
  Natoms     : 16 → 16  (unchanged)
  Coord RMSD : 0.0234 Å  (over 16 matched sites)
```

### 4.4 KPOINTS Diff — 网格/路径差异

```python
@DiffEngine.register("KPOINTS")
class KpointsDiffer(BaseDiffer):
    """
    策略：
    1. 用 pymatgen.io.vasp.Kpoints.from_string() 解析
    2. 自动识别三种模式：
       a. Automatic (Gamma/Monkhorst-Pack) → 比较网格尺寸和偏移
       b. Line-mode → 比较高对称点路径和分点数
       c. Explicit → 比较 k-point 总数
    """
```

**输出示例（网格模式）：**

```
── KPOINTS diff ────────────────────────────────
  Mode   : Gamma → Gamma  (unchanged)
  Grid   : 4×4×4 → 6×6×6  (+50% density per axis)
  Shift  : 0 0 0 → 0 0 0  (unchanged)
  Est. k : 64 → 216 points  (×3.4)
```

**输出示例（路径模式）：**

```
── KPOINTS diff ────────────────────────────────
  Mode      : Line → Line
  Divisions : 20 → 40  (doubled)
  Path      : Γ-X-M-Γ → Γ-X-M-Γ-Z  (added segment Γ-Z)
```

### 4.5 OUTCAR / vasprun.xml — 输出摘要对比

不做全文 diff。仅比较 `output_summary` 中的结构化字段。

```
── Output diff ─────────────────────────────────
  Energy     : -45.2314 → -44.8732 eV  (+0.3582 eV, +0.8%)
  Per atom   : -2.8270 → -2.8046 eV/atom
  Converged  : ✅ → ✅
  Ionic steps: 23 → 31  (+8)
  Max force  : 0.0012 → 0.0008 eV/Å  (improved)
  Bandgap    : 1.82 → 1.95 eV  (+0.13)
  Mag moment : 3.01 → 2.98 μB  (-0.03)
  Warnings   : [] → ["BRMIX: very serious problems"]  ⚠️ NEW WARNING
```

### 4.6 退回策略

对于 `file_type == "OTHER"` 或解析失败的文件，使用标准 unified diff：

```python
@DiffEngine.register("OTHER")
class GenericTextDiffer(BaseDiffer):
    def compute(self, old: bytes, new: bytes, **_) -> DiffResult:
        # 二进制文件检测（含 \x00 字节）→ 仅报告 size 和 hash 变化
        if b"\x00" in old[:8192] or b"\x00" in new[:8192]:
            return DiffResult(summary=f"Binary: {len(old)}B → {len(new)}B")
        # 文本文件 → unified diff
        return DiffResult(summary=unified_diff(old.decode(), new.decode()))
```

---

## 5. 大文件策略

### 5.1 文件分类与追踪决策

```
                          默认行为                    add --force 行为
                  ┌─────────────────────┐     ┌───────────────────────┐
  INCAR           │ ✅ 追踪，存储 blob   │     │  —                    │
  POSCAR          │ ✅ 追踪，存储 blob   │     │  —                    │
  KPOINTS         │ ✅ 追踪，存储 blob   │     │  —                    │
  POTCAR          │ ⚠️ 引用模式 (hash)  │     │  阻断，始终引用        │
  OUTCAR          │ ✅ 追踪，存储 blob   │     │  —                    │
  vasprun.xml     │ ✅ 追踪，blob+gzip  │     │  —                    │
  ──────────────  │─────────────────────│     │───────────────────────│
  WAVECAR         │ ❌ .chemvcsignore   │     │  📌 仅记录 hash+size  │
  CHGCAR          │ ❌ .chemvcsignore   │     │  📌 仅记录 hash+size  │
  CHG             │ ❌ .chemvcsignore   │     │  📌 仅记录 hash+size  │
  PROCAR          │ ❌ .chemvcsignore   │     │  📌 仅记录 hash+size  │
                  └─────────────────────┘     └───────────────────────┘
```

### 5.2 引用指针文件

当 `is_reference = true` 时，blob 内存储的不是原始文件内容，而是一个 JSON 指针：

```jsonc
// .chemvcs/objects/00/11aabb... （POTCAR 引用 blob）
{
  "type": "reference",
  "original_path": "/opt/vasp/pp/PBE/Li_sv/POTCAR",
  "original_size_bytes": 2451678,
  "sha256": "a1b2c3d4e5f6...",
  "elements": [
    {"symbol": "Li", "functional": "PBE", "titel": "PAW_PBE Li_sv 10Sep2004", "hash": "..."},
    {"symbol": "Co", "functional": "PBE", "titel": "PAW_PBE Co 06Sep2000", "hash": "..."},
    {"symbol": "O",  "functional": "PBE", "titel": "PAW_PBE O 08Apr2002",  "hash": "..."}
  ],
  "reason": "VASP license restricts redistribution"
}
```

```jsonc
// .chemvcs/objects/22/33ccdd... （WAVECAR 引用 blob，通过 add --force）
{
  "type": "reference",
  "original_path": "/scratch/liming/LCO/WAVECAR",
  "original_size_bytes": 12884901888,
  "sha256": "deadbeef...",
  "reason": "large_file_hash_only"
}
```

### 5.3 去重

基于 content-addressable 存储天然去重。写入流程：

```python
def write_blob(content: bytes, objects_dir: Path, tmp_dir: Path) -> str:
    """写入 blob，返回 hash。内容相同时跳过写入。"""
    h = hashlib.sha256(content).hexdigest()
    prefix, suffix = h[:2], h[2:]
    target = objects_dir / prefix / suffix
    
    if target.exists():
        return h  # 已存在，跳过（去重）
    
    # 原子写入：先写 tmp，再 rename
    (objects_dir / prefix).mkdir(exist_ok=True)
    tmp_path = tmp_dir / f"blob_{h}"
    tmp_path.write_bytes(content)
    os.rename(str(tmp_path), str(target))  # 同文件系统 → 原子操作
    return h

def write_blob_compressed(content: bytes, objects_dir: Path, 
                          tmp_dir: Path, threshold: int = 200_000_000) -> str:
    """大文件自动 gzip 压缩后存储"""
    if len(content) >= threshold:
        compressed = gzip.compress(content, compresslevel=6)
        h = hashlib.sha256(content).hexdigest()  # hash 基于原始内容
        # blob 文件名后缀 .gz 标识压缩
        target = objects_dir / h[:2] / (h[2:] + ".gz")
        # ... 同上原子写入 ...
        return h
    return write_blob(content, objects_dir, tmp_dir)
```

### 5.4 清理策略（`chemvcs gc`）

MVP 中不实现自动清理，仅提供手动命令：

```
chemvcs gc --dry-run    # 列出不可达 blob（无任何 commit 引用）
chemvcs gc              # 删除不可达 blob
chemvcs gc --aggressive # 重新压缩所有 ≥200MB 的 blob
```

**可达性判定**：从 HEAD 沿 parent_id 遍历所有 commit，收集所有 `blob_hash`；objects/ 中不在集合内的 blob 视为不可达。

---

## 6. HPC 集成（MVP 范围）

### 6.1 SLURM 适配器

MVP 目标：**生成 sbatch 脚本 + 捕获 job ID + 查询状态**。不做任务编排。

```python
# chemvcs/hpc/slurm.py

class SlurmAdapter:
    """SLURM 集成——仅调用系统命令，无守护进程"""
    
    def generate_script(self, template: dict, vasp_command: str) -> str:
        """
        基于模板和当前项目配置生成 sbatch 脚本。
        
        template 来自 .chemvcs/config.json 中的 slurm 字段：
        {
          "partition": "normal",
          "nodes": 1,
          "ntasks_per_node": 24,
          "time": "24:00:00",
          "account": "mat-group",
          "modules": ["vasp/6.3.2", "intel/2023"]
        }
        """
        script = f"""#!/bin/bash
#SBATCH --job-name=chemvcs_{self.project_name}
#SBATCH --partition={template['partition']}
#SBATCH --nodes={template['nodes']}
#SBATCH --ntasks-per-node={template['ntasks_per_node']}
#SBATCH --time={template['time']}
#SBATCH --account={template.get('account', '')}
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# Auto-generated by ChemVCS {__version__}
# Commit: {self.current_head}

{chr(10).join(f'module load {m}' for m in template.get('modules', []))}

{vasp_command}
"""
        return script
    
    def submit(self, script_path: str) -> str:
        """调用 sbatch，返回 job ID"""
        result = subprocess.run(
            ["sbatch", script_path],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            raise SlurmError(f"sbatch failed: {result.stderr}")
        # "Submitted batch job 12345"
        job_id = result.stdout.strip().split()[-1]
        return job_id
    
    def query_status(self, job_id: str) -> dict:
        """调用 squeue，返回状态"""
        result = subprocess.run(
            ["squeue", "-j", job_id, "--format=%T|%M|%R", "--noheader"],
            capture_output=True, text=True, timeout=10
        )
        if not result.stdout.strip():
            return {"state": "COMPLETED_OR_UNKNOWN", "job_id": job_id}
        parts = result.stdout.strip().split("|")
        return {"state": parts[0], "time": parts[1], "reason": parts[2]}
```

### 6.2 环境快照

在 commit 时自动采集，不依赖 SLURM：

```python
# chemvcs/hpc/env_snapshot.py

def capture_environment() -> dict:
    env = {
        "hostname": socket.gethostname(),
        "python_version": platform.python_version(),
        "chemvcs_version": __version__,
        "vasp_version": None,
        "loaded_modules": None,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    
    # 尝试捕获 module list
    try:
        result = subprocess.run(
            ["bash", "-c", "module list 2>&1"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            env["loaded_modules"] = parse_module_list(result.stdout)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass  # module 命令不可用，跳过
    
    return env
```

### 6.3 设计边界

- **不维护任务队列**：`chemvcs submit` 只是 `sbatch` 的包装
- **不做自动重试**：失败后由用户决定
- **不做跨集群同步**：MVP 中项目只存在于一个文件系统
- **不做计算资源估算**：模板由用户在 config.json 中配置

---

## 7. 合规与安全

### 7.1 POTCAR / VASP 许可

```
决策链：
1. `chemvcs add POTCAR` 
   → 检测 file_type == POTCAR（通过文件名匹配 + 内容头部 "PAW_PBE" 特征）
   → 强制 is_reference = true，代码路径中无条件跳过 blob 内容写入
   → 打印警告：⚠️ POTCAR 受 VASP 许可限制，仅记录哈希引用

2. reprodue 时：
   → 打印赝势路径 + 预期 hash
   → 不尝试从任何源下载或拷贝 POTCAR

3. 未来 push/pull 时：
   → transmit manifest 中 POTCAR_REF 条目仅含 hash，不含内容
   → 接收端需自行提供赝势
```

**代码层强制阻断：**

```python
# chemvcs/core/staging.py

def add_file(path: Path, objects_dir: Path, tmp_dir: Path) -> FileEntry:
    file_type = detect_vasp_file_type(path)
    
    if file_type == "POTCAR_REF":
        # 硬阻断：永远不写 POTCAR 的原始内容到 blob
        return _add_as_reference(path, file_type, reason="vasp_license")
    
    if file_type in ("WAVECAR", "CHGCAR", "CHG", "PROCAR"):
        # 默认忽略；--force 时走引用模式
        raise IgnoredFileError(f"{path.name} 在 .chemvcsignore 中，使用 --force 记录哈希引用")
    
    # 正常文件：全量写入 blob
    content = path.read_bytes()
    blob_hash = write_blob_compressed(content, objects_dir, tmp_dir)
    return FileEntry(path=str(path), blob_hash=blob_hash, 
                     size_bytes=len(content), file_type=file_type, is_reference=False)
```

### 7.2 数据隐私

- ChemVCS 不采集、不上传任何用户数据
- 遥测数据写入本地文件，默认关闭
- commit 中的 `author` 字段来自 `$USER@$HOSTNAME`，用户可通过 `config.json` 覆盖为匿名值

---

## 8. 失败模式与恢复

### 8.1 故障分类与应对

| # | 故障场景 | 检测方式 | 影响范围 | 恢复策略 |
|---|----------|----------|----------|----------|
| F1 | **commit 写入中途断电/kill** | `tmp/` 下有残留文件 | 暂存的 blob 丢失 | 启动时清理 `tmp/`；未写入 `commits/` 的 commit 不存在于历史中——安全（无脏状态） |
| F2 | **blob 写入成功但 commit JSON 未写入** | HEAD 未更新，commit 不在链表中 | 孤立 blob | `chemvcs gc` 可清理；数据不丢失（blob 是幂等的） |
| F3 | **commit JSON 写入成功但 SQLite 未更新** | `metadata.db` 缺少该 commit | 索引落后于真相 | `chemvcs rebuild-index` 或下次启动时自动检测 HEAD 与 DB 不一致并修复 |
| F4 | **SQLite 数据库损坏** | `PRAGMA integrity_check` 失败 | 索引不可用 | 删除 `metadata.db`，运行 `chemvcs rebuild-index` 从 commits/ 重建 |
| F5 | **并发写入（两个终端同时 commit）** | `fcntl.flock` 串行化 | 第二个写入者等待 | flock 超时 5s 后报错 `另一个 chemvcs 进程正在写入，请稍后重试` |
| F6 | **共享文件系统延迟（Lustre 元数据缓存）** | 新写入文件 stat() 短暂不可见 | hash 校验误报 | blob 写入后显式 `os.fsync(fd)` 再 close；校验前 re-stat |
| F7 | **blob 文件损坏（bit rot）** | `chemvcs fsck` 校验 hash | 该 commit 的 reproduce 失败 | 报告损坏文件列表；若有冗余备份可手动恢复 |
| F8 | **磁盘满** | write 抛出 `OSError` | 当前 commit 失败 | 预检：commit 前估算所需空间（暂存文件总大小 ×1.1）；空间不足时提前报错 |

### 8.2 写入事务伪代码

```python
def do_commit(repo: Repo, message: str, staged_files: list[Path]):
    """
    写入顺序保证：
    1. blobs 先写（幂等，可重入）
    2. commit JSON 再写（原子 rename）
    3. HEAD 更新（原子 rename）
    4. SQLite 最后更新（可重建，允许失败）
    """
    tmp_dir = repo.chemvcs_dir / "tmp"
    tmp_dir.mkdir(exist_ok=True)
    
    # Step 1: 写入所有 blob（幂等——已存在则跳过）
    file_entries = []
    for path in staged_files:
        entry = add_file(path, repo.objects_dir, tmp_dir)
        file_entries.append(entry)
    
    # Step 2: 构建 commit 对象
    commit_obj = build_commit(
        parent_id=repo.read_head(),
        files=file_entries,
        message=message,
        semantic=parse_semantics(file_entries),
        output=extract_output_summary(file_entries),
        environment=capture_environment(),
    )
    commit_hash = compute_commit_hash(commit_obj)
    commit_obj["id"] = commit_hash
    
    # Step 3: 原子写入 commit JSON
    commit_path = repo.commits_dir / commit_hash[:2] / f"{commit_hash[2:]}.json"
    commit_path.parent.mkdir(exist_ok=True)
    tmp_commit = tmp_dir / f"commit_{commit_hash}.json"
    tmp_commit.write_text(json.dumps(commit_obj, indent=2, ensure_ascii=False))
    os.replace(str(tmp_commit), str(commit_path))  # 原子
    
    # Step 4: 原子更新 HEAD
    tmp_head = tmp_dir / "HEAD.tmp"
    tmp_head.write_text(commit_hash)
    os.replace(str(tmp_head), str(repo.chemvcs_dir / "HEAD"))  # 原子
    
    # Step 5: 更新 SQLite 索引（允许失败——可重建）
    try:
        with DBWriteLock(repo.chemvcs_dir / "db.lock"):
            repo.db.insert_commit(commit_obj)
    except Exception as e:
        # 索引更新失败不影响数据完整性
        logger.warning(f"索引更新失败（可用 rebuild-index 修复）: {e}")
    
    # Step 6: 清理暂存区
    (repo.staging_dir / "manifest.json").write_text("[]")
    
    # Step 7: 清理 tmp
    for f in tmp_dir.iterdir():
        f.unlink(missing_ok=True)
    
    return commit_hash
```

---

## 9. 技术选型与理由

| 决策 | 选择 | 备选方案 | 选择理由 |
|------|------|----------|----------|
| **语言** | Python ≥3.8 | Rust, Go | 目标用户生态（pymatgen/ASE 均 Python）；无需编译；HPC 上 Python 几乎必装 |
| **CLI 框架** | Typer 0.9.x | Click, argparse | Typer 基于 Click，类型提示友好，自动生成 --help；0.9.x 兼容 Python 3.8 |
| **VASP 解析** | pymatgen ≥2024.1 | ASE, 自写 | pymatgen 是 VASP I/O 的事实标准；Incar/Poscar/Kpoints/Outcar/Vasprun 全覆盖 |
| **元数据存储** | SQLite 3 (stdlib) | TinyDB, JSON flat file | 标准库自带，无须安装；支持 SQL 查询；WAL 模式性能好 |
| **哈希算法** | SHA-256 | SHA-1, BLAKE3 | SHA-256 安全性足够且 Python hashlib 原生支持；Git 已从 SHA-1 迁出，不重蹈覆辙 |
| **压缩** | gzip (stdlib) | zstd, lz4 | 标准库自带；对 XML 类文本压缩比 ~10:1 足够；不引入 C 依赖 |
| **文件锁** | `fcntl.flock` | `portalocker`, `filelock` | POSIX 原生，无第三方依赖；Windows 开发环境使用 `msvcrt.locking` 适配 |
| **结构比较** | pymatgen StructureMatcher | 自实现 | 成熟的晶体结构比较算法，处理周期性边界和对称性 |
| **测试** | pytest + pytest-cov | unittest | 社区标准；参数化测试方便覆盖 VASP 边缘用例 |
| **分发** | PyPI sdist+wheel + conda-forge | snap, Homebrew | 计算材料用户 90%+ 使用 pip 或 conda |

### 依赖边界

```
必须依赖（install_requires）：
  - pymatgen >=2024.1.0    # VASP 文件解析（MVP 核心）
  - typer >=0.9.0,<1.0     # CLI 框架
  - rich >=13.0             # 终端格式化输出（diff 着色、表格）

可选依赖（extras_require["dev"]）：
  - pytest, pytest-cov      # 测试
  - mypy, ruff              # 代码质量

不依赖（明确排除）：
  - 数据库（PostgreSQL, MongoDB, Redis）
  - 消息队列（RabbitMQ, Celery）
  - Web 框架（Flask, FastAPI）
  - C 编译扩展
```

### Python 3.8 兼容策略

```python
# 禁用的语法（3.9+）
dict1 | dict2          # → {**dict1, **dict2}
list[int]              # → typing.List[int]
tuple[str, ...]        # → typing.Tuple[str, ...]
match x:               # → if/elif 链
(x := expr)            # walrus 可用（3.8+）

# pyproject.toml
[project]
requires-python = ">=3.8"
```

---

## 10. 未来扩展点

### 10.1 插件系统（Phase 4）

通过 Python entry_points 注册新的计算代码解析器：

```toml
# 第三方插件 pyproject.toml
[project.entry-points."chemvcs.parsers"]
quantum_espresso = "chemvcs_qe:QEParserPlugin"
cp2k = "chemvcs_cp2k:CP2KParserPlugin"
```

```python
# chemvcs/parsers/registry.py
from importlib.metadata import entry_points

def load_parsers():
    """发现并加载所有已安装的解析器插件"""
    eps = entry_points()
    # Python 3.8/3.9 兼容写法
    group = eps.get("chemvcs.parsers", [])
    for ep in group:
        plugin = ep.load()
        DiffEngine.register(plugin.file_type)(plugin.differ_class)
        FileDetector.register(plugin.file_type)(plugin.detector)
```

**插件接口约定：**

```python
class ParserPlugin(Protocol):
    file_type: str                          # 如 "QE_INPUT"
    
    @staticmethod
    def detect(path: Path) -> bool:
        """是否能处理此文件"""
        ...
    
    @staticmethod
    def parse(content: bytes) -> dict:
        """解析为结构化数据"""
        ...
    
    differ_class: type[BaseDiffer]          # 语义 diff 实现
    detector: type[FileDetector]            # 文件类型检测
```

### 10.2 远程同步（Phase 3）

设计预留：

```
# 传输单元是 commit JSON + 引用的 blob
# 协议：SSH + rsync（无需自建服务器）

chemvcs remote add hpc ssh://user@hpc.example.com:/scratch/project
chemvcs push hpc        # rsync .chemvcs/commits/ + .chemvcs/objects/ → remote
chemvcs pull hpc        # rsync remote → local，合并线性历史
```

- `POTCAR_REF` 类型的 blob 在 push 时只传输 JSON 引用（<1 KB），不传内容
- 冲突检测：若 remote HEAD 与 local parent 分叉 → 拒绝 push，要求先 pull

### 10.3 Web UI（Phase 4）

- 可选的 `chemvcs serve` 命令启动本地 HTTP 服务（仅开发/展示用）
- 读取 `metadata.db`，提供 commit 历史浏览、diff 可视化、能量趋势图
- 前端：单页静态 HTML（内嵌），无需 Node.js 构建
- 不改变存储层设计——Web UI 是 SQLite 的只读消费者

### 10.4 架构扩展性总结

```
可扩展维度           扩展机制                 影响核心代码？
─────────────       ────────────────         ─────────────
新计算代码           entry_points 插件        否
新调度器(PBS/LSF)   HPC Adapter 子类         否（新增文件）
远程同步             Transport 层新增         否（新增模块）
Web UI              只读 DB consumer         否
新存储后端           ObjectStore 接口替换     小（接口已抽象）
```

---

*本文档为工程设计起点。实施前需在 Lustre/GPFS/NFS 三种文件系统上验证 §3.3 和 §8 的假设。*
