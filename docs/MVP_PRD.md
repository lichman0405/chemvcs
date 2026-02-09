# ChemVCS MVP PRD

> v1.0 | 2026-02-09  
> 状态：Draft  
> 范围：仅 MVP（Phase 1），不含远程协作、Web UI、非 VASP 代码支持

---

## A. 目标与非目标

### Goals

| # | 目标 | 衡量方式 |
|---|------|----------|
| G1 | 让 VASP 用户在超算上用 ≤3 条命令完成一次"参数快照" | `chemvcs init` → `chemvcs add .` → `chemvcs commit` 端到端 <10s |
| G2 | 提供 INCAR / POSCAR / KPOINTS 的语义级 diff，替代肉眼比对 | diff 输出按参数键名对齐，标注数值变化量和物理含义标签 |
| G3 | 让用户在任意历史 commit 上一键还原完整输入文件集 | `chemvcs reproduce <rev>` 还原后文件 SHA-256 与原始 commit 一致 |
| G4 | 零外部服务依赖，`pip install chemvcs` 后即可在 HPC 登录节点使用 | 安装不触发任何 TCP 端口监听；运行时仅读写本地文件系统 |
| G5 | 自动采集 VASP 输出摘要（能量、收敛状态、警告），减少手动检查 | commit 后自动解析 OUTCAR/vasprun.xml，摘要写入元数据 |

### Non-goals（MVP 明确不做）

| # | 非目标 | 原因 |
|---|--------|------|
| N1 | push / pull / clone 远程同步 | 需要设计传输协议和认证，推迟到 Phase 3 |
| N2 | 分支（branch）与合并（merge） | 语义合并复杂度高，MVP 仅支持线性历史 |
| N3 | Quantum ESPRESSO / CP2K / ORCA 支持 | 每种代码的文件格式不同，MVP 聚焦 VASP 验证核心假设 |
| N4 | Web UI / Dashboard | 目标用户在 SSH 终端工作，CLI 是最短路径 |
| N5 | 与 atomate2 / AiiDA 的集成导入导出 | 需要依赖其生态，MVP 保持独立 |

---

## B. 目标用户与使用场景

### 楔子市场：Persona A — 博士研究生 + VASP + SLURM

**选定画像：**

- 李明，材料科学博士二年级，课题组 ~15 人
- 技能：基本 Linux/Python，知道 Git 但不熟练
- 工具：VASP + VASPKIT + pymatgen
- 环境：学校超算集群（SLURM），月 ~200 核时
- 日常：SSH → 调参 → sbatch → 检查 → 循环

**为什么选 Persona A、不选 B/C/D：**

| Persona | 不选原因 |
|---------|----------|
| **B（博后/高通量）** | 需要 branch、批量操作、多超算同步——这些是 Phase 3+ 特性。如果为 B 做 MVP 会膨胀范围。B 是线性历史验证通过后的自然扩展用户。 |
| **C（PI/课题组长）** | C 本人不直接操作计算工具；其价值来自"学生已经在用"。先让 A 用起来，C 的需求（audit、archive）自动成立。 |
| **D（超算运维）** | D 是分发渠道而非终端用户。当 A 类用户量到达阈值后，D 才有动力将 ChemVCS 纳入 module 体系。 |

**核心使用场景（一句话）：** 博士生在超算上对一个 VASP 项目反复调参，需要记录每次修改、比对差异、在写论文时精确回溯某一版参数。

---

## C. 核心用户旅程

```
场景：李明在超算上做 LiCoO₂ 掺杂计算，需要调 Hubbard U 值

┌─────────────────────────────────────────────────────────────┐
│ 1. 初始化项目                                                │
│    $ cd /scratch/liming/LCO_doping                          │
│    $ chemvcs init                                           │
│    → 创建 .chemvcs/ 目录，检测到 INCAR/POSCAR/KPOINTS       │
│                                                             │
│ 2. 首次快照                                                  │
│    $ chemvcs add .                                          │
│    $ chemvcs commit -m "baseline: PBE+U, U=3.5, 4x4x4"     │
│    → 记录所有输入文件 hash + VASP 语义摘要 + 赝势 hash       │
│                                                             │
│ 3. 提交计算 & 等待结果                                        │
│    $ sbatch run.slurm                                       │
│    （VASP 计算完成，生成 OUTCAR / vasprun.xml）               │
│                                                             │
│ 4. 记录结果                                                  │
│    $ chemvcs add OUTCAR vasprun.xml                         │
│    $ chemvcs commit -m "U=3.5 converged, E=-45.23eV"        │
│    → 自动解析 OUTCAR：总能量、是否收敛、离子步数、警告        │
│                                                             │
│ 5. 调参 & 新快照                                             │
│    $ vim INCAR   # 修改 LDAUU = 4.0                         │
│    $ chemvcs add INCAR                                      │
│    $ chemvcs commit -m "U=4.0 试跑"                         │
│                                                             │
│ 6. 比对差异                                                  │
│    $ chemvcs diff HEAD~2 HEAD                               │
│    → 语义输出：                                              │
│      INCAR:  LDAUU  3.5 → 4.0  (+0.5)                      │
│      OUTCAR: E_total  -45.23eV → -44.87eV  (+0.36eV)       │
│                                                             │
│ 7. 查看历史                                                  │
│    $ chemvcs log                                            │
│    → 按时间倒序显示所有 commit，含语义摘要                    │
│                                                             │
│ 8. 三个月后写论文，复现结果                                    │
│    $ chemvcs log --grep "U=3.5"                             │
│    $ chemvcs reproduce abc1234                              │
│    → 将 commit abc1234 的所有输入文件还原到 ./reproduce/      │
│    → 打印赝势路径 & 版本校验结果                              │
│    → 打印当时记录的 VASP 版本号                               │
└─────────────────────────────────────────────────────────────┘
```

---

## D. 用户故事与验收标准

### US-01：初始化项目（`chemvcs init`）

**故事：** 作为博士生，我想在一个 VASP 计算目录中初始化 ChemVCS，以便开始追踪参数变更。

| 项目 | 内容 |
|------|------|
| **前置条件** | 当前目录存在至少一个 VASP 输入文件（INCAR/POSCAR/POTCAR/KPOINTS 任一）或为空目录 |
| **主流程** | 1. 用户执行 `chemvcs init`<br>2. 系统在当前目录创建 `.chemvcs/` 子目录<br>3. 初始化 SQLite 数据库 `.chemvcs/metadata.db`<br>4. 创建默认 `.chemvcsignore`（含 WAVECAR、CHGCAR、CHG、PROCAR 等默认排除项）<br>5. 扫描当前目录，输出检测到的 VASP 文件列表<br>6. 打印初始化成功信息 |
| **异常流程** | E1：当前目录已有 `.chemvcs/` → 报错 `已初始化，无需重复操作`，退出码 1<br>E2：无写权限 → 报错 `权限不足`，退出码 2<br>E3：磁盘空间不足 → 报错并提示清理 |
| **验收标准** | AC1：执行后 `.chemvcs/metadata.db` 存在且可被 `sqlite3` 打开<br>AC2：`.chemvcsignore` 包含 WAVECAR、CHGCAR、CHG、PROCAR、*.tmp<br>AC3：重复执行返回退出码 1 且不破坏已有数据<br>AC4：整个过程 <2s（SSD）/ <5s（Lustre） |

---

### US-02：快照提交（`chemvcs add` + `chemvcs commit`）

**故事：** 作为博士生，我想把当前的输入文件和可选的输出文件记录为一个不可变快照，附带语义元数据。

| 项目 | 内容 |
|------|------|
| **前置条件** | 项目已 `chemvcs init`；至少有一个文件被修改或新增 |
| **主流程** | 1. `chemvcs add <path>...` 将指定文件加入暂存区（stage）<br>2. 对每个暂存文件计算 SHA-256 hash<br>3. 若文件匹配已知 VASP 类型，解析语义（见 E 节模型）<br>4. `chemvcs commit -m "<msg>"` 创建一条 commit 记录<br>5. 将文件内容以 content-addressable blob 存入 `.chemvcs/objects/`<br>6. 写入 commit 元数据到 SQLite（见 E 节字段）<br>7. 若暂存区含 OUTCAR/vasprun.xml，自动提取输出摘要<br>8. 若检测到 POTCAR 路径，记录赝势哈希指纹（不拷贝内容）<br>9. 打印 commit hash（前 7 位）+ 语义摘要 |
| **异常流程** | E1：暂存区为空时 commit → 报错 `无文件可提交`<br>E2：文件在 `.chemvcsignore` 中 → `add` 时警告并跳过<br>E3：文件 >100 MB 且不在 ignore 中 → 警告 `大文件建议加入 .chemvcsignore`，但仍允许 add<br>E4：POTCAR 被 add → 警告 `POTCAR 含许可内容，仅记录哈希引用` |
| **验收标准** | AC1：commit 后 `.chemvcs/objects/` 下新增 blob 文件，SHA-256 与源文件一致<br>AC2：`metadata.db` 的 commits 表新增一行，所有必需字段非空<br>AC3：同一文件内容多次 commit 不重复存储（去重）<br>AC4：POTCAR 仅写入哈希 + 路径，blob 体积 <1 KB<br>AC5：含 OUTCAR 的 commit 自动填充 `output_summary` 字段 |

---

### US-03：查看历史（`chemvcs log`）

**故事：** 作为博士生，我想查看项目的全部提交历史，快速定位某次参数修改。

| 项目 | 内容 |
|------|------|
| **前置条件** | 项目已初始化且至少有 1 个 commit |
| **主流程** | 1. `chemvcs log` 按时间倒序输出所有 commit<br>2. 每条显示：短 hash、时间戳、commit message、INCAR 关键参数摘要（ENCUT/ISMEAR/LDAUU 等）、输出能量（若有）<br>3. 支持 `--oneline` 精简模式<br>4. 支持 `--grep <keyword>` 过滤 commit message 或参数值 |
| **异常流程** | E1：无 commit → 输出 `暂无提交记录`<br>E2：`--grep` 无匹配 → 输出 `无匹配记录` |
| **验收标准** | AC1：输出顺序为时间倒序<br>AC2：`--oneline` 每条 commit 输出严格 1 行：`<hash> <timestamp> <message>`<br>AC3：`--grep "U=3.5"` 仅返回 message 或 INCAR 摘要含该关键词的 commit<br>AC4：100 条 commit 的 log 输出 <1s |

---

### US-04：语义差异比较（`chemvcs diff`）

**故事：** 作为博士生，我想比较两个版本的参数差异，看到的是"ENCUT 从 400 变为 520"而不是文本行的增删。

| 项目 | 内容 |
|------|------|
| **前置条件** | 项目有 ≥2 个 commit，或指定 1 个 commit 与工作区比较 |
| **主流程** | 1. `chemvcs diff <rev1> <rev2>` 比较两个 commit<br>2. `chemvcs diff <rev>` 比较该 commit 与当前工作区<br>3. `chemvcs diff` （无参数）比较最近 commit 与当前工作区<br>4. 对每种文件类型输出语义 diff：<br>&emsp;• **INCAR**：逐 key 比较，标注变更/新增/删除的参数及其值<br>&emsp;• **POSCAR**：晶格常数变化、原子坐标 RMSD、原子数变化、对称性变化<br>&emsp;• **KPOINTS**：网格密度变化、方案类型变化（Gamma→MP 等）<br>&emsp;• **OUTCAR/vasprun.xml**：总能量差、收敛状态对比、离子步数对比<br>&emsp;• **其他文件**：退回逐行 diff<br>5. 数值变化显示差值（如 `ENCUT: 400 → 520 (+120 eV)`） |
| **异常流程** | E1：指定的 rev 不存在 → 报错 `未找到版本 <rev>`<br>E2：两个 rev 完全相同 → 输出 `两个版本无差异`<br>E3：某文件只在一侧存在 → 标注 `[新增]` 或 `[删除]` |
| **验收标准** | AC1：INCAR diff 输出格式为 `KEY: old_value → new_value (±delta unit)`（数值型参数含差值）<br>AC2：POSCAR diff 输出晶格常数变化精度到 0.001 Å，坐标 RMSD 精度到 0.0001 Å<br>AC3：KPOINTS diff 输出网格尺寸对比（如 `4×4×4 → 6×6×6`）<br>AC4：无已知语义解析器的文件退回标准 unified diff<br>AC5：diff 输出支持 `--format json` 生成机器可读 JSON |

---

### US-05：复现历史计算（`chemvcs reproduce`）

**故事：** 作为博士生，我在写论文时需要精确复现三个月前的某次计算，还原所有输入文件并验证环境一致性。

| 项目 | 内容 |
|------|------|
| **前置条件** | 目标 commit 存在且 blob 完整（无损坏） |
| **主流程** | 1. `chemvcs reproduce <rev>` 在当前目录下创建 `reproduce_<rev>/` 子目录<br>2. 从 `.chemvcs/objects/` 还原所有 input 文件到该目录<br>3. 校验每个文件的 SHA-256 与 commit 记录一致<br>4. 打印赝势引用信息：原始路径 + 预期哈希<br>5. 检查当前赝势路径是否存在且哈希匹配，输出 ✅ 或 ⚠️<br>6. 打印当时记录的环境信息（VASP 版本、编译器等）<br>7. 若当前环境信息可检测，与记录对比并标注差异<br>8. 支持 `--output-dir <path>` 自定义还原目录 |
| **异常流程** | E1：blob 文件损坏（hash 不匹配）→ 报错 `数据完整性校验失败：<filename>`，列出损坏文件，退出码 3<br>E2：赝势路径不存在 → ⚠️ 警告但不阻塞还原<br>E3：目标目录已存在 → 报错 `目录已存在，请删除或使用 --output-dir` |
| **验收标准** | AC1：还原目录中的每个文件的 SHA-256 与 commit 元数据中记录的 hash 完全一致<br>AC2：赝势校验结果区分 ✅ 匹配 / ⚠️ 路径不存在 / ❌ 哈希不匹配<br>AC3：还原操作不修改 `.chemvcs/` 中的任何数据<br>AC4：`--output-dir /tmp/test` 正确写入指定位置<br>AC5：还原 5 个输入文件 + 赝势校验全流程 <3s |

---

## E. 数据与元数据模型

### E1. 存储结构

```
project_root/
├── INCAR, POSCAR, KPOINTS, ...   # 用户工作区
├── .chemvcs/
│   ├── metadata.db               # SQLite，commit 历史 + 索引
│   ├── objects/                   # content-addressable blob 存储
│   │   ├── ab/cd1234...          # SHA-256 前2位为目录，后62位为文件名
│   │   └── ...
│   ├── staging/                  # 暂存区（add 后 commit 前）
│   │   └── manifest.json         # 暂存文件清单
│   └── config.json               # 项目级配置
└── .chemvcsignore                # 排除规则
```

### E2. Commit 记录（`commits` 表）

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `id` | TEXT PK | ✅ | SHA-256 of (parent_id + tree_hash + timestamp + message) |
| `parent_id` | TEXT | ✅ | 父 commit id；首次 commit 为 `NULL` |
| `tree_hash` | TEXT | ✅ | 本次 commit 所有文件 hash 的 Merkle root |
| `message` | TEXT | ✅ | 用户输入的 commit message |
| `timestamp` | TEXT | ✅ | ISO 8601 UTC 时间 |
| `author` | TEXT | ✅ | `$USER@$HOSTNAME` |

### E3. 文件快照（`files` 表）

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `commit_id` | TEXT FK | ✅ | 所属 commit |
| `path` | TEXT | ✅ | 相对于项目根目录的文件路径 |
| `blob_hash` | TEXT | ✅ | 文件内容 SHA-256 |
| `size_bytes` | INTEGER | ✅ | 文件大小 |
| `file_type` | TEXT | ✅ | 枚举：`INCAR` / `POSCAR` / `KPOINTS` / `POTCAR_REF` / `OUTCAR` / `VASPRUN` / `OTHER` |
| `is_reference` | BOOLEAN | ✅ | 是否仅为引用（如 POTCAR），`true` 时 blob 仅含元信息 |

### E4. 语义摘要（`semantic_summary` 表）

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `commit_id` | TEXT FK | ✅ | 所属 commit |
| `incar_params` | JSON | ❌ | INCAR 关键参数快照，如 `{"ENCUT": 520, "ISMEAR": 0, "SIGMA": 0.05, "LDAUU": "3.5 0 0"}` |
| `poscar_formula` | TEXT | ❌ | 化学式，如 `Li4Co3FeO8` |
| `poscar_spacegroup` | TEXT | ❌ | 空间群符号，如 `R-3m` |
| `poscar_lattice` | JSON | ❌ | 晶格常数 `{"a": 2.83, "b": 2.83, "c": 14.05, "alpha": 90, "beta": 90, "gamma": 120}` |
| `poscar_natoms` | INTEGER | ❌ | 总原子数 |
| `kpoints_grid` | TEXT | ❌ | 如 `4x4x4 Gamma` |
| `potcar_refs` | JSON | ❌ | `[{"element": "Li", "hash": "a1b2...", "path": "/opt/vasp/pp/Li_sv/POTCAR"}]` |

### E5. 输出摘要（`output_summary` 表）

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `commit_id` | TEXT FK | ✅ | 所属 commit |
| `total_energy_eV` | REAL | ❌ | 最终总能量 (eV) |
| `energy_per_atom_eV` | REAL | ❌ | 每原子能量 (eV) |
| `is_converged` | BOOLEAN | ❌ | 电子步是否收敛 |
| `ionic_steps` | INTEGER | ❌ | 离子弛豫步数 |
| `max_force_eV_A` | REAL | ❌ | 最大残余力 (eV/Å) |
| `bandgap_eV` | REAL | ❌ | 带隙 (eV)，若输出中可提取 |
| `magnetic_moment` | REAL | ❌ | 总磁矩 (μB) |
| `warnings` | JSON | ❌ | 警告列表，如 `["BRMIX: very serious problems"]` |
| `vasp_version` | TEXT | ❌ | 如 `6.3.2` |
| `elapsed_time_s` | REAL | ❌ | VASP 计算总耗时（秒） |

### E6. 环境快照（`environment` 表）

| 字段 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `commit_id` | TEXT FK | ✅ | 所属 commit |
| `hostname` | TEXT | ✅ | 计算节点主机名 |
| `vasp_version` | TEXT | ❌ | 从 OUTCAR 头部提取 |
| `compiler_info` | TEXT | ❌ | 编译器信息（若 OUTCAR 中有） |
| `loaded_modules` | JSON | ❌ | `module list` 快照，如 `["vasp/6.3.2", "intel/2023"]` |
| `python_version` | TEXT | ✅ | ChemVCS 运行时 Python 版本 |
| `chemvcs_version` | TEXT | ✅ | ChemVCS 自身版本号 |
| `slurm_job_id` | TEXT | ❌ | `$SLURM_JOB_ID`（若在 SLURM 环境中） |

---

## F. 约束与边界

### F1. 大文件不入库策略

| 文件 | 典型大小 | 策略 |
|------|----------|------|
| WAVECAR | 1–50 GB | **默认 `.chemvcsignore`**，不追踪不存储 |
| CHGCAR | 100 MB–2 GB | 默认忽略；用户可显式 `add --force` 纳入 |
| CHG | 50–500 MB | 默认忽略 |
| PROCAR | 100 MB–1 GB | 默认忽略 |
| vasprun.xml | 10–500 MB | **默认追踪**，但仅提取摘要字段存入 `output_summary`；原始文件以 blob 存入 objects（≥200 MB 时 gzip 压缩） |

- 被忽略的文件在 `chemvcs status` 中显示为 `[ignored]`
- `chemvcs add --force WAVECAR` 可覆盖忽略规则，此时记录 hash 指纹 + 文件大小但 **不** 拷贝 blob；元数据标记 `is_reference = true`

### F2. 赝势 / 许可文件策略

- POTCAR 受 VASP 许可证保护，ChemVCS **永远不复制其内容**
- `chemvcs add POTCAR` → 自动转为引用模式：
  - 计算 SHA-256 hash
  - 记录绝对路径
  - 按元素拆分记录（`Li_sv`, `Co`, `O` 各自 hash）
  - blob 文件仅写入 JSON 元信息（<1 KB）
- `chemvcs reproduce` 时打印赝势路径与预期哈希，由用户确认一致性
- 代码中不硬编码任何赝势内容

### F3. HPC 环境约束

| 约束 | 设计应对 |
|------|----------|
| 无 root 权限 | 纯 `pip install --user`；不使用系统级路径 |
| 无 TCP 服务 | SQLite 本地文件；无 daemon 进程 |
| 共享文件系统（Lustre/GPFS/NFS） | SQLite 使用 `journal_mode=WAL`；检测 NFS 时回退 `journal_mode=DELETE`；写操作加 `fcntl.flock` |
| 登录节点资源受限 | CLI 命令峰值内存 <200 MB；commit 操作在 O(文件数×文件大小) 内完成，不做全局索引重建 |
| Python 版本多样 | 最低支持 Python 3.8（CentOS 7 常见版本）；无 C 扩展 |
| 网络访问受限 | 安装可通过离线 wheel；运行时 100% 离线 |

---

## G. 指标与埋点

所有遥测默认 **关闭**，用户通过 `chemvcs config telemetry on` 自愿启用。启用后数据写入本地 `.chemvcs/telemetry.jsonl`，不发送到网络。

| # | 指标名 | 采集方式 | 目标值（T+4 周 Alpha） |
|---|--------|----------|------------------------|
| M1 | `install_to_first_commit_seconds` | 首次 `commit` 的时间戳 − `init` 的时间戳 | ≤300s (5 min) |
| M2 | `commits_per_project_per_week` | 按项目统计周 commit 数 | ≥3（说明用户持续使用） |
| M3 | `diff_invocations_per_week` | `diff` 命令调用次数 | ≥1（验证语义 diff 价值） |
| M4 | `reproduce_success_rate` | reproduce 成功次数 / 总调用次数 | ≥95% |
| M5 | `reproduce_hash_mismatch_count` | reproduce 时文件 hash 不匹配的次数 | 0（完整性目标） |
| M6 | `avg_commit_duration_ms` | `commit` 命令的平均耗时 | ≤5000ms |
| M7 | `large_file_warning_count` | 大文件警告触发次数 | 用于优化默认 ignore 规则 |
| M8 | `error_rate_by_command` | 各命令的非零退出码比例 | ≤5% |

---

## H. 风险与开放问题

| # | 风险/问题 | 类型 | 影响 | 当前应对 | 状态 |
|---|----------|------|------|----------|------|
| R1 | **SQLite 在 Lustre/GPFS 上的并发锁** — WAL 模式在部分共享文件系统上不被支持 | 技术 | 高 | 启动时运行 `PRAGMA journal_mode=WAL` 并捕获失败，回退到 `DELETE` 模式；写操作外层加 `fcntl.flock` | 需在 3 种文件系统（Lustre/GPFS/NFS）上实测 |
| R2 | **content-addressable 存储在共享文件系统上的 rename 原子性** — `os.rename()` 跨文件系统不原子 | 技术 | 高 | blob 先写入 `.chemvcs/tmp/`（同文件系统），再 `os.rename` 到 `objects/`。写入失败时清理临时文件 | 已有设计方案 |
| R3 | **POTCAR 合规风险** — 即使仅存 hash，若用户误操作仍可能泄露 POTCAR 内容 | 合规 | 中 | `add POTCAR` 时强制走引用路径，代码层面阻断 blob 写入；README 和 `--help` 中加入许可声明 | 已有设计方案 |
| R4 | **pymatgen INCAR 解析边界** — 用户可能使用非标准 INCAR 写法（行内注释、多行值） | 技术 | 中 | 先用 pymatgen 解析，失败后退回 key=value 正则；解析失败的参数标记 `[unparsed]` 而非崩溃 | 需收集边缘用例 |
| R5 | **用户不愿多输入一条命令** — `add` + `commit` 两步流程对非 Git 用户有阻力 | 产品 | 高 | 提供 `chemvcs commit -a -m "msg"`（auto-add 所有已追踪文件）；未来可加 sbatch 钩子自动 commit | MVP 中实现 `-a` 快捷方式 |
| R6 | **OUTCAR 解析失败率** — 非正常终止的 VASP 任务产生截断的 OUTCAR | 技术 | 中 | 解析器使用 try-except 包裹每个字段提取；截断文件标记 `is_converged=false, warnings=["truncated output"]` | 需用 10+ 个真实 OUTCAR 样本测试 |
| R7 | **存储膨胀** — 用户频繁 commit 大型 vasprun.xml（200–500 MB）导致 `.chemvcs/objects/` 增长快 | 技术 | 中 | 相同内容 hash 去重；≥200 MB 的 blob 自动 gzip（压缩比约 10:1 for XML）；提供 `chemvcs gc` 清理不可达 blob | MVP 中实现去重 + gzip；gc 延后 |
| R8 | **POSCAR 格式多样性** — VASP 4/5 格式差异（有无元素行）、Cartesian vs. Direct 坐标 | 技术 | 低 | 依赖 pymatgen `Structure.from_file()` 统一处理；diff 基于解析后的 Structure 对象而非原始文本 | pymatgen 已验证支持 |
| R9 | **Python 3.8 兼容性** — 老旧 HPC 系统仍运行 CentOS 7 + Python 3.8 | 技术 | 中 | CI 矩阵包含 3.8/3.9/3.10/3.11/3.12；避免使用 3.9+ 语法（`dict \| dict`、`match-case`）；Typer 需锁定兼容版本 | 需在 CI 中验证 |
| R10 | **"为什么不直接用 Git + 脚本"** — 用户质疑 ChemVCS 相比 Git alias 的增量价值 | 产品 | 高 | 核心差异化：语义 diff（Git 做不到 `LDAUU: 3.5→4.0`）、OUTCAR 自动摘要、赝势哈希校验、reproduce 端到端验证。需在 README 和首次 `init` 时展示对比示例 | 需准备 "ChemVCS vs Git" 对比文档 |

---

*本 PRD 仅覆盖 MVP（T+4 周 Alpha 交付范围）。Branch、远程同步、多代码支持等特性见长期路线图。*
