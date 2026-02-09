# ChemVCS CLI 接口规格

> v1.0 | 2026-02-09  
> 状态：Draft  
> 范围：MVP 命令集（Phase 1 核心版本控制）

---

## 通用约定

### 退出码

| 退出码 | 含义 | 何时使用 |
|--------|------|----------|
| `0` | 成功 | 命令正常完成 |
| `1` | 用户错误 | 参数不合法、项目未初始化、暂存区为空等用户可自行修复的错误 |
| `2` | 系统错误 | 权限不足、磁盘满、文件系统错误等系统级问题 |
| `3` | 数据完整性错误 | blob hash 不匹配、元数据损坏等严重错误 |
| `130` | 用户中断 | Ctrl+C 或 SIGINT |

### 输出格式

- **默认输出**：人类可读的格式化文本（使用 `rich` 库）
- **JSON 输出**：所有命令支持 `--format json` 输出机器可读结果
- **静默模式**：`--quiet` / `-q` 仅输出错误信息
- **详细模式**：`--verbose` / `-v` 输出调试信息

### 颜色方案

```
成功/肯定     → 绿色（rich.green）
警告         → 黄色（rich.yellow）
错误/删除     → 红色（rich.red）
新增         → 青色（rich.cyan）
修改         → 蓝色（rich.blue）
元数据/哈希   → 灰色（rich.dim）
高亮值       → 粗体（rich.bold）
```

### 时间格式

- **人类可读**：`2026-02-09 14:32:01`（本地时区）
- **ISO 8601**：`2026-02-09T14:32:01.123456Z`（存储格式，UTC）
- **相对时间**：`3 days ago`（log 输出中）

---

## 命令详细规格

### `chemvcs init`

**用途**：在当前目录初始化 ChemVCS 项目

**签名**：
```bash
chemvcs init [OPTIONS]
```

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--force` | flag | false | 覆盖已存在的 `.chemvcs/` 目录（危险操作，需二次确认） |
| `--quiet` | flag | false | 静默模式，仅输出错误 |

**行为**：

1. 检查当前目录是否已初始化（`.chemvcs/` 存在）→ 若存在且无 `--force` → 退出码 1
2. 创建 `.chemvcs/` 目录结构（见架构文档 §3.1）
3. 初始化 SQLite 数据库（`metadata.db`）
4. 创建默认 `.chemvcsignore`
5. 扫描当前目录，检测 VASP 文件
6. 输出初始化成功信息 + 检测到的文件列表

**输出示例**：

```
✓ Initialized ChemVCS repository in /scratch/liming/LCO_doping/.chemvcs

Detected VASP files:
  INCAR       (842 bytes)
  POSCAR      (1204 bytes)
  KPOINTS     (56 bytes)
  POTCAR      (2.3 MB, will be tracked as reference)

Next steps:
  chemvcs add .
  chemvcs commit -m "Initial snapshot"

Tip: WAVECAR, CHGCAR are ignored by default (see .chemvcsignore)
```

**JSON 输出** (`--format json`)：

```json
{
  "success": true,
  "repo_path": "/scratch/liming/LCO_doping",
  "chemvcs_dir": "/scratch/liming/LCO_doping/.chemvcs",
  "detected_files": [
    {"path": "INCAR", "size": 842, "type": "INCAR"},
    {"path": "POSCAR", "size": 1204, "type": "POSCAR"},
    {"path": "KPOINTS", "size": 56, "type": "KPOINTS"},
    {"path": "POTCAR", "size": 2451678, "type": "POTCAR_REF"}
  ]
}
```

**错误示例**：

```
✗ Error: Already initialized (.chemvcs/ exists)
  Use --force to reinitialize (WARNING: will delete existing history)
Exit code: 1
```

---

### `chemvcs add`

**用途**：将文件添加到暂存区

**签名**：
```bash
chemvcs add [OPTIONS] <PATH>...
```

**参数**：

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `PATH` | path | ✅ | 文件或目录路径（支持多个）；`.` 表示当前目录所有文件 |

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--force` | flag | false | 覆盖 `.chemvcsignore` 规则（对 WAVECAR 等大文件仍记录为引用） |
| `--dry-run` | flag | false | 仅显示将被添加的文件，不实际修改暂存区 |

**行为**：

1. 检查项目是否已初始化 → 未初始化 → 退出码 1
2. 解析路径（支持通配符，如 `*.vasp`）
3. 对每个文件：
   - 检查是否在 `.chemvcsignore` 中 → 跳过（除非 `--force`）
   - 检测文件类型（INCAR / POSCAR / POTCAR / OTHER）
   - 若 `POTCAR_REF` → 自动转为引用模式 + 打印警告
   - 若大文件（>100 MB）→ 打印警告
   - 计算 SHA-256 hash（不立即写入 blob）
   - 更新 `.chemvcs/staging/manifest.json`
4. 输出已暂存文件列表

**输出示例**：

```
Staging files...
  ✓ INCAR       (842 bytes, SHA: fc1a2b3c...)
  ✓ POSCAR      (1204 bytes, SHA: ab12cd34...)
  ⚠ POTCAR      (2.3 MB, tracked as reference only)
  ⊘ WAVECAR     (ignored, use --force to track as reference)

3 files staged for commit
```

**警告示例**：

```
⚠ Warning: OUTCAR is 345 MB
  Consider adding it to .chemvcsignore if it's reproducible
  Use --force to proceed
```

**JSON 输出**：

```json
{
  "staged": [
    {"path": "INCAR", "size": 842, "hash": "fc1a2b3c...", "type": "INCAR", "is_reference": false},
    {"path": "POSCAR", "size": 1204, "hash": "ab12cd34...", "type": "POSCAR", "is_reference": false},
    {"path": "POTCAR", "size": 2451678, "hash": "00aabbcc...", "type": "POTCAR_REF", "is_reference": true}
  ],
  "ignored": [
    {"path": "WAVECAR", "reason": ".chemvcsignore"}
  ]
}
```

---

### `chemvcs commit`

**用途**：提交暂存区的文件为一个快照

**签名**：
```bash
chemvcs commit [OPTIONS]
```

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `-m, --message` | string | (interactive) | commit 消息（必需） |
| `-a, --all` | flag | false | 自动 add 所有已追踪过的文件（不含新文件） |
| `--allow-empty` | flag | false | 允许空提交（用于里程碑标记） |
| `--author` | string | `$USER@$HOST` | 覆盖默认作者信息 |

**行为**：

1. 检查项目是否已初始化
2. 若 `-a`：自动对所有已追踪文件执行 `add`
3. 检查暂存区 → 为空且无 `--allow-empty` → 退出码 1
4. 若未提供 `-m`：启动编辑器（`$EDITOR` 或 `vim`）编辑消息
5. 对每个暂存文件：
   - 写入 blob 到 `.chemvcs/objects/`（去重）
   - 若文件类型已知 → 解析语义摘要
   - 若 OUTCAR/vasprun.xml → 提取输出摘要
6. 构建 commit 对象（JSON）
7. 计算 commit hash
8. 原子写入 commit JSON 到 `.chemvcs/commits/`
9. 更新 HEAD
10. 更新 SQLite 索引
11. 清空暂存区
12. 输出 commit hash + 摘要

**输出示例**：

```
Committing 3 files...
  [1/3] INCAR       → blob fc1a2b (ENCUT=520, ISMEAR=0, LDAUU=3.5 0 0)
  [2/3] POSCAR      → blob ab12cd (Li4Co4O8, R-3m, 16 atoms)
  [3/3] POTCAR      → ref 00aabb (Li_sv, Co, O)

✓ Committed abc1234
  Author: liming@hpc-login01
  Date:   2026-02-09 14:32:01
  
  PBE+U, U=3.5eV, kpoints 4x4x4
```

**JSON 输出**：

```json
{
  "commit_id": "abc1234567890abcdef1234567890abcdef1234567890abcdef1234567890ab",
  "short_id": "abc1234",
  "parent_id": "dead000...",
  "timestamp": "2026-02-09T14:32:01.123456Z",
  "author": "liming@hpc-login01",
  "message": "PBE+U, U=3.5eV, kpoints 4x4x4",
  "files": 3,
  "semantic_summary": {
    "incar": {"ENCUT": 520, "ISMEAR": 0, "LDAUU": [3.5, 0, 0]},
    "poscar": {"formula": "Li4Co4O8", "spacegroup": "R-3m", "natoms": 16}
  }
}
```

**错误示例**：

```
✗ Error: Nothing staged for commit
  Use 'chemvcs add <file>' to stage files
  Or 'chemvcs commit -a' to commit all tracked files
Exit code: 1
```

---

### `chemvcs log`

**用途**：查看提交历史

**签名**：
```bash
chemvcs log [OPTIONS] [REVISION_RANGE]
```

**参数**：

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `REVISION_RANGE` | string | ❌ | 如 `HEAD~5..HEAD` 或单个 revision（默认：从 HEAD 回溯所有） |

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--oneline` | flag | false | 精简模式，每条 commit 一行 |
| `-n, --max-count` | int | ∞ | 最多显示 N 条 |
| `--grep` | string | - | 过滤 commit message 或语义摘要（不区分大小写） |
| `--since` | string | - | 仅显示指定日期之后的 commit（如 `2026-01-01`） |
| `--until` | string | - | 仅显示指定日期之前的 commit |
| `--author` | string | - | 过滤作者 |
| `--format` | enum | `default` | 输出格式：`default` / `oneline` / `json` |

**行为**：

1. 检查项目是否已初始化
2. 从 HEAD 沿 `parent_id` 回溯
3. 应用过滤条件（--grep / --since / --author）
4. 对每条 commit：
   - 显示短 hash（7 位）
   - 显示时间戳（相对时间 + 绝对时间）
   - 显示作者
   - 显示 message
   - 显示语义摘要关键参数
   - 若有输出摘要 → 显示能量/收敛状态
5. 按时间倒序输出

**输出示例（默认）**：

```
commit abc1234 (HEAD)
Author: liming@hpc-login01
Date:   3 hours ago (2026-02-09 14:32:01)

    PBE+U, U=4.0eV, kpoints 6x6x6

    INCAR:  ENCUT=520, ISMEAR=0, LDAUU=4.0 0 0
    POSCAR: Li4Co3FeO8, R-3m, 16 atoms
    Output: E=-44.87 eV, converged ✓

commit dead000
Author: liming@hpc-login01
Date:   1 day ago (2026-02-08 10:15:42)

    PBE+U, U=3.5eV, kpoints 4x4x4

    INCAR:  ENCUT=520, ISMEAR=0, LDAUU=3.5 0 0
    POSCAR: Li4Co4O8, R-3m, 16 atoms
    Output: E=-45.23 eV, converged ✓
```

**输出示例（--oneline）**：

```
abc1234  2026-02-09 14:32  PBE+U, U=4.0eV, kpoints 6x6x6
dead000  2026-02-08 10:15  PBE+U, U=3.5eV, kpoints 4x4x4
```

**JSON 输出**：

```json
{
  "commits": [
    {
      "id": "abc1234...",
      "short_id": "abc1234",
      "parent_id": "dead000...",
      "author": "liming@hpc-login01",
      "timestamp": "2026-02-09T14:32:01Z",
      "timestamp_relative": "3 hours ago",
      "message": "PBE+U, U=4.0eV, kpoints 6x6x6",
      "semantic_summary": {...},
      "output_summary": {"total_energy_eV": -44.87, "is_converged": true}
    }
  ],
  "total": 2
}
```

---

### `chemvcs diff`

**用途**：比较两个版本或工作区的语义差异

**签名**：
```bash
chemvcs diff [OPTIONS] [REV1] [REV2]
```

**参数**：

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `REV1` | string | ❌ | 源版本（默认：HEAD） |
| `REV2` | string | ❌ | 目标版本（默认：工作区） |

**行为**：

- `chemvcs diff` → HEAD vs 工作区
- `chemvcs diff HEAD~1` → HEAD~1 vs 工作区
- `chemvcs diff v1 v3` → commit v1 vs commit v3

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--files` | string | - | 仅比较指定文件（逗号分隔，如 `INCAR,POSCAR`） |
| `--format` | enum | `semantic` | 输出格式：`semantic` / `unified` / `json` |
| `--no-color` | flag | false | 禁用颜色输出 |
| `--summary` | flag | false | 仅显示变更摘要（不显示详细 diff） |

**行为**：

1. 检查项目是否已初始化
2. 解析 REV1 和 REV2（支持 `HEAD`, `HEAD~N`, commit hash, tag）
3. 获取两个版本的文件列表
4. 对每个文件：
   - 检测文件类型
   - 调用对应的语义 differ（INCAR / POSCAR / KPOINTS / OUTCAR）
   - 若无语义解析器 → 退回 unified diff
5. 输出差异

**输出示例（INCAR 语义 diff）**：

```
── INCAR diff ──────────────────────────────────
  MODIFIED  LDAUU    : 3.5 0 0 → 4.0 0 0  (+0.5 eV on site 1)
  MODIFIED  EDIFF    : 1e-05 → 1e-06  (tighter by 10×)
  ADDED     LWAVE    : .TRUE.
  DELETED   LCHARG   : (was .FALSE.)
  UNCHANGED : ENCUT=520, ISMEAR=0, SIGMA=0.05, ... (18 params)

── POSCAR diff ─────────────────────────────────
  Formula    : Li4Co4O8 → Li4Co3FeO8  (substitution: Co→Fe ×1)
  Spacegroup : R-3m → R-3m  (unchanged)
  Lattice a  : 2.830 → 2.845 Å  (+0.015 Å, +0.5%)
  Lattice c  : 14.050 → 14.120 Å  (+0.070 Å, +0.5%)
  Volume     : 97.45 → 98.12 Å³  (+0.69 Å³, +0.7%)
  Natoms     : 16 → 16  (unchanged)
  Coord RMSD : 0.0234 Å  (over 16 matched sites)

── KPOINTS diff ────────────────────────────────
  Grid       : 4×4×4 → 6×6×6  (+50% per axis)
  Est. k-pts : 64 → 216  (×3.4)

── Output summary diff ─────────────────────────
  Energy     : -45.23 → -44.87 eV  (+0.36 eV, +0.8%)
  Per atom   : -2.827 → -2.805 eV/atom
  Converged  : ✅ → ✅
  Ionic steps: 23 → 31  (+8)
  Max force  : 0.0012 → 0.0008 eV/Å  (-33%, improved)
```

**JSON 输出**：

```json
{
  "rev1": "dead000...",
  "rev2": "abc1234...",
  "files_changed": 3,
  "changes": {
    "INCAR": {
      "type": "semantic",
      "modified": [
        {"key": "LDAUU", "old": [3.5, 0, 0], "new": [4.0, 0, 0], "delta": [0.5, 0, 0]},
        {"key": "EDIFF", "old": 1e-5, "new": 1e-6, "delta_factor": 10}
      ],
      "added": [{"key": "LWAVE", "value": true}],
      "deleted": [{"key": "LCHARG", "old_value": false}]
    },
    "POSCAR": {
      "type": "semantic",
      "formula_change": {"old": "Li4Co4O8", "new": "Li4Co3FeO8"},
      "lattice_a_change": {"old": 2.830, "new": 2.845, "delta": 0.015},
      "coord_rmsd": 0.0234
    }
  }
}
```

---

### `chemvcs reproduce`

**用途**：还原历史 commit 的所有输入文件

**签名**：
```bash
chemvcs reproduce [OPTIONS] <REVISION>
```

**参数**：

| 参数 | 类型 | 必需 | 说明 |
|------|------|------|------|
| `REVISION` | string | ✅ | commit hash 或 tag（支持短 hash，如 `abc1234`） |

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `-o, --output-dir` | path | `reproduce_<rev>/` | 还原目标目录 |
| `--verify-potcar` | flag | true | 检查赝势文件的哈希是否匹配 |
| `--verify-env` | flag | true | 比对环境信息并输出差异 |
| `--no-verify` | flag | false | 跳过所有校验（仅还原文件） |

**行为**：

1. 检查项目是否已初始化
2. 解析 REVISION → 获取 commit 对象
3. 创建输出目录（若已存在 → 报错，除非 `--force`）
4. 对每个文件：
   - 从 `.chemvcs/objects/` 读取 blob
   - 计算 SHA-256 并与 commit 记录对比（完整性校验）
   - 若不匹配 → 退出码 3
   - 写入输出目录
5. 若 `--verify-potcar`：
   - 检查赝势路径是否存在
   - 计算赝势文件的实际 hash
   - 比对并输出结果（✅ / ⚠️ / ❌）
6. 若 `--verify-env`：
   - 检测当前 VASP 版本、module list
   - 与 commit 记录比对
   - 输出差异
7. 打印还原摘要

**输出示例**：

```
Reproducing commit abc1234...

Restoring files to reproduce_abc1234/
  ✓ INCAR       (SHA: fc1a2b3c... verified)
  ✓ POSCAR      (SHA: ab12cd34... verified)
  ✓ KPOINTS     (SHA: ee44ff55... verified)

POTCAR verification:
  ✅ Li_sv      (/opt/vasp/pp/PBE/Li_sv/POTCAR, hash matched)
  ✅ Co         (/opt/vasp/pp/PBE/Co/POTCAR, hash matched)
  ✅ O          (/opt/vasp/pp/PBE/O/POTCAR, hash matched)

Environment check:
  VASP version:    6.3.2 (recorded: 6.3.2) ✓
  Loaded modules:  vasp/6.3.2, intel/2023 (recorded: same) ✓
  Python version:  3.9.7 (recorded: 3.9.7) ✓
  Hostname:        hpc-login01 (recorded: hpc-login01) ✓

✓ Successfully reproduced commit abc1234
  Output: reproduce_abc1234/
  Ready to run VASP in this directory
```

**警告示例**：

```
POTCAR verification:
  ⚠️  Li_sv      (path not found: /opt/vasp/pp/PBE/Li_sv/POTCAR)
  ❌ Co         (hash mismatch: expected a1b2c3..., got e5f6a7...)
  ✅ O          (hash matched)

Environment check:
  ⚠️  VASP version:    6.4.0 (recorded: 6.3.2) — VERSION MISMATCH
  ⚠️  Loaded modules:  vasp/6.4.0 (recorded: vasp/6.3.2) — DIFFERENT
```

**JSON 输出**：

```json
{
  "success": true,
  "commit_id": "abc1234...",
  "output_dir": "reproduce_abc1234",
  "files_restored": 3,
  "integrity_checks": {
    "all_passed": true,
    "files": [
      {"path": "INCAR", "hash_match": true},
      {"path": "POSCAR", "hash_match": true}
    ]
  },
  "potcar_verification": {
    "all_matched": true,
    "elements": [
      {"element": "Li_sv", "path_exists": true, "hash_match": true},
      {"element": "Co", "path_exists": true, "hash_match": true}
    ]
  },
  "environment_verification": {
    "vasp_version_match": true,
    "modules_match": true,
    "differences": []
  }
}
```

---

### `chemvcs status`

**用途**：显示暂存区和工作区状态

**签名**：
```bash
chemvcs status [OPTIONS]
```

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--short` | flag | false | 精简输出（类似 `git status --short`） |

**输出示例**：

```
On commit abc1234
  PBE+U, U=4.0eV, kpoints 6x6x6

Staged for commit:
  modified:  INCAR
  modified:  POSCAR

Unstaged changes:
  modified:  KPOINTS

Ignored files:
  WAVECAR    (2.1 GB)
  CHGCAR     (345 MB)

Tip: Use 'chemvcs add <file>' to stage changes
     Use 'chemvcs commit' to save staged changes
```

---

### `chemvcs fsck`

**用途**：检查数据完整性

**签名**：
```bash
chemvcs fsck [OPTIONS]
```

**选项**：

| 选项 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--repair` | flag | false | 尝试修复可恢复的错误（如重建 SQLite 索引） |

**行为**：

1. 遍历所有 commit JSON
2. 验证每个 commit 的 hash
3. 检查引用的 blob 是否存在
4. 计算 blob 的 SHA-256 并与 commit 记录比对
5. 检查 SQLite 索引是否包含所有 commit
6. 输出完整性报告

**输出示例**：

```
Checking repository integrity...

[1/3] Verifying commits...            ✓ 42 commits OK
[2/3] Verifying blobs...               ✓ 128 blobs OK
[3/3] Verifying SQLite index...        ✓ Index consistent

✓ Repository is healthy
```

**错误示例**：

```
Checking repository integrity...

[1/3] Verifying commits...            ⚠️  1 error
  ✗ Commit dead000: blob fc1a2b not found

[2/3] Verifying blobs...               ✓ 127 blobs OK
[3/3] Verifying SQLite index...        ⚠️  Out of sync

✗ Found 2 issues
  Run 'chemvcs fsck --repair' to attempt fixes

Exit code: 3
```

---

## 全局选项（所有命令通用）

| 选项 | 说明 |
|------|------|
| `--help` | 显示命令帮助 |
| `--version` | 显示 ChemVCS 版本 |
| `--repo <path>` | 指定项目目录（默认：当前目录） |
| `--no-color` | 禁用彩色输出 |
| `--debug` | 输出调试信息（含 Python traceback） |

---

## 环境变量

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `CHEMVCS_REPO` | `.` | 项目根目录 |
| `CHEMVCS_EDITOR` | `$EDITOR` or `vim` | commit message 编辑器 |
| `CHEMVCS_NO_COLOR` | `0` | 设为 `1` 禁用彩色输出 |
| `CHEMVCS_DEBUG` | `0` | 设为 `1` 启用调试模式 |

---

*本规格为 MVP 阶段命令集。Phase 2/3 将增加 `submit` / `monitor` / `push` / `pull` 等命令。*
