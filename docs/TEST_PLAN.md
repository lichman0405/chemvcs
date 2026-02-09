# ChemVCS 测试策略

> v1.0 | 2026-02-09  
> 范围：MVP 测试计划（Phase 1 核心版本控制）

---

## 测试目标

| 层级 | 覆盖率目标 | 度量方式 |
|------|-----------|----------|
| **单元测试** | ≥85% 行覆盖率 | pytest-cov |
| **集成测试** | 100% 核心用户故事覆盖 | 5 个 US（US-01 至 US-05）各 ≥1 个端到端场景 |
| **HPC 环境测试** | 3 种文件系统验证 | Lustre / GPFS / NFS 各手动测试并记录 |
| **性能基准** | 关键操作满足约束 | commit <5s, reproduce <3s, diff <1s |

---

## 1. 单元测试策略

### 1.1 测试框架

```python
# pyproject.toml
[tool.pytest.ini_options]
minversion = "7.0"
testpaths = ["tests"]
python_files = ["test_*.py"]
python_classes = ["Test*"]
python_functions = ["test_*"]
addopts = [
    "--strict-markers",
    "--cov=chemvcs",
    "--cov-report=html",
    "--cov-report=term-missing:skip-covered",
    "--cov-fail-under=85",
]

[tool.coverage.run]
source = ["chemvcs"]
omit = ["*/tests/*", "*/__pycache__/*"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise NotImplementedError",
    "if TYPE_CHECKING:",
    "if __name__ == .__main__.:",
]
```

### 1.2 模块测试范围与关键用例

#### 存储层（`chemvcs/storage/`）

**`object_store.py` — Content-addressable blob 存储**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_write_blob_creates_file` | blob 写入 | `.chemvcs/objects/<hash前2>/<hash后62>` 文件创建成功 |
| `test_write_blob_deduplication` | 去重 | 同内容写入两次只创建一个 blob |
| `test_write_blob_atomic_rename` | 原子性 | 先写 tmp，再 rename；中途失败不产生半成品 |
| `test_write_blob_gzip_large_file` | 大文件压缩 | ≥200MB 文件自动 gzip，后缀 `.gz` |
| `test_read_blob_verifies_hash` | 完整性 | 读取时重算 hash，不匹配抛异常 |
| `test_write_blob_disk_full` | 错误处理 | 磁盘满时抛 `OSError`，清理 tmp 文件 |

**`metadata_db.py` — SQLite 元数据管理**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_init_db_creates_tables` | 表结构初始化 | 6 张表（commits/files/semantic/output/environment/blob_refs）存在 |
| `test_insert_commit_transaction` | 事务 | commit + files 插入在同一事务中，失败时回滚 |
| `test_db_wal_mode_detection` | 文件系统兼容 | 检测 WAL 失败时自动退回 DELETE 模式 |
| `test_concurrent_write_lock` | 并发锁 | 两个进程同时写，第二个等待或超时报错 |
| `test_rebuild_index_from_commits` | 索引重建 | 删除 DB 后从 `commits/` 目录完全重建 |

**`commit_builder.py` — Commit 对象构建**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_compute_commit_hash_deterministic` | hash 可重现 | 相同内容生成相同 hash |
| `test_commit_json_canonical_serialization` | JSON 规范化 | 字典键排序，浮点数精度一致 |
| `test_commit_parent_linkage` | 线性历史 | 新 commit 的 `parent_id` 指向当前 HEAD |

#### 语义解析层（`chemvcs/parsers/`）

**`incar_parser.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_parse_incar_standard_format` | pymatgen 正常路径 | 解析出 key-value dict |
| `test_parse_incar_with_comments` | 注释处理 | 忽略 `#` 和 `!` 注释 |
| `test_parse_incar_multiline_values` | 边缘用例 | `3*0.0` 解析为 `[0.0, 0.0, 0.0]` |
| `test_parse_incar_fallback_regex` | 回退机制 | pymatgen 失败时用正则提取 key=value |
| `test_incar_diff_numeric_change` | 语义 diff | `ENCUT: 400 → 520` 正确标注 `+120 eV` |
| `test_incar_diff_boolean_normalization` | 类型归一化 | `.TRUE.` / `T` / `.T.` 统一为 `True` |

**`poscar_parser.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_parse_poscar_vasp5_format` | VASP 5 格式 | 含元素行，正确解析 |
| `test_parse_poscar_vasp4_format` | VASP 4 格式 | 无元素行，依赖原子数推断 |
| `test_poscar_diff_lattice_change` | 晶格参数 diff | `a: 2.830 → 2.845 Å (+0.015)` |
| `test_poscar_diff_coord_rmsd` | 坐标 RMSD | 计算匹配原子的 RMSD |
| `test_poscar_diff_formula_change` | 化学式变化 | 检测 `Li4Co4O8 → Li4Co3FeO8` |

**`outcar_parser.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_extract_total_energy` | 能量提取 | 从 OUTCAR 最后一次 SCF 提取 `E0=` |
| `test_extract_convergence_status` | 收敛标记 | 检测 `reached required accuracy` |
| `test_extract_truncated_outcar` | 异常文件 | 截断的 OUTCAR 标记 `is_converged=false` |
| `test_extract_warnings` | 警告捕获 | 识别 `BRMIX: very serious problems` 等警告 |

#### CLI 层（`chemvcs/cli/`）

**`init.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_init_creates_directory_structure` | 目录创建 | `.chemvcs/objects/`, `.chemvcs/commits/` 等存在 |
| `test_init_double_init_fails` | 幂等性 | 第二次 init 退出码 1 |
| `test_init_force_overwrites` | --force 选项 | 覆盖已有 `.chemvcs/`，警告用户 |

**`add.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_add_file_to_staging` | 暂存 | `staging/manifest.json` 包含文件条目 |
| `test_add_respects_ignore` | .chemvcsignore | WAVECAR 被跳过 |
| `test_add_potcar_auto_reference` | POTCAR 特殊处理 | `is_reference=true`，打印警告 |

**`commit.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_commit_empty_staging_fails` | 前置检查 | 退出码 1，提示 staging 为空 |
| `test_commit_writes_blobs_and_json` | 端到端 | blob + commit JSON + HEAD 更新 + SQLite 插入 |

**`diff.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_diff_incar_semantic_output` | 语义 diff | 输出包含 `MODIFIED LDAUU: 3.5 → 4.0 (+0.5)` |
| `test_diff_json_format` | JSON 输出 | `--format json` 返回可解析的 JSON |

**`reproduce.py`**

| 测试用例 | 覆盖点 | 预期结果 |
|----------|--------|----------|
| `test_reproduce_restores_files` | 文件还原 | 输出目录包含所有输入文件 |
| `test_reproduce_hash_verification` | 完整性校验 | SHA-256 不匹配时退出码 3 |

### 1.3 Fixture 与测试数据

```python
# tests/conftest.py
import pytest
from pathlib import Path
import tempfile

@pytest.fixture
def temp_repo(tmp_path: Path):
    """创建临时测试仓库"""
    repo = tmp_path / "test_repo"
    repo.mkdir()
    (repo / "INCAR").write_text("ENCUT = 520\nISMEAR = 0\n")
    (repo / "POSCAR").write_text(POSCAR_VASP5_SAMPLE)
    (repo / "KPOINTS").write_text("Automatic\n0\nGamma\n4 4 4\n")
    return repo

@pytest.fixture
def sample_incar():
    return {
        "ENCUT": 520,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "LDAUU": [3.5, 0, 0]
    }

# 真实 VASP 文件样本（匿名化）
POSCAR_VASP5_SAMPLE = """Li4 Co4 O8
1.0
2.830000   0.000000   0.000000
-1.415000   2.450794   0.000000
0.000000   0.000000  14.050000
Li Co O
4 4 8
direct
0.000000  0.000000  0.000000 Li
...
"""
```

---

## 2. 集成测试策略

### 2.1 用户故事端到端测试

每个 US（US-01 至 US-05）对应 ≥1 个集成测试，覆盖完整用户旅程。

#### US-01：初始化项目

```python
# tests/integration/test_init_workflow.py

def test_us01_init_project(tmp_path):
    """
    Given: 一个包含 VASP 文件的空目录
    When:  执行 chemvcs init
    Then:  .chemvcs/ 目录创建成功，检测到 VASP 文件
    """
    # Arrange
    repo = tmp_path / "lco_project"
    repo.mkdir()
    (repo / "INCAR").write_text("ENCUT = 520\n")
    
    # Act
    result = subprocess.run(
        ["chemvcs", "init"],
        cwd=repo,
        capture_output=True,
        text=True
    )
    
    # Assert
    assert result.returncode == 0
    assert (repo / ".chemvcs" / "metadata.db").exists()
    assert "Detected VASP files" in result.stdout
    assert "INCAR" in result.stdout
```

#### US-02：快照提交

```python
def test_us02_commit_workflow(initialized_repo):
    """
    Given: 已初始化的项目，有 INCAR/POSCAR 文件
    When:  执行 add . && commit -m "msg"
    Then:  commit 成功，blob 和 commit JSON 创建，SQLite 更新
    """
    # Act
    subprocess.run(["chemvcs", "add", "."], cwd=initialized_repo, check=True)
    result = subprocess.run(
        ["chemvcs", "commit", "-m", "Initial"],
        cwd=initialized_repo,
        capture_output=True,
        text=True
    )
    
    # Assert
    assert result.returncode == 0
    assert "Committed" in result.stdout
    
    # 验证 HEAD 文件存在且为 40 字符 hex
    head = (initialized_repo / ".chemvcs" / "HEAD").read_text().strip()
    assert len(head) == 64  # SHA-256 hex
    
    # 验证 commit JSON 存在
    commit_path = initialized_repo / ".chemvcs" / "commits" / head[:2] / f"{head[2:]}.json"
    assert commit_path.exists()
    
    # 验证 blob 存在
    objects_count = len(list((initialized_repo / ".chemvcs" / "objects").rglob("*")))
    assert objects_count >= 2  # 至少 INCAR + POSCAR
```

#### US-03：查看历史

```python
def test_us03_log_workflow(repo_with_commits):
    """
    Given: 已有 3 个 commit 的仓库
    When:  执行 chemvcs log
    Then:  按时间倒序显示 3 条记录，含语义摘要
    """
    result = subprocess.run(
        ["chemvcs", "log"],
        cwd=repo_with_commits,
        capture_output=True,
        text=True
    )
    
    assert result.returncode == 0
    assert result.stdout.count("commit") == 3
    assert "INCAR:" in result.stdout
    assert "POSCAR:" in result.stdout
```

#### US-04：语义差异比较

```python
def test_us04_diff_workflow(repo_with_two_commits):
    """
    Given: 仓库有两个 commit，INCAR 中 LDAUU 参数不同
    When:  执行 chemvcs diff HEAD~1 HEAD
    Then:  输出语义 diff，标注 LDAUU 变化量
    """
    result = subprocess.run(
        ["chemvcs", "diff", "HEAD~1", "HEAD"],
        cwd=repo_with_two_commits,
        capture_output=True,
        text=True
    )
    
    assert result.returncode == 0
    assert "INCAR diff" in result.stdout
    assert "LDAUU" in result.stdout
    assert "→" in result.stdout  # 变化箭头
    assert "eV" in result.stdout  # 单位
```

#### US-05：复现历史计算

```python
def test_us05_reproduce_workflow(repo_with_commits):
    """
    Given: 仓库有历史 commit
    When:  执行 chemvcs reproduce <commit_hash>
    Then:  文件还原成功，hash 校验通过，赝势路径打印
    """
    # 获取第一个 commit hash
    log_result = subprocess.run(
        ["chemvcs", "log", "--format", "json", "-n", "1"],
        cwd=repo_with_commits,
        capture_output=True,
        text=True
    )
    commits = json.loads(log_result.stdout)["commits"]
    target_hash = commits[0]["short_id"]
    
    # Act
    result = subprocess.run(
        ["chemvcs", "reproduce", target_hash],
        cwd=repo_with_commits,
        capture_output=True,
        text=True
    )
    
    # Assert
    assert result.returncode == 0
    assert f"reproduce_{target_hash}" in result.stdout
    
    # 验证文件存在
    output_dir = repo_with_commits / f"reproduce_{target_hash}"
    assert output_dir.exists()
    assert (output_dir / "INCAR").exists()
    
    # 验证 hash
    assert "verified" in result.stdout.lower()
```

### 2.2 错误场景测试

| 场景 | 预期行为 | 退出码 |
|------|----------|--------|
| 未初始化时执行 `add` | 输出 "项目未初始化"，提示先 `init` | 1 |
| staging 为空时 `commit` | 输出 "无文件可提交" | 1 |
| `reproduce` 时 blob 损坏 | 输出损坏文件列表，退出 | 3 |
| 磁盘满时 `commit` | 输出磁盘空间不足，清理 tmp 文件 | 2 |
| 并发 `commit`（两个进程） | 第二个进程等待锁或超时报错 | 2 |

---

## 3. HPC 环境测试

### 3.1 测试矩阵

| 文件系统 | 测试站点 | 调度器 | SQLite 模式 | 状态 |
|----------|---------|--------|------------|------|
| Lustre | 天河二号 / 神威太湖之光 | SLURM | WAL → DELETE（回退） | 待验证 |
| GPFS | 国家超算济南 | PBS | WAL 或 DELETE | 待验证 |
| NFS | 学校小型集群 | SLURM | DELETE（NFS 不支持 WAL） | 待验证 |

### 3.2 测试清单（每种文件系统）

```bash
# 在真实 HPC 环境执行的测试脚本
#!/bin/bash
# tests/hpc/test_on_lustre.sh

set -e

echo "=== ChemVCS HPC 环境测试 ==="
echo "文件系统: $(df -T . | tail -1 | awk '{print $2}')"
echo "主机: $(hostname)"

# 1. 安装测试
pip install --user chemvcs
chemvcs --version || exit 1

# 2. 基本流程测试
cd /scratch/$USER/chemvcs_test
chemvcs init
echo "ENCUT = 520" > INCAR
chemvcs add INCAR
chemvcs commit -m "test commit" || exit 1

# 3. 并发测试（两个 shell 同时 commit）
(sleep 1; chemvcs commit -m "concurrent 1") &
(sleep 1; chemvcs commit -m "concurrent 2") &
wait
# 预期：一个成功，一个失败或等待

# 4. SQLite 模式检测
sqlite3 .chemvcs/metadata.db "PRAGMA journal_mode;" | tee journal_mode.txt
# 记录是 WAL 还是 DELETE

# 5. 大文件测试（若有空间）
dd if=/dev/zero of=large.dat bs=1M count=500
chemvcs add large.dat --force
chemvcs commit -m "large file"
ls -lh .chemvcs/objects/*/* | grep large

echo "✓ HPC 环境测试通过"
```

### 3.3 问题记录模板

| 日期 | 文件系统 | 问题描述 | 复现步骤 | 临时解决方案 | 长期修复计划 |
|------|---------|---------|---------|-------------|-------------|
| 2026-02-15 | Lustre | WAL 模式导致 `disk I/O error` | 在 Lustre 上执行 commit | 自动回退到 DELETE 模式 | 在文档中明确说明 |

---

## 4. 性能基准测试

### 4.1 基准场景与目标

| 操作 | 场景 | 目标延迟 | 测量方式 |
|------|------|---------|----------|
| `init` | 空目录 | <2s (SSD) / <5s (Lustre) | `time chemvcs init` |
| `add` | 5 个 VASP 文件（含 POTCAR） | <3s | `time chemvcs add .` |
| `commit` | 5 个输入文件 | <5s (SSD) / <10s (Lustre) | `time chemvcs commit -m "msg"` |
| `log` | 100 条 commit | <1s | `time chemvcs log -n 100` |
| `diff` | 两个 commit 的 INCAR/POSCAR | <1s | `time chemvcs diff v1 v2` |
| `reproduce` | 5 个输入文件 | <3s | `time chemvcs reproduce <hash>` |

### 4.2 性能测试脚本

```python
# tests/performance/benchmark.py
import time
import subprocess
from pathlib import Path

def benchmark_commit(repo: Path, iterations: int = 10):
    """基准测试 commit 操作"""
    durations = []
    for i in range(iterations):
        # 修改文件触发新 commit
        (repo / "INCAR").write_text(f"ENCUT = {520 + i}\n")
        subprocess.run(["chemvcs", "add", "INCAR"], cwd=repo, check=True)
        
        start = time.perf_counter()
        subprocess.run(
            ["chemvcs", "commit", "-m", f"iter {i}"],
            cwd=repo,
            check=True
        )
        end = time.perf_counter()
        durations.append(end - start)
    
    print(f"Commit avg: {sum(durations)/len(durations):.3f}s")
    print(f"Commit p95: {sorted(durations)[int(0.95*len(durations))]:.3f}s")
    assert max(durations) < 5.0, "Commit timeout (>5s)"
```

### 4.3 内存基准

```bash
# 测试 commit 操作的峰值内存
/usr/bin/time -v chemvcs commit -m "test" 2>&1 | grep "Maximum resident"
# 目标：<200 MB
```

---

## 5. 测试自动化与 CI

### 5.1 GitHub Actions 工作流

```yaml
# .github/workflows/test.yml
name: Test

on: [push, pull_request]

jobs:
  unit-tests:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: ['3.8', '3.9', '3.10', '3.11', '3.12']
    
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}
      
      - name: Install dependencies
        run: |
          pip install -e ".[dev]"
      
      - name: Run unit tests
        run: |
          pytest tests/unit --cov=chemvcs --cov-report=xml
      
      - name: Upload coverage
        uses: codecov/codecov-action@v3
        with:
          files: ./coverage.xml

  integration-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      
      - name: Install
        run: pip install -e ".[dev]"
      
      - name: Run integration tests
        run: pytest tests/integration -v

  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      
      - name: Lint with ruff
        run: |
          pip install ruff
          ruff check chemvcs/
      
      - name: Type check with mypy
        run: |
          pip install mypy
          mypy chemvcs/ --strict
```

### 5.2 Pre-commit Hooks

```yaml
# .pre-commit-config.yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.9
    hooks:
      - id: ruff
        args: [--fix, --exit-non-zero-on-fix]
  
  - repo: https://github.com/pre-commit/mirrors-mypy
    rev: v1.8.0
    hooks:
      - id: mypy
        additional_dependencies: [types-all]
  
  - repo: local
    hooks:
      - id: pytest-quick
        name: Quick unit tests
        entry: pytest tests/unit -x --tb=short
        language: system
        pass_filenames: false
        always_run: true
```

---

## 6. 测试数据管理

### 6.1 测试数据来源

| 数据类型 | 来源 | 匿名化策略 |
|---------|------|-----------|
| VASP 输入文件 | Materials Project API | 公开数据，无需匿名化 |
| VASP 输出文件 | 真实计算结果（合作者提供） | 去除机构名、路径；仅保留数值 |
| 边缘用例 | 手动构造 | 包含非标准格式、截断文件、超大值 |

### 6.2 测试数据目录

```
tests/
├── fixtures/
│   ├── vasp_inputs/
│   │   ├── standard/
│   │   │   ├── INCAR
│   │   │   ├── POSCAR
│   │   │   └── KPOINTS
│   │   ├── edge_cases/
│   │   │   ├── INCAR_with_comments
│   │   │   ├── POSCAR_vasp4
│   │   │   └── POSCAR_cartesian
│   ├── vasp_outputs/
│   │   ├── OUTCAR_converged
│   │   ├── OUTCAR_truncated
│   │   └── vasprun.xml.gz
│   └── repos/
│       └── sample_repo.tar.gz  # 预构建的测试仓库
```

---

## 7. 回归测试

### 7.1 回归测试触发条件

- 每次 PR 合并到 `main` 分支
- 每周自动运行（GitHub Actions scheduled workflow）
- 发布前手动触发（Release Candidate 验证）

### 7.2 回归套件内容

```python
# tests/regression/test_backward_compatibility.py

def test_v0_1_repos_readable_by_v0_2():
    """确保新版本能读取旧版本创建的仓库"""
    # 使用归档的 v0.1 仓库
    old_repo = extract_fixture("repos/v0.1_sample.tar.gz")
    result = subprocess.run(
        ["chemvcs", "log"],
        cwd=old_repo,
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "commit" in result.stdout
```

---

## 8. 缺陷管理

### 8.1 缺陷分级

| 等级 | 定义 | SLA |
|------|------|-----|
| P0 | 数据丢失 / 数据损坏 / 无法使用核心功能 | 24h 内修复 |
| P1 | 性能严重下降 / 用户体验严重受损 | 1 周内修复 |
| P2 | 边缘用例失败 / 非核心功能异常 | 2 周内修复或降级为改进 |
| P3 | 文档错误 / 优化建议 | 下个版本修复 |

### 8.2 Issue 模板

```markdown
## Bug Report

**Environment:**
- ChemVCS version: 
- Python version: 
- OS: 
- File system (if HPC): 

**Steps to Reproduce:**
1. 
2. 
3. 

**Expected Behavior:**

**Actual Behavior:**

**Error Message:**
```

---

*本测试策略覆盖 MVP 范围。随着 Phase 2/3 功能（HPC 集成、远程同步）的加入，将持续扩展测试矩阵。*
