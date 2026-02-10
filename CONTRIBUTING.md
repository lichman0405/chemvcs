# ChemVCS 开发规范

> v1.0 | 2026-02-09  
> 适用范围：所有贡献者（核心团队 + 社区）

---

## 1. 开发环境设置

### 1.1 前置要求

- Python ≥3.8（推荐 3.10 或 3.11）
- Git ≥2.30
- 操作系统：Linux / macOS / Windows（WSL2 或 PowerShell）

### 1.2 初次设置

```bash
# 1. Fork 仓库（如果你是外部贡献者）
# 访问 https://github.com/lichman0405/chemvcs 并 fork

# 2. Clone 到本地
git clone https://github.com/<your-username>/chemvcs.git
cd chemvcs

# 3. 添加上游远程仓库
git remote add upstream https://github.com/lichman0405/chemvcs.git

# 4. 创建虚拟环境（推荐）
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# 5. 安装开发依赖
pip install -e ".[dev]"

# 6. 安装 pre-commit hooks
pre-commit install

# 7. 验证安装
chemvcs --version
pytest --version
ruff --version
```

### 1.3 开发依赖

```toml
# pyproject.toml 中的 dev 依赖
[project.optional-dependencies]
dev = [
    "pytest >=7.4.0",
    "pytest-cov >=4.1.0",
    "pytest-timeout >=2.1.0",
    "ruff >=0.1.9",
    "mypy >=1.8.0",
    "black >=23.12.0",
    "pre-commit >=3.5.0",
    "ipython >=8.12.0",  # REPL 调试
    "hypothesis >=6.92.0",  # 属性测试
]
```

---

## 2. 代码风格

### 2.1 自动化工具配置

#### Ruff（Linter + Formatter）

```toml
# pyproject.toml
[tool.ruff]
target-version = "py38"
line-length = 100
indent-width = 4

[tool.ruff.lint]
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # pyflakes
    "I",   # isort
    "B",   # flake8-bugbear
    "C4",  # flake8-comprehensions
    "UP",  # pyupgrade
    "SIM", # flake8-simplify
]
ignore = [
    "E501",  # line too long (由 formatter 处理)
    "B008",  # do not perform function calls in argument defaults
]

[tool.ruff.format]
quote-style = "double"
indent-style = "space"
skip-magic-trailing-comma = false

[tool.ruff.lint.isort]
known-first-party = ["chemvcs"]
```

**使用方式：**

```bash
# 检查代码
ruff check chemvcs/

# 自动修复
ruff check --fix chemvcs/

# 格式化
ruff format chemvcs/
```

#### MyPy（类型检查）

```toml
# pyproject.toml
[tool.mypy]
python_version = "3.8"
warn_return_any = true
warn_unused_configs = true
disallow_untyped_defs = true
disallow_any_generics = true
check_untyped_defs = true
no_implicit_optional = true
warn_redundant_casts = true
warn_unused_ignores = true
warn_no_return = true
strict_equality = true

[[tool.mypy.overrides]]
module = "pymatgen.*"
ignore_missing_imports = true
```

**使用方式：**

```bash
mypy chemvcs/
```

### 2.2 命名约定

| 类型 | 规则 | 示例 |
|------|------|------|
| **模块** | 小写 + 下划线 | `object_store.py`, `incar_parser.py` |
| **类** | PascalCase | `ObjectStore`, `IncarParser` |
| **函数/方法** | 小写 + 下划线 | `write_blob()`, `parse_incar()` |
| **常量** | 大写 + 下划线 | `DEFAULT_IGNORE_PATTERNS`, `MAX_FILE_SIZE` |
| **私有成员** | 前缀 `_` | `_compute_hash()`, `_internal_state` |
| **类型变量** | PascalCase + `T` | `FileT`, `CommitT` |

### 2.3 类型注解

**所有公开 API 必须包含类型注解：**

```python
# ✅ 正确
def write_blob(content: bytes, objects_dir: Path) -> str:
    """写入 blob，返回 SHA-256 hash。
    
    Args:
        content: 文件内容
        objects_dir: objects/ 目录路径
    
    Returns:
        40 字符的 SHA-256 hex string
    
    Raises:
        OSError: 写入失败（权限/磁盘满）
    """
    ...

# ❌ 错误（缺少类型注解）
def write_blob(content, objects_dir):
    ...
```

**使用 `typing` 模块：**

```python
from typing import Optional, Dict, List, Tuple, Union
from pathlib import Path

def parse_incar(path: Path) -> Dict[str, Union[int, float, bool, str, List[float]]]:
    ...

def get_commit(commit_id: str) -> Optional[Dict]:
    """找不到时返回 None"""
    ...
```

### 2.4 Docstring 规范（Google Style）

```python
class IncarParser:
    """INCAR 文件的语义解析器。
    
    使用 pymatgen 作为主解析器，失败时退回正则表达式。
    
    Attributes:
        fallback_enabled: 是否启用正则回退模式
        param_meta: 参数元数据（类型、单位等）
    
    Example:
        >>> parser = IncarParser()
        >>> params = parser.parse(Path("INCAR"))
        >>> params["ENCUT"]
        520.0
    """
    
    def parse(self, path: Path) -> Dict[str, Any]:
        """解析 INCAR 文件。
        
        Args:
            path: INCAR 文件路径
        
        Returns:
            参数字典，键为参数名，值为归一化后的值
        
        Raises:
            FileNotFoundError: 文件不存在
            ValueError: 文件内容无法解析
        """
        ...
```

---

## 3. 分支策略

### 3.1 分支模型

```
main (protected)          ← 稳定版本，每次发布打 tag
  ↑
  merge ← dev             ← 开发分支，集成所有功能
           ↑
           merge ← feature/issue-123-incar-parser
                ← feature/issue-124-diff-poscar
                ← fix/issue-125-potcar-hash
```

### 3.2 分支命名规则

| 类型 | 格式 | 示例 |
|------|------|------|
| 功能开发 | `feature/<issue-号>-<简短描述>` | `feature/42-semantic-diff` |
| Bug 修复 | `fix/<issue-号>-<简短描述>` | `fix/55-potcar-hash-mismatch` |
| 文档 | `docs/<描述>` | `docs/update-cli-spec` |
| 重构 | `refactor/<描述>` | `refactor/storage-layer` |
| 发布准备 | `release/v<版本号>` | `release/v0.1.0` |

### 3.3 工作流程

```bash
# 1. 更新本地 dev 分支
git checkout dev
git pull upstream dev

# 2. 创建功能分支（从 dev 拉取）
git checkout -b feature/42-semantic-diff

# 3. 开发 + commit
# ... 编码 ...
git add chemvcs/parsers/incar_diff.py
git commit -m "feat: implement INCAR semantic diff

- Add IncarDiffer class
- Support numeric delta calculation
- Add unit tests (coverage: 92%)

Closes #42"

# 4. 本地测试
pytest tests/unit/parsers/
ruff check chemvcs/
mypy chemvcs/

# 5. 推送到 fork
git push origin feature/42-semantic-diff

# 6. 在 GitHub 上创建 PR（target: upstream/dev）
```

---

## 4. Commit 规范

### 4.1 Commit Message 格式（Conventional Commits）

```
<type>(<scope>): <subject>

<body>

<footer>
```

**Type（必需）：**

| Type | 说明 | 示例 |
|------|------|------|
| `feat` | 新功能 | `feat(cli): add chemvcs reproduce command` |
| `fix` | Bug 修复 | `fix(storage): handle disk full error in write_blob` |
| `docs` | 文档 | `docs: update CLI_SPEC with --format option` |
| `test` | 测试 | `test(parsers): add edge case for INCAR comments` |
| `refactor` | 重构 | `refactor(diff): extract common diff logic` |
| `perf` | 性能优化 | `perf(commit): cache blob hash to avoid recalc` |
| `chore` | 构建/工具 | `chore: upgrade ruff to 0.2.0` |

**Scope（可选）：** 受影响的模块，如 `cli`, `storage`, `parsers`, `diff`

**Subject：** 简短描述（≤50 字符），使用祈使句（add, fix, update）

**Body（可选）：** 详细说明（wrap at 72 chars）

**Footer（可选）：** Issue 引用，Breaking Changes

**示例：**

```
feat(parsers): add KPOINTS parser with grid mode support

- Implement KpointsParser class using pymatgen
- Support Automatic (Gamma/Monkhorst-Pack) mode
- Add diff logic for grid size comparison
- Test coverage: 88%

Closes #42
```

### 4.2 Commit 原子性

- 每个 commit 应该是**独立可测试的**（不破坏 CI）
- 单一职责：一个 commit 只做一件事（实现功能 vs 重构 vs 修复 bug）
- 避免"WIP"、"fix typo"等无信息 commit（可在本地 squash 后再推送）

---

## 5. Pull Request 流程

### 5.1 PR 检查清单

**提交 PR 前必须确认：**

- [ ] 所有测试通过（`pytest tests/`）
- [ ] 代码风格检查通过（`ruff check`）
- [ ] 类型检查通过（`mypy chemvcs/`）
- [ ] 新增代码有单元测试（覆盖率 ≥85%）
- [ ] 更新了相关文档（若改动了 API）
- [ ] PR 描述清晰且关联了 Issue（`Closes #123`）

### 5.2 PR 模板

```markdown
## Description
<!-- 简要描述这个 PR 做了什么 -->

Closes #<issue-number>

## Changes
<!-- 列出主要变更点 -->

- Added X
- Fixed Y
- Refactored Z

## Testing
<!-- 如何测试这些变更？ -->

- [ ] Unit tests added (coverage: __%)
- [ ] Integration test: <scenario>
- [ ] Manual testing on HPC: <cluster-name>

## Checklist
- [ ] Tests pass locally
- [ ] Ruff/mypy checks pass
- [ ] Documentation updated
- [ ] Breaking changes noted (if any)

## Screenshots/Output (if applicable)
```

### 5.3 Code Review 要求

**所有 PR 必须经过 ≥1 人 review 才能合并。**

**Reviewer checklist：**

- [ ] 代码逻辑正确
- [ ] 边缘情况处理（空文件、大文件、并发等）
- [ ] 错误处理完善（有意义的错误信息）
- [ ] 性能无明显退化
- [ ] 测试覆盖关键路径
- [ ] 代码可读性（命名、注释）

**PR 反馈响应时间：**

- 核心团队：24 小时内首次响应
- 社区贡献者：48 小时内首次响应

---

## 6. 测试要求

### 6.1 单元测试覆盖率

- **最低要求**：新增代码覆盖率 ≥85%
- **目标**：核心模块覆盖率 ≥90%（`storage/`, `parsers/`, `cli/`）

### 6.2 测试文件组织

```
tests/
├── unit/
│   ├── storage/
│   │   ├── test_object_store.py
│   │   ├── test_metadata_db.py
│   │   └── test_commit_builder.py
│   ├── parsers/
│   │   ├── test_incar_parser.py
│   │   ├── test_poscar_parser.py
│   │   └── test_outcar_extractor.py
│   └── cli/
│       ├── test_init.py
│       ├── test_add.py
│       └── test_commit.py
├── integration/
│   ├── test_init_workflow.py
│   ├── test_commit_workflow.py
│   └── test_reproduce_workflow.py
├── fixtures/
│   └── ...
└── conftest.py
```

### 6.3 测试命名规则

```python
# 测试函数命名：test_<被测函数>_<场景>_<预期结果>
def test_write_blob_deduplication_skips_existing():
    """测试 write_blob 在内容已存在时跳过写入"""
    ...

def test_parse_incar_with_comments_ignores_them():
    """测试 parse_incar 在遇到注释时正确忽略"""
    ...
```

### 6.4 运行测试

```bash
# 运行所有测试
pytest

# 运行特定模块
pytest tests/unit/parsers/

# 带覆盖率报告
pytest --cov=chemvcs --cov-report=html

# 快速模式（跳过慢速测试）
pytest -m "not slow"

# 详细模式（失败时显示完整 traceback）
pytest -vv
```

---

## 7. 文档规范

### 7.1 文档类型

| 类型 | 位置 | 受众 |
|------|------|------|
| 用户文档 | `docs/` | 终端用户 |
| API 文档 | Docstring（自动生成） | 开发者 |
| 架构文档 | `docs/TECH_ARCHITECTURE.md` | 贡献者 |
| 开发规范 | `CONTRIBUTING.md`（本文） | 贡献者 |

### 7.2 更新文档的时机

**必须同步更新文档：**

- 新增 CLI 命令 → 更新 `CLI_SPEC.md`
- 修改退出码 → 更新 `CLI_SPEC.md` 的退出码表
- 新增存储字段 → 更新 `TECH_ARCHITECTURE.md` 的表结构
- Breaking change → 在 `CHANGELOG.md` 中特别标注

---

## 8. 发布流程

### 8.1 版本号规范（Semantic Versioning）

```
v<MAJOR>.<MINOR>.<PATCH>

v0.1.0 → 初次 MVP 发布
v0.1.1 → Bug 修复
v0.2.0 → 新增 HPC 集成功能
v1.0.0 → 稳定版（生产就绪）
```

### 8.2 Release Checklist

```markdown
- [ ] 所有 CI 通过（dev 分支）
- [ ] 版本号更新（`chemvcs/__init__.py`, `pyproject.toml`）
- [ ] CHANGELOG.md 更新
- [ ] 创建 release/v<version> 分支
- [ ] 在 HPC 环境手动测试（至少 2 种文件系统）
- [ ] 合并 release 分支到 main
- [ ] 在 main 上打 tag：`git tag -a v0.1.0 -m "Release v0.1.0"`
- [ ] 推送 tag：`git push upstream v0.1.0`
- [ ] 构建并上传到 PyPI：`python -m build && twine upload dist/*`
- [ ] 在 GitHub 上创建 Release（附带 changelog）
- [ ] 公告（邮件列表 / 论坛 / Twitter）
```

---

## 9. 问题追踪

### 9.1 Issue 模板

#### Bug Report

```markdown
**Environment:**
- ChemVCS version: 
- Python version: 
- OS: 
- File system (if HPC): 

**Steps to Reproduce:**
1. 
2. 

**Expected Behavior:**

**Actual Behavior:**

**Error Message/Traceback:**
```

#### Feature Request

```markdown
**Problem Statement:**
<!-- 描述你想解决的问题 -->

**Proposed Solution:**
<!-- 你建议的解决方案 -->

**Alternatives Considered:**
<!-- 其他方案 -->

**Additional Context:**
<!-- 任何其他信息 -->
```

### 9.2 Issue 标签

| 标签 | 说明 |
|------|------|
| `bug` | Bug 报告 |
| `enhancement` | 新功能请求 |
| `documentation` | 文档改进 |
| `good first issue` | 适合新贡献者 |
| `help wanted` | 需要社区帮助 |
| `priority: high` | 高优先级 |
| `wontfix` | 不会修复 |

---

## 10. 社区行为准则

### 10.1 核心原则

- **友好与尊重**：欢迎所有技能水平的贡献者
- **建设性反馈**：批评代码而非人；提出问题时附带解决方案
- **耐心**：社区成员可能来自不同时区和背景
- **包容性**：避免使用排他性语言（如"显而易见"、"任何人都知道"）

### 10.2 不当行为

以下行为不被接受：

- 人身攻击、侮辱、贬低性评论
- 未经许可发布他人私人信息
- 骚扰行为（公开或私下）
- 其他违反职业道德的行为

**报告机制：** 发送邮件至 <maintainer-email>（保密处理）

---

## 11. 常见问题（FAQ for Contributors）

### Q1: 我是新手，从哪里开始？

**A:** 查找标记为 `good first issue` 的 Issue，这些通常是独立且范围明确的任务。

### Q2: 我需要在真实 HPC 环境测试吗？

**A:** 对于核心功能（如 blob 存储、commit），本地测试 + 单元测试即可。HPC 特定功能（如 SLURM 集成）可在 PR 中说明"未在 HPC 测试"，由维护者安排。

### Q3: 我的 PR 被要求修改，但我不同意怎么办？

**A:** 在 PR 评论中礼貌地说明你的理由。如果仍有分歧，可标记 `@chemvcs/maintainers` 请求第三方意见。

### Q4: 如何设置远程调试（在 HPC 上）？

**A:** 使用 `chemvcs --debug` 启用详细日志。或在代码中插入 `import pdb; pdb.set_trace()` 断点（需交互式 shell）。

### Q5: 测试失败但本地无法复现怎么办？

**A:** 检查是否与操作系统相关（Linux vs Windows 路径）。可在 PR 中标注 `needs-investigation`，CI 日志通常有详细信息。

---

## 12. 联系方式

- **GitHub Issues**: <https://github.com/lichman0405/chemvcs/issues>
- **Discussion Forum**: <https://github.com/lichman0405/chemvcs/discussions>
- **Email**: shadow.li981@gmail.com（技术问题优先使用 GitHub）

---

**感谢你的贡献！** 🎉

所有贡献者将被列入 `CONTRIBUTORS.md` 文件。我们致力于构建一个友好且高效的计算材料开源社区。
