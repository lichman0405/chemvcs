# ChemVCS — Claude Code 工作说明

## 项目简介

ChemVCS 是一个面向计算材料科学的 Git 式版本控制工具，支持 VASP、LAMMPS、ORCA 三种计算软件的语义化 diff 和快照管理。

- **版本**：0.1.0 (Alpha)
- **Python**：3.10+
- **入口点**：`chemvcs = "chemvcs.cli.main:app"`

## 目录结构

```
src/chemvcs/
├── cli/main.py          # 所有 CLI 命令（Typer）
├── core/staging.py      # StagingManager（暂存区管理）
├── storage/
│   ├── object_store.py  # SHA-256 内容寻址对象存储
│   ├── commit_builder.py# 创建 commit JSON 文件
│   └── metadata_db.py   # SQLite 索引（可从 commits/ 重建，非源码真值）
├── parsers/
│   ├── base_parser.py   # BaseParser 抽象类，DiffEntry
│   ├── diff_engine.py   # 按文件名路由到对应 parser
│   ├── incar_parser.py, kpoints_parser.py, outcar_parser.py
│   ├── lammps_input_parser.py, lammps_data_parser.py, lammps_log_parser.py
│   └── orca_input_parser.py, orca_output_parser.py
├── plugins/
│   ├── base.py          # Plugin, ValidatorPlugin 基类
│   └── manager.py       # PluginManager（发现、加载、执行）
└── constants.py         # 路径常量、文件类型枚举
```

`.chemvcs/` 仓库内部结构：

```
.chemvcs/
├── objects/<hash[:2]>/<hash[2:]>  # blob 对象（SHA-256 内容寻址，≥200MB 自动 gzip）
├── commits/<commit_hash>          # commit JSON 文件（含 parent 哈希，链表结构）
├── metadata.db                    # SQLite 索引（派生，可从 commits/ 重建，不是源码真值）
├── index                          # 暂存区（JSON）
└── HEAD                           # 当前 commit 哈希（纯文本）
```

## 常用命令

```bash
# 安装（开发模式）
pip install -e ".[dev]"

# 运行全部测试（最低覆盖率 87%）
pytest

# 运行 linter
ruff check src/ tests/

# 运行类型检查
mypy src/

# 格式化代码
ruff format src/ tests/

# 查看 HTML 覆盖率报告
pytest --cov=chemvcs --cov-report=html
```

## 编码规范

- 行长度：100 字符（ruff 配置）
- 类型注解：全面（mypy strict 模式，`disallow_untyped_defs=true`）
- 格式：ruff format（Black 兼容）
- 导入排序：isort（通过 ruff I 规则）
- 错误处理：使用项目已有的异常类（`BlobNotFoundError`、`CommitBuilderError` 等）
- 测试：单元测试放 `tests/unit/`，集成测试放 `tests/integration/`
- 不要在未修改的代码上添加 docstring 或注释
- 不要过度工程化，只做被要求或明确必要的改动

---

## 任务一：Bug 审查与修复（优先执行，在任何新功能之前完成）

请全面审查现有代码库中的 bug，重点关注以下模块：

1. **`storage/object_store.py`** — 原子写入（tmp → rename）是否在 Windows 路径上正确？压缩/解压逻辑是否健壮？哈希校验是否覆盖全部读写路径？
2. **`storage/commit_builder.py`** — parent 哈希继承逻辑是否正确？文件快照是否完整继承父 commit 的文件列表？
3. **`storage/metadata_db.py`** — SQLite WAL 模式切换逻辑（NFS 回退）是否安全？事务边界是否正确？重建逻辑是否幂等？
4. **`core/staging.py`** — `.chemvcsignore` 匹配逻辑、index JSON 序列化/反序列化是否存在边界情况？
5. **`cli/main.py`** — 所有命令的错误处理是否正确使用 `constants.py` 中定义的 exit code？Rich 输出是否有未捕获的异常路径？
6. **`parsers/`** — 各 parser 的 regex fallback 逻辑、文件类型检测边界情况、`DiffEntry` 构建是否有类型不一致？
7. **`plugins/manager.py`** — 插件加载失败时的异常处理是否安全（不应崩溃主流程）？

修复后验收：
```bash
pytest          # 全部通过，覆盖率 ≥87%
ruff check src/ # 无报错
mypy src/       # 无报错
```

---

## 任务二：Phase 3 — 远端协作（GitHub 风格 push/pull）

### 设计原则

chemvcs 的存储天然支持远端同步：
- `objects/` 是纯内容寻址（SHA-256），天然幂等，rsync 传输安全
- `commits/` 是 JSON 文件，每个 commit 含 `parent` 哈希（链表结构），可遍历历史
- `metadata.db` 是派生索引，**不需要同步**，在目标机器上从 `commits/` 重建
- `index`（暂存区）是本地工作状态，**不同步**

**传输方式**：rsync over SSH（subprocess 调用系统 `rsync`/`ssh`，不引入 `paramiko` 依赖）

**冲突策略**：Fast-forward only（与 GitHub 默认行为完全一致）
- push 前检查：远端 HEAD 必须是本地 HEAD 的祖先
- 若非快进 → 拒绝 push，提示用户先 pull

### 新增配置文件

**`.chemvcs/remotes.toml`**（由 `chemvcs remote add` 写入，用户不应手动编辑）：

```toml
[origin]
url = "user@host:/path/to/remote/repo"
```

**Python 3.10 兼容的 TOML 读取**（在所有需要读取 TOML 的模块中使用）：

```python
try:
    import tomllib          # Python 3.11+（标准库）
except ImportError:
    import tomli as tomllib  # type: ignore[no-redef]
```

需在 `pyproject.toml` 的 `dependencies` 中添加：
```toml
"tomli>=2.0.0; python_version < '3.11'",
```

### 新增代码模块

```
src/chemvcs/
└── remote/
    ├── __init__.py
    ├── config.py    # RemoteConfig：读写 .chemvcs/remotes.toml
    └── manager.py   # RemoteManager：SSH/rsync 封装
```

#### `remote/config.py`

提供以下函数（模块级，不需要实例化类）：

- `list_remotes(chemvcs_dir: Path) -> Dict[str, str]` — 返回 `{name: url}`
- `add_remote(chemvcs_dir: Path, name: str, url: str) -> None` — 若已存在抛 `RemoteAlreadyExistsError`
- `remove_remote(chemvcs_dir: Path, name: str) -> None` — 若不存在抛 `RemoteNotFoundError`
- `get_remote_url(chemvcs_dir: Path, name: str) -> str` — 若不存在抛 `RemoteNotFoundError`

自定义异常：`RemoteNotFoundError`、`RemoteAlreadyExistsError`（继承 `ValueError`）

#### `remote/manager.py`

```python
class RemoteManager:
    def __init__(self, url: str) -> None:
        # 解析 "user@host:path" 格式
        # 存储 self.user_host（"user@host"）和 self.remote_path

    def check_connection(self) -> bool:
        # ssh user@host 'echo ok'，返回是否成功

    def get_remote_head(self) -> Optional[str]:
        # SSH 读取远端 .chemvcs/HEAD，不存在返回 None

    def init_remote_if_needed(self) -> None:
        # SSH 在远端创建 .chemvcs/objects/ 和 .chemvcs/commits/ 目录（首次 push）

    def push_objects(self, local_chemvcs_dir: Path) -> None:
        # rsync local objects/ → remote（--ignore-existing，不覆盖已有 blob）

    def push_commits(self, local_chemvcs_dir: Path) -> None:
        # rsync local commits/ → remote（--ignore-existing）

    def fetch_objects(self, local_chemvcs_dir: Path) -> None:
        # rsync remote objects/ → local（--ignore-existing）

    def fetch_commits(self, local_chemvcs_dir: Path) -> None:
        # rsync remote commits/ → local（--ignore-existing）

    def update_remote_head(self, commit_hash: str) -> None:
        # SSH 写入远端 .chemvcs/HEAD（echo hash > path）
```

### 修改现有文件

**`src/chemvcs/constants.py`** — 添加：
```python
REMOTES_FILE = "remotes.toml"
```

**`src/chemvcs/cli/main.py`** — 添加命令：

```
remote_app = typer.Typer(name="remote", help="Manage remote repositories")
app.add_typer(remote_app)

chemvcs remote add <name> <url>    # 注册 remote
chemvcs remote list                 # 列出所有 remote，格式：name  url
chemvcs remote remove <name>        # 删除 remote

@app.command()
def push(remote: str) -> None:
    # 见下方 push 逻辑

@app.command()
def pull(remote: str) -> None:
    # 见下方 pull 逻辑
```

### push 详细逻辑

```
1. 查找 .chemvcs/（向上遍历，找不到报错）
2. 读取本地 HEAD（若无 commit 则报错）
3. 从 remotes.toml 获取 remote URL
4. RemoteManager(url).check_connection()，失败则报错
5. get_remote_head() → remote_head
6. 快进检查：
   a. remote_head is None → 允许（首次 push）
   b. remote_head == local_head → 提示"Already up to date"，退出
   c. 遍历本地 commit 链，若 remote_head 在祖先中 → 允许（快进）
   d. 否则 → 报错："Remote has diverged. Run: chemvcs pull <remote>"
7. init_remote_if_needed()
8. push_objects(chemvcs_dir)
9. push_commits(chemvcs_dir)
10. update_remote_head(local_head)
11. 打印：Pushed to <remote> (<url>)
```

### pull 详细逻辑

```
1. 查找 .chemvcs/
2. 读取本地 HEAD（可为 None，即空仓库）
3. 从 remotes.toml 获取 remote URL
4. RemoteManager(url).check_connection()
5. get_remote_head() → remote_head，若 None 则报错"Remote is empty"
6. local_head == remote_head → 提示"Already up to date"，退出
7. fetch_objects(chemvcs_dir)
8. fetch_commits(chemvcs_dir)
9. 快进检查（从远端视角）：
   a. local_head is None → 直接更新 HEAD 到 remote_head
   b. 遍历远端 commit 链（已拉到本地），若 local_head 在祖先中 → 快进，更新 HEAD
   c. 否则 → 报错："Local has diverged. Push your commits first."
10. 更新本地 HEAD
11. 重建 metadata.db（调用现有的 MetadataDB 重建机制）
12. 打印：Pulled N commits from <remote>
```

### 测试要求

- `tests/unit/remote/test_config.py` — RemoteConfig 增删查，包括错误场景
- `tests/unit/remote/test_manager.py` — 用 `unittest.mock.patch("subprocess.run")` mock，测试 SSH/rsync 命令参数构造正确性
- `tests/integration/test_push_pull.py` — 用两个临时目录模拟 local/remote 仓库，mock SSH/rsync 为本地文件复制，测试完整 push → pull 流程
- 快进检查逻辑单独覆盖：空远端、正常快进、已是最新、分叉拒绝四种场景

### 验收示例

```bash
# 注册远端
chemvcs remote add origin user@server:/data/repos/my_calc

# 推送
chemvcs push origin
# Pushed to origin (user@server:/data/repos/my_calc)

# 查看 remote 列表
chemvcs remote list
# origin    user@server:/data/repos/my_calc

# 另一台机器拉取（先 init）
chemvcs init
chemvcs pull origin
# Pulled 5 commits from origin

# 非快进时的拒绝
chemvcs push origin
# Error: Remote has diverged. Run: chemvcs pull origin
```

---

## 任务执行顺序

1. **任务一**：Bug 审查 → 修复 → `pytest` 全绿
2. **任务二**：按 Phase 3 设计实现 `remote/`、`push`、`pull` 命令 → 测试全绿

