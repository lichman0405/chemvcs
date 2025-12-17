# ChemVCS 命令大全（超详细）

本文专门讲 `chemvcs`（CLI 客户端）每个命令怎么用、参数含义、常见输出与故障排查。

适用人群：运维/研发/CI；尤其适合“只知道要跑什么命令，但不熟 ChemVCS”的同学。

- 运维部署/测试/排障总指南见：[Tutorial.md](Tutorial.md)
- 命令实现与 Help 输出源代码见：[go/cmd/chemvcs/main.go](../go/cmd/chemvcs/main.go)

---

## 1. 约定与通用规则

### 1.1 仓库与工作目录

绝大多数命令默认在“当前目录”打开仓库：
- 你需要在 ChemVCS 仓库目录（包含 `.chemvcs/` 的目录）里执行命令。
- 如果当前目录不是仓库，通常会报类似 “failed to open repository” 的错误。

### 1.2 哈希/ID 的基本要求

- ChemVCS 对象与快照使用 SHA-256 哈希（通常是 64 个十六进制字符）。
- 多个命令会把 `run-hash` 或 `snapshot` 直接截取 `[:8]` 用于显示。
  - 所以参数至少要 ≥ 8 个字符；生产使用建议传完整 64 字符。

### 1.3 Exit Code（对 CI 很重要）

从 `chemvcs` 主程序逻辑可得：
- `0`：成功
- `1`：命令执行失败（Error: ...）
- `2`：用法错误（缺少子命令/参数、未知命令等）

### 1.4 获取帮助

- `chemvcs help`
- `chemvcs --help`
- `chemvcs -h`

注意：当前实现是“总览式 help”，不会像 Git 那样给每个子命令单独的 `--help` 文本；每个子命令的更细参数以本文与源码为准。

---

## 2. 环境变量（常用）

### 2.1 作者信息（commit 默认作者）

- `CHEMVCS_AUTHOR_NAME`：默认作者名
- `CHEMVCS_AUTHOR_EMAIL`：默认作者邮箱

当两者都设置时，作者会拼成：`<name> <email>` 形态（详见 [go/cmd/chemvcs/main.go](../go/cmd/chemvcs/main.go) 的 `getAuthor()`）。

### 2.2 远端鉴权 Token（push/pull/fetch 与远端 HPC 网关通用）

令牌选择优先级（详见 `remoteTokenFor()`）：

1) `CHEMVCS_REMOTE_TOKEN_<REMOTE_NAME>`
2) `CHEMVCS_REMOTE_TOKEN`

其中 `<REMOTE_NAME>` 会被归一化：
- 转为大写
- 非字母数字字符会变成 `_`

示例：
- remote 名为 `origin` → `CHEMVCS_REMOTE_TOKEN_ORIGIN`
- remote 名为 `slurm-prod` → `CHEMVCS_REMOTE_TOKEN_SLURM_PROD`

---

## 3. 命令总览（与 help 一致）

```text
init, commit, log, branch, checkout, status, merge,
remote add,
push, pull, fetch,
inspect-object, list-objects,
submit, jobs, retrieve, cancel, watch,
pack, gc, fsck,
version, help
```

---

## 4. 基础版本控制命令

### 4.1 init —— 初始化仓库

用途：在指定路径（默认当前目录）创建 `.chemvcs/`。

用法：
- `chemvcs init [path]`

参数：
- `path`：可选。仓库目录路径。

示例：

```bash
mkdir demo
cd demo
chemvcs init .
```

常见问题：
- 报 “directory does not exist”：你给的 `path` 不存在。

---

### 4.2 status —— 查看工作区变更

用途：把当前工作区与 HEAD 快照做对比，列出 Added/Modified/Deleted。

用法：
- `chemvcs status`

输出说明：
- 如果没有任何提交：打印 `No commits yet`
- 如果工作区干净：打印 `No changes (working directory clean)`
- 变更行会带颜色（A/M/D）。
  - 在不支持 ANSI 颜色的终端里可能看到转义序列。

重要注意：
- status 会扫描工作区并把“扫描得到的根对象”临时写入对象库（用于比较）。这是预期行为。

---

### 4.3 commit —— 创建快照

用途：扫描工作区并创建新快照（snapshot）。

用法：
- `chemvcs commit -m <message> [--author="Name <email>"]`

参数：
- `-m`：必填。提交信息。
- `--author`：可选。覆盖环境变量作者。

示例：

```bash
echo hello > a.txt
chemvcs commit -m "add a.txt"
```

常见问题：
- 报 “commit message required (use -m)”：忘了 `-m`。

---

### 4.4 log —— 查看历史

用途：按时间倒序显示最近快照。

用法：
- `chemvcs log [-n <count>] [--oneline]`

参数：
- `-n`：默认 20。
- `--oneline`：简洁格式（只显示短 hash + message）。

示例：

```bash
chemvcs log
chemvcs log -n 5 --oneline
```

---

### 4.5 branch —— 分支管理（当前只支持 list / create）

用途：
- 无参数：列出分支
- 有参数：创建新分支

用法：
- `chemvcs branch`
- `chemvcs branch <name>`

限制：
- 分支名不能包含空格
- 分支名不能包含 `/`

示例：

```bash
chemvcs branch
chemvcs branch exp
```

输出：
- 列表中当前分支带 `*`。

---

### 4.6 checkout —— 切换分支或快照

用途：
- 切到分支：更新 HEAD 指向分支
- 切到快照：进入 detached HEAD

用法：
- `chemvcs checkout <branch|snapshot>`

重要注意（非常关键）：
- checkout 会“恢复工作区”到目标快照对应的内容。
- 如果你工作区有未提交修改，可能被覆盖/丢失。
  - 建议先 `chemvcs status`，再决定是否 `chemvcs commit`。

示例：

```bash
chemvcs checkout main
chemvcs checkout a1b2c3d4e5f6...
```

---

### 4.7 merge —— 合并分支

用途：把指定分支合并到当前分支。

用法：
- `chemvcs merge <branch>`

限制：
- 如果 HEAD 是 detached，merge 会直接报错：`cannot merge: HEAD is detached`。

冲突处理：
- 发生冲突时命令会列出冲突文件并返回错误：`merge aborted due to conflicts`。

成功后：
- 会更新 HEAD，并恢复工作区到新 HEAD。

---

## 5. 远端（remote/push/pull/fetch）

ChemVCS 的远端同步是通过 `chemvcs-server` 的 HTTP API 完成的。

参考文档：
- 远端协议：[docs/05-remote-protocol-and-server.md](../docs/05-remote-protocol-and-server.md)
- 运维部署：[Tutorial.md](Tutorial.md)

### 5.1 remote add —— 添加远端

当前只支持一个子命令：`add`。

用法：
- `chemvcs remote add <name> <url>`

建议 URL 形态（repo-scoped）：
- `http://<host>:<port>/chemvcs/v1/repos/<owner>/<repo>`

示例：

```bash
chemvcs remote add origin http://127.0.0.1:8080/chemvcs/v1/repos/owner/repo
```

常见问题：
- `unknown remote subcommand`：目前没有 `remote rm` / `remote list`。

---

### 5.2 push —— 推送分支

用法：
- `chemvcs push <remote> <branch>`

行为：
- 会把本地 `refs/heads/<branch>` 推送到远端同名引用。
- 当前实现 `Force=false`，因此默认只允许快进（fast-forward）。

鉴权：
- 通过 `CHEMVCS_REMOTE_TOKEN` 或 `CHEMVCS_REMOTE_TOKEN_<REMOTE>` 提供 bearer token。

常见问题：
- 401：token 未提供/错误
- 403：token 没被授权访问该 repo
- 404：服务端没有该仓库目录或缺 `.chemvcs`
- not a fast-forward：远端领先且发生分叉

---

### 5.3 pull —— 拉取分支并更新工作区

用法：
- `chemvcs pull <remote> <branch>`

重要注意：
- pull 完成后会直接把工作区恢复到新 HEAD。
- 若本地有未提交修改，可能被覆盖。

建议流程：
- `chemvcs status` 确认工作区干净
- 再 `chemvcs pull ...`

---

### 5.4 fetch —— 仅拉取对象，不改工作区

用法：
- `chemvcs fetch <remote> <branch>`

输出：
- 成功时会打印获取到的 snapshot 短 hash。

---

## 6. 对象查看与调试

### 6.1 inspect-object —— 查看对象内容

用法：
- `chemvcs inspect-object <hash> [--format=text|json]`

参数：
- `--format`：默认 `text`。

示例：

```bash
chemvcs inspect-object <hash>
chemvcs inspect-object <hash> --format=json
```

---

### 6.2 list-objects —— 列出对象

用法：
- `chemvcs list-objects [--type=<type>] [--format=text|json]`

参数：
- `--type`：按对象类型过滤（具体类型取决于对象模型）。
- `--format`：默认 `text`。

示例：

```bash
chemvcs list-objects
chemvcs list-objects --type=blob
chemvcs list-objects --format=json
```

---

## 7. 维护命令（pack/gc/fsck）

参考：对象与 pack 机制说明见 [docs/03-object-model-and-storage.md](../docs/03-object-model-and-storage.md)

### 7.1 pack —— 打包 loose objects

用途：把 loose objects 打包成 packfile，减少文件数量并提升部分场景 IO。

用法：
- `chemvcs pack [--all] [--keep-loose] [-v]`

参数：
- `--all`：打包所有 loose（包括不可达对象）。
- `--keep-loose`：打包后保留 loose（更安全但占空间）。
- `-v`：打印详细过程。

建议：
- 第一次跑建议加 `--keep-loose`，确认无误后再去掉。

---

### 7.2 gc —— 垃圾回收

用途：清理不可达对象。

用法：
- `chemvcs gc [--prune=<age>] [--dry-run] [-v]`

参数：
- `--prune`：默认 `2w`。支持单位：`s/m/h/d/w`，或 `now`。
- `--dry-run`：只输出将删除什么，不实际删除。
- `-v`：详细输出。

示例：

```bash
chemvcs gc --dry-run
chemvcs gc --prune=168h
chemvcs gc --prune=now
```

安全建议：
- 先 `chemvcs fsck --full` 再 `gc`。
- 遇到疑似损坏先停 `gc`。

---

### 7.3 fsck —— 完整性检查

用途：校验仓库一致性。

用法：
- `chemvcs fsck [--full] [-v]`

参数：
- `--full`：包含 packfile 校验（更慢但更完整）。
- `-v`：详细输出。

退出码：
- 若未开启 `-v` 且存在 error，命令会直接 `os.Exit(1)`。

---

## 8. HPC 命令（本地模式与远端网关模式）

HPC 总指南：
- 用户视角：[docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md)
- 远端网关实现细节：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)

### 8.1 统一概念：run-hash / job-id

- `run-hash`：你某次计算/工作流运行对应的标识（通常是一次 run 对象的 hash）。
- `job-id`：调度系统返回的作业 ID（数字或字符串）。

多个命令支持用 `run-hash` 或 `job-id` 作为标识。

---

### 8.2 submit —— 提交作业

用法：
- 本地模式：`chemvcs submit <run-hash> <script> [--capture-env=true|false]`
- 远端网关：`chemvcs submit --remote=<name> <run-hash> <script> [--capture-env=true|false]`

参数：
- `--capture-env`：默认 true。
- `--remote`：通过 `chemvcs-server` 远端网关提交。

注意：
- `script` 必须存在；远端模式会读取脚本内容并通过 HTTP 发送。

---

### 8.3 jobs —— 列出作业

用法：
- `chemvcs jobs [--status=<status>] [-v] [--remote=<name>] [<run-hash|job-id>]`

参数：
- `--status`：例如 `RUNNING`、`COMPLETED`。
- `-v`：显示更多字段（Submitted/Updated/Workdir/ExitCode/Reason 等）。
- `--remote`：从远端网关拉取。

过滤规则：
- 传 `job-id`（纯数字）会按 jobId 精确匹配。
- 传 `run-hash` 前缀会按 runHash 前缀匹配；若匹配多个会报 ambiguous。

---

### 8.4 watch —— 轮询直到完成

用法：
- `chemvcs watch <run-hash|job-id> [--interval=<sec>] [--timeout=<sec>] [--remote=<name>]`

参数：
- `--interval`：默认 30 秒。
- `--timeout`：默认 0（无限）。

远端模式行为：
- 状态变化时打印 `Status: ...`。
- 终止状态（COMPLETED/FAILED/CANCELLED）直接返回。

---

### 8.5 cancel —— 取消作业

用法：
- `chemvcs cancel <run-hash|job-id> [--remote=<name>]`

---

### 8.6 retrieve —— 取回结果（可选自动 commit）

用法：
- 本地模式：
  - `chemvcs retrieve <run-hash> [--patterns=<p1,p2,...>] [--dest=<path>] [--commit] [--commit-message=<msg>]`
- 远端网关：
  - `chemvcs retrieve --remote=<name> <run-hash> ...`

参数：
- `--patterns`：逗号分隔，例如 `*.out,*.log`。
- `--dest`：默认 `.`。
- `--commit`：取回后创建一次新的 snapshot。
- `--commit-message`：可选；不填则生成 `Retrieve results for run <run-hash>`。

安全限制（避免把文件写到仓库外）：
- 当 `--commit` 开启时，`--dest` 必须在仓库目录内；否则会报错。

远端模式说明：
- 远端会返回 zip；客户端解压到 `--dest`。

---

## 9. 其它命令

### 9.1 version

用法：
- `chemvcs version`

### 9.2 help

用法：
- `chemvcs help`
- `chemvcs --help`
- `chemvcs -h`

---

## 10. 运维与排障速查（命令级）

- `Error: failed to open repository`：当前目录不是仓库；确认 `.chemvcs/` 存在。
- `Unknown command: ...`：命令拼写错误；先跑 `chemvcs help`。
- `usage: ...`：参数数量不对；按本文用法补齐。
- 远端 401/403/404：按 [Tutorial.md](Tutorial.md) 的鉴权与 repo-root 检查。
- HPC 远端模式失败：优先在服务端用服务账号手动运行 `sbatch/squeue/sacct/scancel` 验证环境。
