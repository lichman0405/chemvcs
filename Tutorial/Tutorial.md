# ChemVCS 运维部署 / 测试 / 使用 / 排障（超详细教程）

面向对象：初级运维/初级研发（能看懂命令行即可）。

目标：让你在没有“项目背景知识”的情况下，也能把 ChemVCS 从源码构建出来，部署 `chemvcs-server`，完成本地与远端的基本使用、HPC（本地/远端网关）使用，并能按清单进行故障排查。

本文以仓库现状为准，并引用项目内已有文档作为“进一步参考”。

- 项目概览与功能面：[README.md](../README.md)
- 功能状态（最权威的“现在能做什么”）：[FEATURE_STATUS.md](../FEATURE_STATUS.md)
- 远端协议/服务设计（API 语义）：[docs/05-remote-protocol-and-server.md](../docs/05-remote-protocol-and-server.md)
- 对象模型与本地存储（.chemvcs 布局）：[docs/03-object-model-and-storage.md](../docs/03-object-model-and-storage.md)
- HPC 使用指南（CLI/远端网关）：[docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md)
- 远端 HPC 网关实现说明（最贴近实现）：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)
- HPC 示例工作流（可当验收用例）：[examples/hpc-workflow/README.md](../examples/hpc-workflow/README.md)
- ChemVCS 与 Git 差异理解：[docs/10-chemvcs-vs-git.md](../docs/10-chemvcs-vs-git.md)
- ChemVCS 命令大全（超详细）：[ChemVCS-Commands.md](ChemVCS-Commands.md)

---

## 0. 你需要先知道的几件事（非常重要）

1) ChemVCS 有两个 Go 二进制：
- `chemvcs`：CLI 客户端（本地操作仓库、push/pull/fetch、维护命令、HPC 命令等）。
- `chemvcs-server`：HTTP 服务端（托管远端仓库 + 远端 HPC 网关）。

2) 远端仓库的 URL 形态与目录结构（别踩坑）：
- 服务端 API 前缀固定是 `/chemvcs/v1/`。
- 一个远端仓库的“repoId”是 `owner/repo`。
- 服务端的仓库目录必须存在且包含 `.chemvcs`，路径是：`<repo-root>/owner/repo/.chemvcs`。
  - 服务端不会自动“创建远端仓库”。如果目录里没有 `.chemvcs`，会报 `repository not found`。
    - 这点在实现里是硬约束（见服务端在 [go/internal/server/server.go](../go/internal/server/server.go) 打开的逻辑）。

3) 远端同步（push/pull/fetch）与远端 HPC 网关（no SSH）使用同一个 base URL：
- 推荐直接把 remote URL 配成 repo-scoped：
  - `http://<host>:<port>/chemvcs/v1/repos/<owner>/<repo>`
- 客户端会自动从这个 URL 推断 repoId，并把 base URL 归一化成 API 根（见远端客户端构造逻辑，参考 [go/internal/remote/client.go](../go/internal/remote/client.go)）。

4) 版本命令：
- `chemvcs version`（不是 `chemvcs --version`）
- `chemvcs-server -version`

5) HPC 分两种模式：
- 本地模式：`chemvcs` 在“能直接运行调度器命令”的机器上工作（例如 SLURM login node）。这条路径需要 Python 包 `chemvcs_py`。
- 远端网关模式（no SSH）：你的笔记本只发 HTTP 给 `chemvcs-server`，服务端在 SLURM 主机上执行 `sbatch/squeue/sacct/scancel`（MVP 目前为 SLURM）。这条路径不需要在服务端运行 Python（见 [docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md) 的关键约束）。

---

## 1. 部署拓扑（先选你要哪种部署）

### 1.1 最简（单机开发/自测）
- 同一台机器：构建 `chemvcs`，本地用 `.chemvcs` 仓库。
- 适用：功能验证、跑单元测试。

### 1.2 远端仓库（团队协作）
- 一台服务器运行 `chemvcs-server`（存储多个仓库）。
- 多个客户端通过 `remote add` + `push/pull/fetch` 同步。
- 适用：实验室/组内共享仓库。

### 1.3 HPC 本地模式
- 在 SLURM/PBS/LSF 可用的机器（一般是登录节点）安装：
  - `chemvcs` + Python 包 `chemvcs_py`
- 用 `chemvcs submit/jobs/watch/retrieve` 直接操作作业。
- 参考：[docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md)

### 1.4 HPC 远端网关模式（推荐给运维场景）
- 在 SLURM 主机上部署 `chemvcs-server`，并确保服务账号能执行调度命令。
- 客户端（笔记本）只需 `chemvcs`，通过 `--remote=<name>` 使用 HTTP 网关。
- 参考：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)

---

## 2. 前置条件与依赖清单（按角色）

### 2.1 构建机（编译二进制）
- Go：建议 1.19+（见 [README.md](../README.md)）
- Git：用于拉取源码（可选但强烈建议）

### 2.2 客户端（运行 chemvcs）
- Windows / Linux 均可
- 如果只做 VCS 与远端同步：不需要 Python
- 如果要“本地模式”HPC：
  - Python 3.8+（见 [docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md)）
  - 安装 `chemvcs_py`
  - 对应调度器命令（如 SLURM：`sbatch/squeue/sacct/scancel`）必须在 PATH 中

### 2.3 服务端（运行 chemvcs-server）
- 建议 Linux（更常见的运维场景），Windows 也可运行
- 存储：用于 repo-root（包含多个仓库目录）
- 网络：对外暴露 TCP 端口（默认 8080）
- 鉴权：建议开启 bearer token（见 [docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md) 的部署说明）
- 如果开启远端 HPC 网关（SLURM）：
  - `chemvcs-server` 运行账号必须能执行 SLURM 命令

---

## 3. 从源码构建与安装（二进制）

> 仓库结构与设计背景可先看 [README.md](../README.md) 与 [docs/02-architecture-overview.md](../docs/02-architecture-overview.md)。

### 3.1 构建 chemvcs（CLI）

在仓库根目录：

1) 进入 Go 模块目录：

```bash
cd go
```

2) 构建：

```bash
go build -o chemvcs ./cmd/chemvcs
```

- Windows 下会生成 `chemvcs.exe`

3) 验证：

```bash
./chemvcs version
```

### 3.2 构建 chemvcs-server（服务端）

在 `go/` 目录：

```bash
go build -o chemvcs-server ./cmd/chemvcs-server
```

验证：

```bash
./chemvcs-server -version
```

### 3.3 安装到 PATH（建议）

- Linux：把二进制放到 `/usr/local/bin/` 或者在 `~/.local/bin/` 并加入 PATH。
- Windows：把 `chemvcs.exe` 所在目录加入“系统环境变量 PATH”。

---

## 4. Python 层安装（仅在需要本地 HPC/领域对象时）

Python 层说明与安装参考：[python/README.md](../python/README.md)

在仓库根目录：

```bash
cd python
pip install -e .
```

验证模块导入（示例）：

```bash
python -c "from chemvcs_py.hpc import SlurmAdapter; print('OK')"
```

注意：
- 如果你只部署远端 `chemvcs-server` 作为 HPC 网关（no SSH），服务端不需要 Python。
- 如果你要在登录节点直接运行 `chemvcs submit`（本地模式），需要 Python。

---

## 5. 部署 chemvcs-server（远端仓库 + 可选远端 HPC 网关）

### 5.1 目录规划（repo-root）

服务端需要一个“仓库根目录”用于存放多个仓库。假设你选：
- Linux：`/srv/chemvcs/repos`
- Windows：`D:\chemvcs\repos`

远端仓库目录结构必须是：

```text
<repo-root>/
  owner1/
    repoA/
      .chemvcs/
  owner2/
    repoB/
      .chemvcs/
```

`.chemvcs` 的布局参考：[docs/03-object-model-and-storage.md](../docs/03-object-model-and-storage.md)

### 5.2 创建远端仓库（必须做，否则 404）

> 关键点：服务端不会自动创建仓库，必须提前初始化。

以 Linux 为例（Windows 类似，改路径即可）：

1) 创建目录：

```bash
mkdir -p /srv/chemvcs/repos/owner/repo
cd /srv/chemvcs/repos/owner/repo
```

2) 初始化：

```bash
/path/to/chemvcs init .
```

完成后应存在：

```text
/srv/chemvcs/repos/owner/repo/.chemvcs
```

### 5.3 启动命令与参数（以实现为准）

服务端参数来自实现：[go/cmd/chemvcs-server/main.go](../go/cmd/chemvcs-server/main.go)

最简启动（不带鉴权，不建议生产）：

```bash
./chemvcs-server -port 8080 -repo-root /srv/chemvcs/repos
```

推荐启动（开启 token 鉴权 + repo scope）：

```bash
./chemvcs-server \
  -port 8080 \
  -repo-root /srv/chemvcs/repos \
  -auth-token "change-me" \
  -auth-repos "owner/repo"
```

可选：admin token（用于列出 repos，且绕过 repo scoping）：

```bash
./chemvcs-server \
  -port 8080 \
  -repo-root /srv/chemvcs/repos \
  -auth-token "change-me" \
  -auth-repos "owner/repo" \
  -admin-token "change-me-admin"
```

也可使用环境变量（更适合 systemd / Windows 服务）：
- `CHEMVCS_SERVER_AUTH_TOKEN`
- `CHEMVCS_SERVER_AUTH_REPOS`
- `CHEMVCS_SERVER_ADMIN_TOKEN`

### 5.4 Linux systemd 部署（推荐）

1) 创建服务用户（建议）：

```bash
sudo useradd --system --home /srv/chemvcs --shell /usr/sbin/nologin chemvcs
sudo mkdir -p /srv/chemvcs/repos
sudo chown -R chemvcs:chemvcs /srv/chemvcs
```

2) 写入环境文件（权限建议 0600）：

文件：`/etc/chemvcs-server.env`

```bash
CHEMVCS_SERVER_AUTH_TOKEN=change-me
CHEMVCS_SERVER_AUTH_REPOS=owner/repo
# CHEMVCS_SERVER_ADMIN_TOKEN=change-me-admin
```

3) systemd unit：

文件：`/etc/systemd/system/chemvcs-server.service`

```ini
[Unit]
Description=ChemVCS Server
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
User=chemvcs
Group=chemvcs
WorkingDirectory=/srv/chemvcs
EnvironmentFile=-/etc/chemvcs-server.env
ExecStart=/usr/local/bin/chemvcs-server -port 8080 -repo-root /srv/chemvcs/repos
Restart=on-failure
RestartSec=2

[Install]
WantedBy=multi-user.target
```

4) 启动：

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now chemvcs-server
sudo systemctl status chemvcs-server
```

5) 查看日志：

```bash
sudo journalctl -u chemvcs-server -n 200 --no-pager
```

### 5.5 Windows 部署（最简可运维方案）

不引入第三方服务管理器的前提下，推荐用“任务计划程序”以“开机自启 + 崩溃重启”的方式跑。

1) 准备目录：
- 程序：`C:\Program Files\chemvcs\chemvcs-server.exe`
- 仓库：`D:\chemvcs\repos\owner\repo\.chemvcs`

2) 配置环境变量（系统级）：
- `CHEMVCS_SERVER_AUTH_TOKEN`
- `CHEMVCS_SERVER_AUTH_REPOS`

3) 任务计划程序：
- 触发器：开机
- 操作：启动程序
  - Program/script：`C:\Program Files\chemvcs\chemvcs-server.exe`
  - Arguments：`-port 8080 -repo-root D:\chemvcs\repos`
- “如果任务失败，重新启动”：建议开启

4) 日志：
- Windows 任务计划程序可配置输出到文件；或用包装脚本重定向 stdout/stderr。

### 5.6 TLS 与反向代理（可选但强烈建议）

实现里 `chemvcs-server` 直接提供 HTTP（不带 TLS）。生产建议：
- `chemvcs-server` 仅监听内网或 localhost
- 用 Nginx/Caddy 做 TLS 终止、限流、请求体大小限制

远端 HPC retrieve 会返回 zip，建议在代理层设置合适的 `client_max_body_size`（并结合你们的输出大小）。

更多部署安全建议参考：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)

---

## 6. 客户端基本使用（本地仓库）

### 6.1 初始化仓库

```bash
mkdir myproj
cd myproj
chemvcs init .
```

### 6.2 提交与查看历史

```bash
echo "hello" > a.txt
chemvcs status
chemvcs commit -m "add a.txt"
chemvcs log
```

### 6.3 分支、切换、合并

```bash
chemvcs branch
chemvcs branch exp
chemvcs checkout exp

echo "change" >> a.txt
chemvcs commit -m "exp change"

chemvcs checkout main
chemvcs merge exp
```

### 6.4 维护命令（pack/gc/fsck）

命令列表可见 CLI help（实现中在 [go/cmd/chemvcs/main.go](../go/cmd/chemvcs/main.go)）。更详细的逐命令说明见：[ChemVCS-Commands.md](ChemVCS-Commands.md)

- `chemvcs pack`：将 loose objects 打包为 packfile（节省 inode/提升 IO 效率）
- `chemvcs gc`：清理不可达对象（支持保留期/dry-run）
- `chemvcs fsck`：校验对象与 pack 完整性

建议运维侧的常用顺序（在停机窗口或低峰期）：

```bash
chemvcs fsck --full
chemvcs pack
chemvcs gc --prune=168h
chemvcs fsck --full
```

如果 `fsck` 报错，优先停止继续 `gc`，先定位损坏来源。

---

## 7. 远端仓库：配置与同步（push/pull/fetch）

### 7.1 服务端准备（再次强调）

在 push 之前，服务端必须已存在 `owner/repo` 目录且包含 `.chemvcs`。

### 7.2 客户端配置 remote

推荐使用 repo-scoped URL：

```bash
chemvcs remote add origin http://<host>:8080/chemvcs/v1/repos/owner/repo
```

### 7.3 客户端配置 token（如果服务端开启鉴权）

两种方式任选其一：

1) 全局 token：

```bash
export CHEMVCS_REMOTE_TOKEN="change-me"
```

2) 针对某个 remote（remote 名称会被转成大写与下划线）：

```bash
export CHEMVCS_REMOTE_TOKEN_ORIGIN="change-me"
```

### 7.4 push / pull / fetch

```bash
chemvcs push origin main
chemvcs fetch origin main
chemvcs pull origin main
```

常见失败与原因：
- `401 unauthorized`：没带 token 或 token 错误
- `403 forbidden`：token 没被授权访问该 repo（服务端的 `-auth-repos` 未包含）
- `404 not_found`：远端目录不存在或没有 `.chemvcs`
- `push rejected: not a fast-forward`：远端分支领先且发生分叉（MVP 默认只允许 fast-forward）

---

## 8. HPC 本地模式（在登录节点/能跑调度命令的机器上）

详细用法以 [docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md) 为准，这里给运维友好的“最短闭环”。

### 8.1 前置检查

- `chemvcs` 可运行：`chemvcs version`
- Python 包已安装：参考 [python/README.md](../python/README.md)
- 调度器命令可用（以 SLURM 为例）：

```bash
which sbatch
which squeue
which sacct
which scancel
```

### 8.2 提交/查看/取回

在仓库中：

```bash
chemvcs submit <run-hash> vasp.slurm
chemvcs jobs
chemvcs watch <run-hash|job-id> --interval 10 --timeout 3600
chemvcs retrieve <run-hash> --dest=./results --patterns="*.out,*.log" --commit --commit-message="Import results"
```

### 8.3 常见故障

- `sbatch: command not found`：PATH 或模块环境没加载
- 作业状态一直 UNKNOWN：`squeue/sacct` 输出解析失败或命令超时（先手动运行命令看是否卡住）
- retrieve 没有文件：pattern 匹配不到；或输出不在仓库工作目录下

---

## 9. HPC 远端网关模式（no SSH，推荐运维部署形态）

这部分以实现文档为准：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)

### 9.1 服务端要求（SLURM 主机）

- `chemvcs-server` 运行账号能执行：`sbatch/squeue/sacct/scancel`
- repo-root 下的 `owner/repo` 仓库存在且包含 `.chemvcs`
- 开启鉴权（强烈建议）

启动示例：

```bash
chemvcs-server -port 8080 -repo-root /srv/chemvcs/repos -auth-token "change-me" -auth-repos "owner/repo"
```

### 9.2 客户端配置

```bash
chemvcs remote add slurm http://<host>:8080/chemvcs/v1/repos/owner/repo
export CHEMVCS_REMOTE_TOKEN_SLURM="change-me"
```

### 9.3 通过网关提交/查看/取回

```bash
chemvcs submit   --remote=slurm <run-hash> vasp.slurm
chemvcs jobs     --remote=slurm
chemvcs watch    --remote=slurm <run-hash|job-id> --interval=10
chemvcs cancel   --remote=slurm <run-hash|job-id>
chemvcs retrieve --remote=slurm <run-hash> --patterns="*.out,*.log" --dest=./results --commit
```

### 9.4 运维侧排查要点

- 认证失败：先看服务端日志里 `principal=... status=401/403`；检查 `CHEMVCS_SERVER_AUTH_*`。
- 作业提交失败：用服务账号在服务端手动运行 `sbatch` 验证权限与环境。
- retrieve 失败：服务端会做 repo 边界保护（不会打包 `.chemvcs`，也不会允许越界路径）。

---

## 10. 测试与验收（运维最小可交付清单）

### 10.1 单元测试（建议在提交/升级前跑）

Go 测试：

```bash
cd go
go test ./...
```

Python 测试：

```bash
cd python
pytest -q tests
```

### 10.2 远端服务冒烟测试

1) 服务端启动后，看监听：

```bash
# Linux
ss -lntp | grep 8080

# Windows PowerShell
Get-NetTCPConnection -LocalPort 8080
```

2) 访问 repo info（需要 token 则加 Authorization）：

```bash
curl -H "Authorization: Bearer change-me" \
  http://<host>:8080/chemvcs/v1/repos/owner/repo
```

期望返回 JSON，包含 `id` / `default_branch`。

3) 客户端 push/pull 冒烟：
- 本地建一个小仓库，commit 一个文件
- `remote add` 指向服务端
- `push origin main` 成功

---

## 11. 故障排查手册（按现象查）

### 11.1 现象：服务端启动就退出

排查步骤：
1) 检查参数：`-repo-root` 路径是否存在/是否有权限
2) 看 stderr/stdout 日志（systemd 用 `journalctl`）
3) 常见原因：端口被占用、repo-root 无权限

### 11.2 现象：客户端报 401 unauthorized

原因：缺少或错误的 Authorization。

排查步骤：
1) 确认服务端是否开启鉴权（设置了 `-auth-token` 或 `CHEMVCS_SERVER_AUTH_TOKEN`）
2) 客户端设置 `CHEMVCS_REMOTE_TOKEN` 或 `CHEMVCS_REMOTE_TOKEN_<REMOTE>`
3) 复核 token 值是否一致（注意空格/引号）

### 11.3 现象：客户端报 403 forbidden

原因：token 正确，但未授权访问该 repo。

排查步骤：
1) 服务端 `-auth-repos` / `CHEMVCS_SERVER_AUTH_REPOS` 是否包含 `owner/repo` 或 `*`
2) 注意大小写与分隔符（逗号分隔）

### 11.4 现象：客户端报 404 repository not found

原因：服务端 repo-root 下不存在该仓库或缺少 `.chemvcs`。

排查步骤：
1) 在服务端检查目录：`<repo-root>/owner/repo/.chemvcs` 是否存在
2) 如果没有：进入 `<repo-root>/owner/repo` 执行 `chemvcs init .`

### 11.5 现象：push 被拒绝 not a fast-forward

原因：远端分支领先，且本地历史不是其后代（发生分叉）。

排查建议：
1) 先 `chemvcs fetch origin main`
2) 再 `chemvcs pull origin main`（若本地也有提交，需人工处理分叉；MVP 默认不支持强制覆盖）

### 11.6 现象：fsck 报错

建议流程：
1) 立即停止 `gc`（避免误删加重损坏）
2) 记录错误对象 hash
3) 优先排查磁盘/网络存储问题（断电、坏盘、NFS 抖动）
4) 如果是 pack 校验失败，优先重新 `pack`（在保留 loose 的情况下）再验证

### 11.7 现象：远端 HPC retrieve 失败或文件不全

排查步骤：
1) 先在服务端确认输出文件确实存在于仓库工作目录（不是在别的 scratch 路径）
2) 检查 patterns：例如 `*.out,*.log` 是否匹配
3) 如果服务端返回 zip，但客户端解压失败：检查磁盘空间、目标目录权限

---

## 12. 参考资料（项目内可跳转）

- 项目入口与快速开始：[README.md](../README.md)
- 当前功能面（建议运维直接以此为准）：[FEATURE_STATUS.md](../FEATURE_STATUS.md)
- 本地仓库存储规范（理解 .chemvcs）：[docs/03-object-model-and-storage.md](../docs/03-object-model-and-storage.md)
- 远端协议（理解 push/pull/fetch 的 HTTP 语义）：[docs/05-remote-protocol-and-server.md](../docs/05-remote-protocol-and-server.md)
- HPC 用户指南（CLI + 网关）：[docs/10-hpc-user-guide.md](../docs/10-hpc-user-guide.md)
- 远端 HPC 网关实现说明（运维必读）：[docs/11-remote-hpc-design.md](../docs/11-remote-hpc-design.md)
- HPC 示例流程（可作为验收脚本）：[examples/hpc-workflow/README.md](../examples/hpc-workflow/README.md)
- ChemVCS 与 Git 差异理解：[docs/10-chemvcs-vs-git.md](../docs/10-chemvcs-vs-git.md)
