# ChemVCS Demo 手动演示指南

> **场景**: Si(硅) 体相 DFT 收敛性测试  
> **目标**: 逐步展示 ChemVCS 的所有 MVP 功能

---

## 准备工作

```powershell
# 创建一个临时演示目录
mkdir $env:TEMP\si_convergence
cd $env:TEMP\si_convergence
```

> 💡 以下所有命令都在这个目录下执行。

---

## Step 1: 初始化仓库

**讲解**: ChemVCS 为计算化学项目提供版本控制。第一步是初始化仓库。

```powershell
chemvcs init
```

**预期输出**: 显示仓库初始化成功，包含 `.chemvcs/` 目录路径和后续使用提示。

---

## Step 2: 创建初始 VASP 输入文件

**讲解**: 我们要对 Si 体相做 DFT 计算。先准备三个标准 VASP 输入文件：
- **POSCAR** — 硅的金刚石结构 (晶格常数 5.43 Å)
- **INCAR** — 计算参数 (ENCUT=400 eV, ISMEAR=0)
- **KPOINTS** — 4×4×4 Gamma k 点网格

```powershell
# 从预制文件复制
# ⚠️ 请将下面的路径替换为你的实际 chemvcs 项目路径
$DEMO = "C:\Users\lishi\code\chemvcs\demo\vasp_files"

Copy-Item "$DEMO\step1_initial\POSCAR" .
Copy-Item "$DEMO\step1_initial\INCAR"  .
Copy-Item "$DEMO\step1_initial\KPOINTS" .
```

可以看一下文件内容：

```powershell
type INCAR
```

**预期输出**:
```
# Si bulk - initial SCF calculation
SYSTEM  = Si-diamond
PREC    = Accurate
ENCUT   = 400
ISMEAR  = 0
SIGMA   = 0.05
...
```

---

## Step 3: 将文件添加到暂存区

**讲解**: `chemvcs add` 会计算文件内容哈希并识别 VASP 文件类型。

```powershell
chemvcs add POSCAR INCAR KPOINTS
```

**预期输出**: 显示 3 个文件被添加，每个文件显示大小和类型 (POSCAR / INCAR / KPOINTS)。

---

## Step 4: 查看仓库状态

**讲解**: 类似 `git status`，显示哪些文件已暂存、准备提交。

```powershell
chemvcs status
```

**预期输出**: 显示 3 个待提交的文件及其哈希值。

---

## Step 5: 创建初始提交

**讲解**: 将当前文件快照保存到版本历史中。

```powershell
chemvcs commit -m "Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)"
```

**预期输出**: 显示提交哈希、作者、时间。注意这是 root commit（没有父提交）。

---

## Step 6: 模拟计算完成 — 添加 OUTCAR

**讲解**: 假设 VASP 计算已完成，得到 OUTCAR。总能量 E = **-10.8265 eV**。我们把所有文件一起提交，保存完整的计算快照。

```powershell
Copy-Item "$DEMO\step2_scf_done\OUTCAR" .

# 查看能量结果
Select-String "energy.*sigma" OUTCAR
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
chemvcs commit -m "SCF completed: E=-10.8265 eV (ENCUT=400)"
```

**预期输出**: 4 个文件提交成功。

---

## Step 7: ENCUT 收敛性测试 (400 → 520 eV) ⭐

**讲解**: 能量截断 ENCUT 是 DFT 中最关键的参数。我们将 ENCUT 从 400 提高到 520 eV，看能量是否收敛。

```powershell
# 更新 INCAR 和 OUTCAR
Copy-Item "$DEMO\step3_encut_conv\INCAR"  . -Force
Copy-Item "$DEMO\step3_encut_conv\OUTCAR" . -Force

# 看看 INCAR 的变化
type INCAR
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
```

```powershell
chemvcs commit -m "ENCUT convergence: 400->520 eV, E=-10.8452 eV"
```

**⭐ 重点观察**: 提交输出中会出现 **Semantic Changes** 部分：
```
Semantic Changes:

  INCAR:
    ‼️  1 critical change(s)
    Key changes:
      ~ ENCUT: 400 → 520
```

**讲解要点**: ChemVCS 自动检测到 ENCUT 是 **critical** 级别的参数变化！这比 `git diff` 显示的纯文本差异有意义得多。

---

## Step 8: K-point 收敛性测试 (4×4×4 → 8×8×8) ⭐⭐

**讲解**: 接下来测试 k 点网格密度。同时调整 ISMEAR (0→1) 和 SIGMA (0.05→0.1)。

```powershell
# 更新 INCAR, KPOINTS, OUTCAR
Copy-Item "$DEMO\step4_kpoint_conv\INCAR"   . -Force
Copy-Item "$DEMO\step4_kpoint_conv\KPOINTS" . -Force
Copy-Item "$DEMO\step4_kpoint_conv\OUTCAR"  . -Force
```

```powershell
chemvcs add POSCAR INCAR KPOINTS OUTCAR
```

```powershell
chemvcs commit -m "K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV"
```

**⭐⭐ 重点观察**: 同时检测到 INCAR 和 KPOINTS 的语义变化：
```
Semantic Changes:

  INCAR:
    ‼️  5 critical change(s)
    ⚠️  1 major change(s)
    Key changes:
      ~ ISMEAR: 0 → 1
      ~ SIGMA: 0.05 → 0.1
      ~ EDIFF: 1e-06 → 1e-07
    ... and 3 more change(s)

  KPOINTS:
    ‼️  1 critical change(s)
    Key changes:
      ~ grid: [4, 4, 4] → [8, 8, 8]
```

**讲解要点**:
- INCAR 中 5 个 critical 级别变化 + 1 个 major 变化
- KPOINTS 的 k 点网格变化被自动标记为 critical
- 帮助研究者快速判断：这次修改会显著影响计算结果

---

## Step 9: 查看提交历史

**讲解**: 查看完整的计算版本历史，记录了每次参数调整。

```powershell
chemvcs log
```

简洁模式：

```powershell
chemvcs log --oneline
```

**预期输出** (类似):
```
d77f4cb K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV
3fcec3e ENCUT convergence: 400->520 eV, E=-10.8452 eV
9efe155 SCF completed: E=-10.8265 eV (ENCUT=400)
41c41d0 Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)
```

**讲解要点**: 可以清楚看到能量从 -10.8265 → -10.8452 → -10.8489 eV 逐步收敛。

---

## Step 10: 语义 Diff ⭐⭐⭐

**讲解**: 这是 ChemVCS 的核心特色——对 VASP 文件做语义级别的比较。

### 10a. 详细 diff（默认）

```powershell
chemvcs diff
```

**预期输出**: 显示 HEAD 和上一个提交的完整语义差异，包括每个参数的 ‼️ critical / ⚠️ major 标记。

### 10b. 摘要模式

```powershell
chemvcs diff --summary
```

**预期输出**: 只显示变化数量统计。

### 10c. JSON 输出（机器可读）

```powershell
chemvcs diff --format json
```

**讲解要点**: JSON 格式可以被其他工具解析，方便做自动化分析,比如批量比较收敛性测试的参数变化。

---

## Step 11: 复现历史计算 ⭐⭐⭐

**讲解**: 最核心的功能之一 —— 给定提交哈希，完整复现该版本的计算输入。

```powershell
# 先找到初始提交的哈希
chemvcs log --oneline
```

```powershell
# 用初始提交的哈希来复现 (替换为你实际看到的哈希)
chemvcs reproduce <初始提交哈希> -o reproduce_initial
```

```powershell
# 查看复现的文件
Get-ChildItem reproduce_initial

# 验证 INCAR 是初始版本 (ENCUT=400)
type reproduce_initial\INCAR
```

**讲解要点**:
- 完整导出该版本的所有输入文件到新目录
- 确保计算 100% 可重复
- 对论文审稿、结果验证至关重要

---

## Step 12: 总结

| 功能 | 命令 | 解决的问题 |
|------|------|-----------|
| 初始化 | `chemvcs init` | 为计算项目建立版本控制 |
| 追踪文件 | `chemvcs add` | 自动识别 VASP 文件类型 |
| 查看状态 | `chemvcs status` | 了解工作区的变化 |
| 提交快照 | `chemvcs commit` | 保存计算参数+结果的完整快照 |
| 语义 Diff | `chemvcs diff` | **理解**参数变化，不只是看文本差异 |
| 版本历史 | `chemvcs log` | 追踪计算演化过程 |
| 复现计算 | `chemvcs reproduce` | 确保计算可重复 |

### 能量收敛过程

```
ENCUT=400, K=4×4×4  →  E = -10.8265 eV
ENCUT=520, K=4×4×4  →  E = -10.8452 eV  (ΔE = 18.7 meV)
ENCUT=520, K=8×8×8  →  E = -10.8489 eV  (ΔE =  3.7 meV → 收敛!)
```

---

## 清理

```powershell
cd ~
Remove-Item $env:TEMP\si_convergence -Recurse -Force
```

---

## 文件结构说明

```
demo/
├── GUIDE.md                      ← 本文件（操作手册）
├── run_demo.py                   ← 自动化演示脚本（备用）
└── vasp_files/                   ← 各阶段预制的 VASP 文件
    ├── step1_initial/            ← 初始输入
    │   ├── POSCAR                   Si 金刚石结构
    │   ├── INCAR                    ENCUT=400, ISMEAR=0
    │   └── KPOINTS                  4×4×4 Gamma
    ├── step2_scf_done/           ← 计算完成后新增
    │   └── OUTCAR                   E=-10.8265 eV
    ├── step3_encut_conv/         ← ENCUT 收敛测试
    │   ├── INCAR                    ENCUT=520 (变化!)
    │   └── OUTCAR                   E=-10.8452 eV
    └── step4_kpoint_conv/        ← K-point 收敛测试
        ├── INCAR                    ISMEAR=1, SIGMA=0.1 (多处变化!)
        ├── KPOINTS                  8×8×8 (变化!)
        └── OUTCAR                   E=-10.8489 eV
```
