#!/usr/bin/env python3
"""
ChemVCS Interactive Demo
========================

Simulates a real VASP workflow: Si bulk ENCUT & k-point convergence test.

Demonstrates:
  - chemvcs init
  - chemvcs add / status / commit
  - Semantic diff for INCAR & KPOINTS
  - chemvcs log / diff / reproduce

Usage:
    python demo/run_demo.py
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

# ── ANSI Colors ──────────────────────────────────────────────────────────────

BOLD = "\033[1m"
DIM = "\033[2m"
GREEN = "\033[92m"
CYAN = "\033[96m"
YELLOW = "\033[93m"
MAGENTA = "\033[95m"
RED = "\033[91m"
RESET = "\033[0m"
UNDERLINE = "\033[4m"

# Detect Windows console
if platform.system() == "Windows":
    os.system("")  # enable ANSI on Windows


# ── Helpers ──────────────────────────────────────────────────────────────────

def banner(text: str) -> None:
    """Print a prominent section banner."""
    width = 72
    print()
    print(f"{MAGENTA}{'═' * width}")
    print(f"  {BOLD}{text}{RESET}{MAGENTA}")
    print(f"{'═' * width}{RESET}")
    print()


def narrate(text: str) -> None:
    """Print narrative / explanation text."""
    for line in textwrap.wrap(text, width=70):
        print(f"  {DIM}{line}{RESET}")
    print()


def show_command(cmd: str) -> None:
    """Display the command that is about to run."""
    print(f"  {GREEN}${RESET} {BOLD}{cmd}{RESET}")


def pause() -> None:
    """Wait for user to press Enter."""
    try:
        input(f"\n  {DIM}[ 按 Enter 继续 ]{RESET}")
    except (EOFError, KeyboardInterrupt):
        print(f"\n{YELLOW}Demo terminated.{RESET}")
        sys.exit(0)
    print()


def run(cmd: str, cwd: Path, *, show: bool = True) -> str:
    """Run a shell command and return stdout."""
    if show:
        show_command(cmd)
        print()

    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"

    result = subprocess.run(
        cmd,
        shell=True,
        cwd=cwd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        env=env,
    )
    output = result.stdout
    if result.stderr and result.returncode != 0:
        output += result.stderr
    if output.strip():
        for line in output.rstrip("\n").split("\n"):
            print(f"  {line}")
        print()
    return output


def write_file(path: Path, content: str) -> None:
    """Write a file and print a short indicator."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"    {CYAN}✎{RESET} 写入 {path.name}  ({len(content)} bytes)")


# ── VASP File Templates ─────────────────────────────────────────────────────

POSCAR_SI = """\
Si2 Diamond
5.43
  0.5  0.5  0.0
  0.0  0.5  0.5
  0.5  0.0  0.5
Si
2
Direct
  0.00  0.00  0.00
  0.25  0.25  0.25
"""

INCAR_V1 = """\
# Si bulk - initial SCF calculation
SYSTEM  = Si-diamond
PREC    = Accurate
ENCUT   = 400
ISMEAR  = 0
SIGMA   = 0.05
EDIFF   = 1E-6
NELM    = 100
LREAL   = .FALSE.
LWAVE   = .TRUE.
LCHARG  = .TRUE.
"""

KPOINTS_V1 = """\
Automatic mesh
0
Gamma
4 4 4
0 0 0
"""

OUTCAR_V1 = """\
 running on    4 total cores
 vasp.6.3.2 2022

 POSCAR found type information on POSCAR: Si
 POSCAR found :  2 types and       2 ions

 energy  without entropy=     -10.82654312  energy(sigma->0) =     -10.82654312

 General timing and accounting informations for this job:
  Total CPU time used (sec):       42.358
  Elapsed time (sec):              45.123
  Maximum memory used (kb):      125408
"""

INCAR_V2 = """\
# Si bulk - ENCUT convergence test
SYSTEM  = Si-diamond
PREC    = Accurate
ENCUT   = 520
ISMEAR  = 0
SIGMA   = 0.05
EDIFF   = 1E-6
NELM    = 100
LREAL   = .FALSE.
LWAVE   = .TRUE.
LCHARG  = .TRUE.
"""

OUTCAR_V2 = """\
 running on    4 total cores
 vasp.6.3.2 2022

 POSCAR found type information on POSCAR: Si
 POSCAR found :  2 types and       2 ions

 energy  without entropy=     -10.84521078  energy(sigma->0) =     -10.84521078

 General timing and accounting informations for this job:
  Total CPU time used (sec):       58.214
  Elapsed time (sec):              62.001
  Maximum memory used (kb):      138720
"""

KPOINTS_V2 = """\
Automatic mesh
0
Gamma
8 8 8
0 0 0
"""

INCAR_V3 = """\
# Si bulk - final converged parameters
SYSTEM  = Si-diamond
PREC    = Accurate
ENCUT   = 520
ISMEAR  = 1
SIGMA   = 0.1
EDIFF   = 1E-7
NELM    = 200
LREAL   = .FALSE.
LWAVE   = .FALSE.
LCHARG  = .TRUE.
NSW     = 0
"""

OUTCAR_V3 = """\
 running on    4 total cores
 vasp.6.3.2 2022

 POSCAR found type information on POSCAR: Si
 POSCAR found :  2 types and       2 ions

 energy  without entropy=     -10.84893256  energy(sigma->0) =     -10.84890143

 General timing and accounting informations for this job:
  Total CPU time used (sec):      124.875
  Elapsed time (sec):             130.422
  Maximum memory used (kb):      189504
"""


# ── Demo Steps ───────────────────────────────────────────────────────────────

def step_0_intro(demo_dir: Path) -> None:
    """Show introduction."""
    print(f"""
{BOLD}{CYAN}
   ╔═══════════════════════════════════════════════════════════════╗
   ║                                                               ║
   ║          ChemVCS  ─  化学计算版本控制系统                     ║
   ║          Interactive Demo                                     ║
   ║                                                               ║
   ╚═══════════════════════════════════════════════════════════════╝
{RESET}
  {DIM}场景: Si(硅) 体相 DFT 收敛性测试{RESET}
  {DIM}工作流: 初始计算 → ENCUT 收敛 → k-point 收敛 → 复现{RESET}

  {DIM}Demo 工作目录: {demo_dir}{RESET}
""")


def step_1_init(demo_dir: Path) -> None:
    """Initialize ChemVCS repository."""
    banner("Step 1: 初始化 ChemVCS 仓库")
    narrate(
        "首先，我们在项目目录中初始化一个 ChemVCS 仓库。"
        "这会创建 .chemvcs/ 目录，用于存储对象、提交记录和元数据。"
    )
    run("chemvcs init", demo_dir)


def step_2_create_inputs(demo_dir: Path) -> None:
    """Create initial VASP input files."""
    banner("Step 2: 创建 VASP 输入文件")
    narrate(
        "我们正在进行 Si 体相的 DFT 计算。创建标准 VASP 输入文件："
        "POSCAR (晶体结构)、INCAR (计算参数)、KPOINTS (k点网格)。"
        "初始参数: ENCUT=400 eV, K-points 4×4×4 Gamma 网格。"
    )
    write_file(demo_dir / "POSCAR", POSCAR_SI)
    write_file(demo_dir / "INCAR", INCAR_V1)
    write_file(demo_dir / "KPOINTS", KPOINTS_V1)
    print()

    narrate("查看 INCAR 的内容:")
    run("type INCAR" if platform.system() == "Windows" else "cat INCAR", demo_dir)


def step_3_add_commit_initial(demo_dir: Path) -> None:
    """Add and commit initial files."""
    banner("Step 3: 添加文件并创建初始提交")
    narrate(
        "使用 chemvcs add 将文件添加到暂存区。"
        "ChemVCS 会自动识别 VASP 文件类型并计算内容哈希。"
    )
    run("chemvcs add POSCAR INCAR KPOINTS", demo_dir)
    pause()

    narrate("查看当前状态 — 显示已暂存的文件:")
    run("chemvcs status", demo_dir)
    pause()

    narrate("创建初始提交:")
    run('chemvcs commit -m "Initial setup: Si bulk SCF (ENCUT=400, K=4x4x4)"', demo_dir)


def step_4_add_output(demo_dir: Path) -> None:
    """Simulate calculation completion and add output."""
    banner("Step 4: 模拟计算完成，提交结果")
    narrate(
        "假设 VASP 计算已完成。我们得到了 OUTCAR 输出文件，"
        "其中包含总能量 E = -10.8265 eV。将所有文件一起提交，"
        "确保每个版本都包含完整的计算快照。"
    )
    write_file(demo_dir / "OUTCAR", OUTCAR_V1)
    print()

    run("chemvcs add POSCAR INCAR KPOINTS OUTCAR", demo_dir)
    run('chemvcs commit -m "SCF completed: E=-10.8265 eV (ENCUT=400)"', demo_dir)


def step_5_encut_convergence(demo_dir: Path) -> None:
    """Modify ENCUT for convergence test."""
    banner("Step 5: ENCUT 收敛性测试 (400 → 520 eV)")
    narrate(
        "能量截断(ENCUT)是 DFT 计算中最关键的参数之一。"
        "我们将 ENCUT 从 400 eV 提高到 520 eV，重新提交。"
        "注意观察 ChemVCS 的语义 diff —— 它会告诉你 ENCUT 是一个 critical 参数！"
    )

    write_file(demo_dir / "INCAR", INCAR_V2)
    write_file(demo_dir / "OUTCAR", OUTCAR_V2)
    print()

    run("chemvcs add POSCAR INCAR KPOINTS OUTCAR", demo_dir)
    pause()

    narrate("提交 — 注意观察 Semantic Changes 部分:")
    run('chemvcs commit -m "ENCUT convergence: 400->520 eV, E=-10.8452 eV"', demo_dir)


def step_6_kpoint_convergence(demo_dir: Path) -> None:
    """Modify KPOINTS for convergence test."""
    banner("Step 6: K-point 收敛性测试 (4×4×4 → 8×8×8)")
    narrate(
        "接下来测试 k 点网格密度。将 4×4×4 Gamma 网格增加到 8×8×8。"
        "同时调整 ISMEAR 和 SIGMA 以适应金属/半导体系统。"
        "ChemVCS 会同时检测 INCAR 和 KPOINTS 的语义变化。"
    )

    write_file(demo_dir / "INCAR", INCAR_V3)
    write_file(demo_dir / "KPOINTS", KPOINTS_V2)
    write_file(demo_dir / "OUTCAR", OUTCAR_V3)
    print()

    run("chemvcs add POSCAR INCAR KPOINTS OUTCAR", demo_dir)
    pause()

    narrate("提交 — 观察 INCAR 和 KPOINTS 的多个语义变化:")
    run(
        'chemvcs commit -m "K-point convergence: 8x8x8, ISMEAR=1, E=-10.8489 eV"',
        demo_dir,
    )


def step_7_log(demo_dir: Path) -> None:
    """Show commit history."""
    banner("Step 7: 查看提交历史")
    narrate(
        "使用 chemvcs log 查看完整的计算版本历史。"
        "每个提交都记录了参数快照和结果。"
    )
    run("chemvcs log", demo_dir)
    pause()

    narrate("也可以用 --oneline 格式查看简洁版本:")
    run("chemvcs log --oneline", demo_dir)


def step_8_diff(demo_dir: Path) -> None:
    """Show semantic diff between revisions."""
    banner("Step 8: 语义 Diff — 比较版本差异")
    narrate(
        "chemvcs diff 会对 VASP 文件进行语义分析，"
        "自动标记 critical/major/minor 级别的参数变化，"
        "帮助研究者快速理解两次计算之间的区别。"
    )
    run("chemvcs diff", demo_dir)
    pause()

    narrate("使用 --summary 查看精简摘要:")
    run("chemvcs diff --summary", demo_dir)
    pause()

    narrate("使用 --format json 输出机器可读格式:")
    run("chemvcs diff --format json", demo_dir)


def step_9_reproduce(demo_dir: Path) -> str:
    """Reproduce an earlier calculation."""
    banner("Step 9: 复现历史计算")
    narrate(
        "这是 ChemVCS 最核心的功能之一：复现！"
        "给定某个提交的哈希值，chemvcs reproduce 会将该版本的所有输入文件"
        "完整导出到指定目录，确保计算完全可重复。"
    )

    # Get the first commit hash (Initial setup)
    result = run("chemvcs log --oneline", demo_dir, show=False)
    lines = [l.strip() for l in result.strip().split("\n") if l.strip()]
    # Get the last line (oldest commit = Initial setup)
    first_commit_line = lines[-1] if lines else ""
    first_hash = first_commit_line.split()[0] if first_commit_line else ""

    if first_hash:
        narrate(
            f"让我们复现最初的计算设置 (commit {first_hash[:7]}...)。"
            "这会将 ENCUT=400, K=4×4×4 的初始文件导出到新目录。"
        )
        run(
            f"chemvcs reproduce {first_hash} -o reproduce_initial",
            demo_dir,
        )
        pause()

        narrate("查看复现出的文件:")
        if platform.system() == "Windows":
            run("type reproduce_initial\\INCAR", demo_dir)
        else:
            run("cat reproduce_initial/INCAR", demo_dir)
    else:
        narrate("(无法获取初始提交哈希)")

    return first_hash


def step_10_conclusion(demo_dir: Path) -> None:
    """Show conclusion."""
    print(f"""
{BOLD}{CYAN}
   ╔═══════════════════════════════════════════════════════════════╗
   ║                                                               ║
   ║          Demo 完成!                                           ║
   ║                                                               ║
   ╚═══════════════════════════════════════════════════════════════╝
{RESET}

  {BOLD}ChemVCS MVP 功能总结:{RESET}

    {GREEN}✓{RESET} {BOLD}chemvcs init{RESET}       — 初始化仓库
    {GREEN}✓{RESET} {BOLD}chemvcs add{RESET}        — 添加文件到暂存区 (自动识别VASP文件类型)
    {GREEN}✓{RESET} {BOLD}chemvcs status{RESET}     — 查看工作区状态
    {GREEN}✓{RESET} {BOLD}chemvcs commit{RESET}     — 创建提交 (含语义 diff 摘要)
    {GREEN}✓{RESET} {BOLD}chemvcs log{RESET}        — 查看提交历史
    {GREEN}✓{RESET} {BOLD}chemvcs diff{RESET}       — 语义差异分析 (critical/major/minor)
    {GREEN}✓{RESET} {BOLD}chemvcs reproduce{RESET}  — 复现历史计算

  {BOLD}语义分析支持:{RESET}

    {CYAN}•{RESET} INCAR  — 参数级别变化检测 (ENCUT, ISMEAR, SIGMA...)
    {CYAN}•{RESET} KPOINTS — k点网格/类型变化检测

  {DIM}Demo 工作目录: {demo_dir}{RESET}
  {DIM}你可以进入该目录继续体验 chemvcs 命令。{RESET}
""")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    """Run the interactive demo."""
    # Create demo workspace
    demo_dir = Path(tempfile.mkdtemp(prefix="chemvcs_demo_"))

    try:
        step_0_intro(demo_dir)
        pause()

        step_1_init(demo_dir)
        pause()

        step_2_create_inputs(demo_dir)
        pause()

        step_3_add_commit_initial(demo_dir)
        pause()

        step_4_add_output(demo_dir)
        pause()

        step_5_encut_convergence(demo_dir)
        pause()

        step_6_kpoint_convergence(demo_dir)
        pause()

        step_7_log(demo_dir)
        pause()

        step_8_diff(demo_dir)
        pause()

        step_9_reproduce(demo_dir)
        pause()

        step_10_conclusion(demo_dir)

    except KeyboardInterrupt:
        print(f"\n\n{YELLOW}Demo interrupted.{RESET}")
        print(f"{DIM}Workspace: {demo_dir}{RESET}\n")


if __name__ == "__main__":
    main()
