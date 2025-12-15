"""Environment provenance capture for computational reproducibility."""

import os
import re
import subprocess
from typing import Dict, List, Optional


class EnvironmentCapture:
    """Capture computational environment for provenance tracking."""
    
    @staticmethod
    def capture_modules() -> List[str]:
        """
        Capture loaded environment modules.
        
        Returns:
            List of module names with versions (e.g., ['vasp/6.3.0', 'intel/2021.4'])
        """
        try:
            # Environment modules output to stderr
            result = subprocess.run(
                ['module', 'list', '-t'],
                capture_output=True,
                text=True,
                stderr=subprocess.STDOUT,
                timeout=5
            )
            
            modules = []
            for line in result.stdout.strip().split('\n'):
                line = line.strip()
                # Skip header lines
                if line and not line.startswith('Currently') and not line.startswith('No'):
                    modules.append(line)
            
            return modules
        except (FileNotFoundError, subprocess.TimeoutExpired):
            # Module command not available
            return []
    
    @staticmethod
    def capture_env_vars(keys: Optional[List[str]] = None) -> Dict[str, str]:
        """
        Capture environment variables.
        
        Args:
            keys: Specific variables to capture. If None, captures common ones.
            
        Returns:
            Dictionary of environment variables
        """
        if keys is None:
            # Common HPC/computational environment variables
            keys = [
                'OMP_NUM_THREADS',
                'MKL_NUM_THREADS',
                'SLURM_JOB_ID',
                'SLURM_JOB_NAME',
                'SLURM_NTASKS',
                'PBS_JOBID',
                'PBS_JOBNAME',
                'PATH',
                'LD_LIBRARY_PATH',
            ]
        
        env_vars = {}
        for key in keys:
            value = os.environ.get(key)
            if value:
                env_vars[key] = value
        
        return env_vars
    
    @staticmethod
    def parse_slurm_script(script_path: str) -> Dict[str, any]:
        """
        Extract resource requirements from SLURM job script.
        
        Args:
            script_path: Path to SLURM script
            
        Returns:
            Dictionary with resource requirements
        """
        resources = {}
        
        try:
            with open(script_path, 'r') as f:
                content = f.read()
            
            # Parse #SBATCH directives
            for line in content.split('\n'):
                line = line.strip()
                if not line.startswith('#SBATCH'):
                    continue
                
                # Remove #SBATCH prefix
                directive = line[7:].strip()
                
                # Parse different directives
                if match := re.match(r'--nodes?=(\d+)', directive):
                    resources['nodes'] = int(match.group(1))
                elif match := re.match(r'--ntasks-per-node=(\d+)', directive):
                    resources['ntasks_per_node'] = int(match.group(1))
                elif match := re.match(r'--ntasks?=(\d+)', directive):
                    resources['ntasks'] = int(match.group(1))
                elif match := re.match(r'--cpus-per-task=(\d+)', directive):
                    resources['cpus_per_task'] = int(match.group(1))
                elif match := re.match(r'--time=(.+)', directive):
                    resources['walltime'] = match.group(1)
                elif match := re.match(r'--mem=(.+)', directive):
                    resources['memory'] = match.group(1)
                elif match := re.match(r'--mem-per-cpu=(.+)', directive):
                    resources['memory_per_cpu'] = match.group(1)
                elif match := re.match(r'--partition=(.+)', directive):
                    resources['partition'] = match.group(1)
                elif match := re.match(r'--job-name=(.+)', directive):
                    resources['job_name'] = match.group(1)
        
        except FileNotFoundError:
            pass  # Script doesn't exist, return empty resources
        
        return resources
    
    @staticmethod
    def parse_pbs_script(script_path: str) -> Dict[str, any]:
        """
        Extract resource requirements from PBS/Torque job script.
        
        Args:
            script_path: Path to PBS script
            
        Returns:
            Dictionary with resource requirements
        """
        resources = {}
        
        try:
            with open(script_path, 'r') as f:
                content = f.read()
            
            # Parse #PBS directives
            for line in content.split('\n'):
                line = line.strip()
                if not line.startswith('#PBS'):
                    continue
                
                # Remove #PBS prefix
                directive = line[4:].strip()
                
                # Parse different directives
                if match := re.match(r'-l nodes=(\d+):ppn=(\d+)', directive):
                    resources['nodes'] = int(match.group(1))
                    resources['ppn'] = int(match.group(2))
                elif match := re.match(r'-l walltime=(.+)', directive):
                    resources['walltime'] = match.group(1)
                elif match := re.match(r'-l mem=(.+)', directive):
                    resources['memory'] = match.group(1)
                elif match := re.match(r'-q (.+)', directive):
                    resources['queue'] = match.group(1)
                elif match := re.match(r'-N (.+)', directive):
                    resources['job_name'] = match.group(1)
        
        except FileNotFoundError:
            pass
        
        return resources
    
    @staticmethod
    def save_script_snapshot(script_path: str) -> str:
        """
        Read and return job script content.
        
        Args:
            script_path: Path to job script
            
        Returns:
            Full script content as string
        """
        try:
            with open(script_path, 'r') as f:
                return f.read()
        except FileNotFoundError:
            return ""
    
    @staticmethod
    def detect_script_type(script_path: str) -> Optional[str]:
        """
        Detect scheduler type from script content.
        
        Args:
            script_path: Path to job script
            
        Returns:
            Scheduler type ('slurm', 'pbs', 'sge') or None
        """
        try:
            with open(script_path, 'r') as f:
                content = f.read()
            
            if '#SBATCH' in content:
                return 'slurm'
            elif '#PBS' in content:
                return 'pbs'
            elif '#$' in content:
                return 'sge'
        except FileNotFoundError:
            pass
        
        return None
