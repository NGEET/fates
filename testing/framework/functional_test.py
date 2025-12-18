from abc import ABC, abstractmethod
from pathlib import Path
from framework.fates_test import FatesTest
from framework.utils.general import str_to_bool, str_to_list, copy_file


class FunctionalTest(FatesTest):
    """Class for running FATES functional tests"""

    def __init__(self, config: dict):
        super().__init__(name=config['name'], test_dir=config['test_dir'],
                         test_exe=config['test_exe'])
        
        self.out_file = config['out_file']
        self.use_param_file = str_to_bool(config['use_param_file'])
        self.other_args = str_to_list(config['other_args'])
        self.plot = True
        
    def run_cmd(self, param_file) -> list[str]:
        """Returns the run command"""
        base_cmd = [f"./{self.test_exe}"]
        if self.use_param_file:
            base_cmd.append(str(param_file))
        base_cmd.extend(self.other_args)
        
    def run(self, build_dir: Path, run_dir: Path, param_file: Path):
        """Execute a functional test"""
        
        # find executable
        exe_path = self.find_build(build_dir)
        
        # copy file to run directory
        copy_file(exe_path, run_dir)
        local_exe = run_dir / self.test_exe
        
        cmd = run_cmd()

        
        return self.execute_shell(cmd, run_dir)
        
    @abstractmethod
    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        pass
