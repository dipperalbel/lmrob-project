import subprocess
import os
import sys
from setuptools import setup
from distutils.command.build_py import build_py as _build_py
class build_py(_build_py):
    """Specialized Python source builder."""
    def run(self):
        command = "./tools/lmrob_lib_so/compile_lmrob_so.sh"
        p = subprocess.Popen(command, shell=True)
        p.wait()
        if not os.path.exists("lmrob/liblmrob.so"):
            raise Exception("liblmrob.so not found. Verify the requirements")
        super().run()
    

    

setup(
    name='Lmrob',
    packages=['lmrob'],
    package_dir={'lmrob': 'lmrob'},
    package_data={'lmrob': ['liblmrob.so']},
    version='1.0',
    author='',
    cmdclass={'build_py': build_py},
    install_requires=[package for package in open("requirements.txt")],
    )

