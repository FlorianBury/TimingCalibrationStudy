#!/usr/bin/env python 
# Copyright (c) 2022, Université catholique de Louvain 
# 
# Distributed under the GNU General Public License, see accompanying file LICENSE 
# or https://gitlab.cern.ch/cp3-cms/pyplotit for details. 

from setuptools import setup 

setup(
    entry_points = {
        'console_scripts': [
        'BXRun=timeCal.utils.BXRun:main',
        'plotScan=timeCal.PlotScan.plotScan:main',
        'runCalibration=timeCal.Calibration.runCalibration:main',
        ],
    },
    install_requires = [
        'numpy', 
        'pandas', 
        'enlighten',
        'pyyaml',
        'ipython',
        'tk',
    ]
)

