#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: josue_nataren
"""

import os
import subprocess
import time

def test_accuracy(capfd):
    result = subprocess.run(['python', 'EKGclassification.py'], capture_output=True, text=True)
    output = result.stdout
    
    assert 'Training accuracy higher than 90%' in output, "Training FAILED"
    assert 'Validation accuracy higher than 90%' in output, "Validation FAILED"
    assert 'Testing accuracy higher than 90%' in output, "Testing FAILED"
    print("Accuracy test PASSED!")

def test_runtime(capfd):
    start_time = time.time()
    subprocess.run(['python', 'EKGclassification.py'])
    end_time = time.time()

    assert end_time - start_time < 600, "Runtime FAILED"
    print("Runtime test PASSED!")
    
