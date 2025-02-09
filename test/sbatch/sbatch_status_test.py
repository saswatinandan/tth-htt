#!/usr/bin/env python
from tthAnalysis.HiggsToTauTau.sbatchManager import sbatchManager, \
  sbatchManagerTimeoutError, sbatchManagerMemoryError, sbatchManagerSyntaxError, sbatchManagerSigBusError
from tthAnalysis.HiggsToTauTau.jobTools import create_if_not_exists

import os
import unittest
import shutil
import uuid
import subprocess

'''Tests sbatchManager submission and failure detection with a set erroneous jobs

The aim of this test is to verify that sbatchManager can deduce some of the typical errors a job might create.

The test cases are:
1) 1 job produces SIGBUS error
2) 1 job that consumes more memory than preset memory limit
3) 1 job that runs longer than the preset time limit
4) 1 job that has a syntax error in the batch code
Each job should raise an error in the sbatchManager instance that is specific to the created problem at hand. If any of
the jobs do not create such errors, the tests are bound to fail.

In order to run the test, navigate to the directory where this file resides and execute it:

$ ./sbatch_status_test.py

If you see the following message

> Ran 4 tests in ...s
>
> OK

then the tests passed.

'''

# use a directory which is universally available by cluster nodes and able to support writing in append mode, i.e. /home
testDir = os.path.expanduser('~/test_status_sbatch')

class SbatchStatusTestCase(unittest.TestCase):

  def __init__(self, *args, **kwargs):
    super(SbatchStatusTestCase, self).__init__(*args, **kwargs)
    self.manager        = None
    self.queue          = 'short'
    self.pool_id        = uuid.uuid4()
    self.verbose        = True
    self.log_completion = True
    self.workingDir     = os.getcwd()

  def setUp(self):
    self.manager = sbatchManager(pool_id = self.pool_id, verbose = self.verbose)
    self.manager.setWorkingDir(self.workingDir)
    self.manager.setLogFileDir(testDir)
    self.manager.queue = self.queue

    self.manager.log_completion = self.log_completion

  def tearDown(self):
    del self.manager

  def testBusError(self):
    # prepare executable that raises SIGBUS error
    cpp_source_name = os.path.join(testDir, 'bus_error.c')
    cmd = os.path.splitext(cpp_source_name)[0]
    with open(cpp_source_name, 'w') as cpp_source:
      cpp_source.write('#include <signal.h>\nint main() { raise(SIGBUS); }')
    compile_cmd_str = 'gcc -o %s %s' % (cmd, cpp_source_name)
    compile_cmd = subprocess.Popen(compile_cmd_str.split())
    compile_cmd.communicate()

    self.manager.poll_interval = 1
    self.manager.submitJob(
      inputFiles             = [],
      executable             = cmd,
      command_line_parameter = "",
      outputFilePath         = "",
      outputFiles            = [],
      scriptFile             = os.path.join(testDir, 'dummy_bus.sh'),
    )
    # if passes, true negative; otherwise true positive
    self.assertRaises(sbatchManagerSigBusError, self.manager.waitForJobs)

  def testMemory(self):
    cmd = 'python -c "`echo -e "a = []\\nwhile True:\\n\\ta.append(\' \' * 1024 * 1024)"`"'
    self.manager.max_mem = '2M'
    self.manager.poll_interval = 1

    self.manager.submitJob(
      inputFiles             = [],
      executable             = cmd,
      command_line_parameter = "",
      outputFilePath         = "",
      outputFiles            = [],
      scriptFile             = os.path.join(testDir, 'dummy_memory.sh'),
    )
    # if passes, true negative; otherwise true positive
    self.assertRaises(sbatchManagerMemoryError, self.manager.waitForJobs)

  def testTimeout(self):
    cmd = "python -c 'import time; time.sleep(100)'"
    self.manager.sbatchArgs = '--time=1'
    self.manager.poll_interval = 5

    self.manager.submitJob(
      inputFiles             = [],
      executable             = cmd,
      command_line_parameter = "",
      outputFilePath         = "",
      outputFiles            = [],
      scriptFile             = os.path.join(testDir, 'dummy_timeout.sh'),
    )
    # if passes, true negative; otherwise true positive
    self.assertRaises(sbatchManagerTimeoutError, self.manager.waitForJobs)

  def testSyntax(self):
    cmd = 'echo "missing quote'
    self.manager.poll_interval = 1

    self.manager.submitJob(
      inputFiles             = [],
      executable             = cmd,
      command_line_parameter = "",
      outputFilePath         = "",
      outputFiles            = [],
      scriptFile             = os.path.join(testDir, 'dummy_syntax.sh'),
    )
    # if passes, true negative; otherwise true positive
    self.assertRaises(sbatchManagerSyntaxError, self.manager.waitForJobs)

def suite():
  testSuite = unittest.TestSuite()
  testSuite.addTest(unittest.makeSuite(SbatchStatusTestCase))
  return testSuite

create_if_not_exists(testDir)

suite_instance = suite()
runner = unittest.TextTestRunner()
runner.run(suite_instance)

if os.path.isdir(testDir):
  shutil.rmtree(testDir)
