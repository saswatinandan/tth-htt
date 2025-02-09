from tthAnalysis.HiggsToTauTau.safe_root import ROOT
from tthAnalysis.HiggsToTauTau.jobTools import run_cmd
from tthAnalysis.HiggsToTauTau.analysisTools import createFile
from tthAnalysis.HiggsToTauTau.common import logging

import os
import jinja2
import math
import uuid

MIN_FILE_SIZE_DEFAULT = 5000
MAX_NOF_RUNNING_JOBS_DEFAULT = 5000
MAX_INPUTS_PER_JOB_DEFAULT = 10

def find_duplicates(input_file_names):
  input_file_names_occurrence = {}
  for input_file_name in input_file_names:
    if input_file_name not in input_file_names_occurrence:
      input_file_names_occurrence[input_file_name] = 0
    input_file_names_occurrence[input_file_name] += 1
  input_file_names_duplicates = [
    input_file_name for input_file_name in input_file_names_occurrence \
    if input_file_names_occurrence[input_file_name] > 1
  ]
  return input_file_names_duplicates

def createScript_sbatch(
    sbatch_script_file_name,
    executable,
    command_line_parameters,
    input_file_names,
    output_file_names,
    script_file_names,
    log_file_names        = None,
    keep_logs             = False,
    working_dir           = None,
    max_num_jobs          = -1,
    cvmfs_error_log       = None,
    pool_id               = '',
    cmssw_base_dir        = None,
    verbose               = False,
    job_template_file     = 'sbatch-node.sh.template',
    dry_run               = False,
    validate_outputs      = True,
    min_file_size         = MIN_FILE_SIZE_DEFAULT,
    max_num_submittedJobs = MAX_NOF_RUNNING_JOBS_DEFAULT,
    use_home              = False,
    copy_output_file      = True,
  ):
    """Creates the python script necessary to submit analysis and/or Ntuple production jobs to the batch system
    """
    if not working_dir:
        working_dir = os.getcwd()

    if not cmssw_base_dir:
        cmssw_base_dir = os.environ['CMSSW_BASE']

    if not pool_id:
        raise ValueError('pool_id is empty')
    sbatch_analyze_lines, num_jobs = generate_sbatch_lines(
        executable              = executable,
        command_line_parameters = command_line_parameters,
        input_file_names        = input_file_names,
        output_file_names       = output_file_names,
        script_file_names       = script_file_names,
        log_file_names          = log_file_names,
        keep_logs               = keep_logs,
        working_dir             = working_dir,
        max_num_jobs            = max_num_jobs,
        cvmfs_error_log         = cvmfs_error_log,
        pool_id                 = pool_id,
        cmssw_base_dir          = cmssw_base_dir,
        verbose                 = verbose,
        job_template_file       = job_template_file,
        dry_run                 = dry_run,
        validate_outputs        = validate_outputs,
        min_file_size           = min_file_size,
        max_num_submittedJobs   = max_num_submittedJobs,
        use_home                = use_home,
        copy_output_file        = copy_output_file,
    )
    createFile(sbatch_script_file_name, sbatch_analyze_lines)
    return num_jobs

def generate_sbatch_lines(
    executable,
    command_line_parameters,
    input_file_names,
    output_file_names,
    script_file_names,
    log_file_names,
    keep_logs,
    working_dir,
    max_num_jobs,
    cvmfs_error_log       = None,
    pool_id               = '',
    cmssw_base_dir        = None,
    verbose               = False,
    job_template_file     = 'sbatch-node.sh.template',
    dry_run               = False,
    validate_outputs      = True,
    min_file_size         = MIN_FILE_SIZE_DEFAULT,
    max_num_submittedJobs = MAX_NOF_RUNNING_JOBS_DEFAULT,
    use_home              = False,
    copy_output_file      = True,
  ):
    if not pool_id:
        raise ValueError('pool_id is empty')
    lines_sbatch = [
        "from tthAnalysis.HiggsToTauTau.sbatchManager import sbatchManager",
        "",
        "m = sbatchManager('%s', verbose = %s, dry_run = %s, use_home = %s, min_file_size = %d, max_num_submittedJobs = %d)" % \
          (pool_id, verbose, dry_run, use_home, min_file_size, max_num_submittedJobs),
        "m.setWorkingDir('%s')"     % working_dir,
        "m.setcmssw_base_dir('%s')" % cmssw_base_dir,
        "m.log_completion = %s"     % verbose,
    ]

    num_jobs = 0
    for key_file, command_line_parameter in command_line_parameters.items():
        log_file_name = None
        if log_file_names:
            log_file_name = log_file_names[key_file]
        if num_jobs <= max_num_jobs or max_num_jobs <= 0:
            keep_logs_arg = False
            if keep_logs:
              if type(keep_logs) == dict and keep_logs[key_file]:
                keep_logs_arg = True
              elif type(keep_logs) == bool:
                keep_logs_arg = keep_logs
            sbatch_line = generate_sbatch_line(
                executable = executable,
                command_line_parameter = command_line_parameter,
                input_file_names       = input_file_names[key_file],
                output_file_name       = output_file_names[key_file],
                script_file_name       = script_file_names[key_file],
                log_file_name          = log_file_name,
                keep_logs              = keep_logs_arg,
                cvmfs_error_log        = cvmfs_error_log,
                job_template_file      = job_template_file,
                validate_outputs       = validate_outputs,
                min_file_size          = min_file_size,
                copy_output_file       = copy_output_file,
            )
            if sbatch_line:
                lines_sbatch.append(sbatch_line)
                num_jobs = num_jobs + 1
    if max_num_jobs > 0 and num_jobs > max_num_jobs:
        logging.warning(
          "number of jobs = %i exceeds limit of %i --> skipping submission of %i jobs !!" % \
          (num_jobs, max_num_jobs, num_jobs - max_num_jobs)
        )

    lines_sbatch.append("m.waitForJobs(%s)" %validate_outputs)
    return lines_sbatch, num_jobs

def is_file_ok(output_file_name, validate_outputs = True, min_file_size = MIN_FILE_SIZE_DEFAULT):
  if not (output_file_name and os.path.exists(output_file_name)):
    return False

  logging.info("Output file %s already exists" % output_file_name)

  if not output_file_name.lower().endswith('.root'):
    return True

  command = "rm %s" % output_file_name
  ret_value = False
  if min_file_size > 0:
    output_file_size = os.stat(output_file_name).st_size
    if output_file_size > min_file_size:
      if not validate_outputs:
        ret_value = True
    else:
      logging.info(
        "Deleting output file and resubmitting job because the output file has size of %d bytes, smaller than %d bytes, indicating that a failure in writing the output file occured." % (output_file_size, min_file_size)
      )

  if validate_outputs:
    root_tfile = ROOT.TFile(output_file_name, "read")
    if not root_tfile:
      logging.info("Not a valid ROOT file, deleting it")
    else:
      if root_tfile.IsZombie():
        logging.info("Output file is corrupted, deleting file and resubmitting job")
      else:
        # Let's open the file via bash as well to see if ROOT tries to recover the file
        open_cmd = "root -b -l -q %s 2>&1 > /dev/null | grep 'trying to recover' | wc -l" % output_file_name
        open_out = run_cmd(open_cmd)
        if open_out.rstrip('\n') != '0':
          logging.info("Output file is probably corrupted, deleting file and resubmitting job")
        else:
          ret_value = True
      root_tfile.Close()

  # no other checks were performed -> assume that the file is OK
  if min_file_size <= 0 and not validate_outputs:
    ret_value = True

  if not ret_value:
    run_cmd(command)

  return ret_value

def generate_sbatch_line(
    executable,
    command_line_parameter,
    input_file_names,
    output_file_name,
    script_file_name,
    log_file_name     = None,
    keep_logs         = False,
    cvmfs_error_log   = None,
    min_file_size     = MIN_FILE_SIZE_DEFAULT,
    job_template_file = 'sbatch-node.sh.template',
    validate_outputs  = True,
    copy_output_file = True,
  ):
    #------------------------------------------------------------------------------------------------
    # CV: This check (which considerably slows down the execution of 'python tthAnalyzeRun_XXX.py') is no longer needed,
    #     as the template scripts in https://github.com/HEP-KBFI/tth-htt/tree/master/python/templates
    #     perform this check before starting to process the batch job
    #     and skip the batch job in case the output file already exists.
    #
    #if is_file_ok(output_file_name, validate_outputs, min_file_size):
    #  return None
    #------------------------------------------------------------------------------------------------

    if log_file_name and os.path.exists(log_file_name):
        time           = None
        hostname       = None
        is_cvmfs_error = False
        with open(log_file_name, 'r') as log_file:
          for line in log_file:
              if line.find("Time") != -1:
                  time = line.split(':')[1].strip()
              if line.find("Hostname") != -1:
                  hostname = line.split(':')[1].strip()
              if line.find("Transport endpoint is not connected") != -1:
                  is_cvmfs_error = True
        if is_cvmfs_error:
            logging.error("Problem with cvmfs access reported in log file = '%s':" % log_file_name)
            logging.error(" host = '%s': time = %s" % (hostname, time))
            if cvmfs_error_log:
                if not hostname in cvmfs_error_log.keys():
                    cvmfs_error_log[hostname] = []
                cvmfs_error_log[hostname].append(time)

    if type(input_file_names) is str:
        input_file_names = [ input_file_names ]

    input_file_names_duplicates = find_duplicates(input_file_names)
    if input_file_names_duplicates:
      raise RuntimeError(
        "Found duplicate input files to produce output file %s: %s" % \
        (output_file_name, ", ".join(input_file_names_duplicates))
      )

    submissionStatement = "m.submitJob(\n"                            \
      "  inputFiles             = {input_file_names},\n"              \
      "  executable             = '{executable}',\n"                  \
      "  command_line_parameter = '{command_line_parameter}',\n"      \
      "  outputFilePath         = '{output_file_name}',\n"            \
      "  outputFiles            = {output_file_basename},\n"          \
      "  scriptFile             = '{script_file_name}',\n"            \
      "  logFile                = '{log_file_name}',\n"               \
      "  keep_logs              = {keep_logs},\n"                     \
      "  skipIfOutputFileExists = {skipIfOutputFileExists},\n"        \
      "  job_template_file      = '{job_template_file}',\n"           \
      "  copy_output_file       = {copy_output_file},\n"              \
      "  validate_output        = {validate_outputs},\n"              \
      ")".format(
        input_file_names       = input_file_names,
        executable             = executable,
        command_line_parameter = command_line_parameter,
        output_file_name       = os.path.dirname(output_file_name) if output_file_name else '',
        output_file_basename   = [ os.path.basename(output_file_name) ] if output_file_name else [],
        script_file_name       = script_file_name,
        log_file_name          = log_file_name,
        keep_logs              = keep_logs,
        skipIfOutputFileExists = False,
        job_template_file      = job_template_file,
        copy_output_file       = copy_output_file,
        validate_outputs       = validate_outputs,
    )
    return submissionStatement

def get_num_jobs(nof_input_files, max_input_files_per_job):
    assert(nof_input_files > 0)
    assert(max_input_files_per_job > 0)
    nof_jobs = [ int(math.ceil(float(nof_input_files) / max_input_files_per_job)) ]
    while nof_jobs[-1] > 1:
        nof_jobs.append(int(math.ceil(float(nof_jobs[-1]) / max_input_files_per_job)))
    return sum(nof_jobs)

def createScript_sbatch_hadd_nonBlocking(
        sbatch_script_file_name,
        inputFiles,
        outputFiles,
        script_file_name,
        log_file_name           = None,
        working_dir             = None,
        auxDirName              = '',
        verbose                 = False,
        dry_run                 = False,
        max_input_files_per_job = MAX_INPUTS_PER_JOB_DEFAULT,
        use_home                = False,
        min_file_size           = MIN_FILE_SIZE_DEFAULT,
        max_num_submittedJobs   = MAX_NOF_RUNNING_JOBS_DEFAULT,
        max_mem                 = '',
        validate_output         = True,
      ):

    header = """
from tthAnalysis.HiggsToTauTau.sbatchManager import sbatchManager
from tthAnalysis.HiggsToTauTau.ClusterHistogramAggregatorNonBlocking import ClusterHistogramAggregatorNonBlocking

import time

sbatch_managers = []
cluster_histogram_aggregators = []

    """
    job_template = """
m_{{idx}} = sbatchManager(
  '{{pool_id}}', 
  verbose               = {{verbose}}, 
  dry_run               = {{dry_run}}, 
  use_home              = {{use_home}}, 
  min_file_size         = {{min_file_size}},
  max_num_submittedJobs = {{max_num_submittedJobs}},
)
m_{{idx}}.setWorkingDir('{{working_dir}}')
m_{{idx}}.log_completion = {{verbose}}
{% if max_mem|length %}m_{{idx}}.max_mem = '{{max_mem}}'{% endif %}
sbatch_managers.append(m_{{idx}})

cluster_histogram_aggregator_{{ idx }} = ClusterHistogramAggregatorNonBlocking(
  input_files             = {{input_file_names}},
  final_output_file       = '{{output_file_name}}',
  max_input_files_per_job = {{max_input_files_per_job}},
  sbatch_manager          = m_{{ idx }},
  auxDirName              = '{{auxDirName}}',
  script_file_name        = '{{script_file_name}}',
  log_file_name           = '{{log_file_name}}',
  validate_output         = {{validate_output}},
)
cluster_histogram_aggregator_{{idx}}.create_jobs()
cluster_histogram_aggregators.append(cluster_histogram_aggregator_{{idx}})
"""

    footer = """
aggregator_status = { jobIdx : False for jobIdx in range(len(cluster_histogram_aggregators)) }
while True:
  are_all_finished = True
  for jobIdx in aggregator_status:
    if aggregator_status[jobIdx]:
      continue
    else:
      aggregator_status[jobIdx] = cluster_histogram_aggregators[jobIdx].is_done()
      if not aggregator_status[jobIdx]:
        are_all_finished = False
  if are_all_finished:
    break
  else:
    time.sleep(60)
""" 
    script_str = "ClusterHistogramAggregator"
    
    content = []
    nof_jobs = 0
    for idxKey, key in enumerate(outputFiles.keys()):
        input_file_names_duplicates = find_duplicates(inputFiles[key])
        if input_file_names_duplicates:
          raise RuntimeError(
            "Found duplicate input files to produce output file %s: %s" % \
            (outputFiles[key], ", ".join(input_file_names_duplicates))
          )
        assert(script_file_name.find(script_str) != -1)
        script_file_name_key = script_file_name.replace(script_str, "{}_{}".format(key, script_str))
        assert (log_file_name.find(script_str) != -1)
        log_file_name_key = log_file_name.replace(script_str, "{}_{}".format(key, script_str))
        template_vars = {
            'working_dir'             : working_dir,
            'input_file_names'        : inputFiles[key],
            'output_file_name'        : outputFiles[key],
            'auxDirName'              : auxDirName,
            'pool_id'                 : uuid.uuid4(),
            'max_input_files_per_job' : max_input_files_per_job,
            'script_file_name'        : script_file_name_key,
            'log_file_name'           : log_file_name_key,
            'verbose'                 : verbose,
            'dry_run'                 : dry_run,
            'use_home'                : use_home,
            'min_file_size'           : min_file_size,
            'max_num_submittedJobs'   : max_num_submittedJobs,
            'idx'                     : idxKey,
            'max_mem'                 : max_mem,
            'validate_output'         : validate_output,
        }
        job_code = jinja2.Template(job_template).render(**template_vars)
        content.append(job_code)
        nof_jobs += get_num_jobs(len(inputFiles[key]), max_input_files_per_job)

    script_content = header + '\n'.join(content) + footer
    createFile(sbatch_script_file_name, script_content.splitlines())

    return nof_jobs

def createScript_sbatch_hadd(
    sbatch_script_file_name,
    input_file_names,
    output_file_name,
    script_file_name,
    log_file_name           = None,
    working_dir             = None,
    waitForJobs             = True,
    auxDirName              = '',
    pool_id                 = '',
    verbose                 = False,
    dry_run                 = False,
    max_input_files_per_job = MAX_INPUTS_PER_JOB_DEFAULT,
    use_home                = False,
    min_file_size           = MIN_FILE_SIZE_DEFAULT,
    max_num_submittedJobs   = MAX_NOF_RUNNING_JOBS_DEFAULT,
  ):
    """Creates the python script necessary to submit 'hadd' jobs to the batch system
    """
    if not working_dir:
        working_dir = os.getcwd()
    if not pool_id:
        raise ValueError('pool_id is empty')
    sbatch_hadd_lines, num_jobs = generate_sbatch_lines_hadd(
        input_file_names        = input_file_names,
        output_file_name        = output_file_name,
        script_file_name        = script_file_name,
        log_file_name           = log_file_name,
        working_dir             = working_dir,
        waitForJobs             = waitForJobs,
        auxDirName              = auxDirName,
        pool_id                 = pool_id,
        verbose                 = verbose,
        dry_run                 = dry_run,
        max_input_files_per_job = max_input_files_per_job,
        use_home                = use_home,
        min_file_size           = min_file_size,
        max_num_submittedJobs   = max_num_submittedJobs,
    )
    createFile(sbatch_script_file_name, sbatch_hadd_lines)
    return num_jobs

def generate_sbatch_lines_hadd(
    input_file_names,
    output_file_name,
    script_file_name,
    log_file_name,
    working_dir,
    waitForJobs             = True,
    auxDirName              = '',
    pool_id                 = '',
    verbose                 = False,
    dry_run                 = False,
    max_input_files_per_job = MAX_INPUTS_PER_JOB_DEFAULT,
    use_home                = False,
    min_file_size           = MIN_FILE_SIZE_DEFAULT,
    max_num_submittedJobs   = MAX_NOF_RUNNING_JOBS_DEFAULT,
  ):
    template_vars = {
        'working_dir'             : working_dir,
        'input_file_names'        : input_file_names,
        'output_file_name'        : output_file_name,
        'waitForJobs'             : waitForJobs,
        'auxDirName'              : auxDirName,
        'pool_id'                 : pool_id,
        'max_input_files_per_job' : max_input_files_per_job,
        'script_file_name'        : script_file_name,
        'log_file_name'           : log_file_name,
        'verbose'                 : verbose,
        'dry_run'                 : dry_run,
        'use_home'                : use_home,
        'min_file_size'           : min_file_size,
        'max_num_submittedJobs'   : max_num_submittedJobs,
    }
    if not pool_id:
        raise ValueError('pool_id is empty')

    input_file_names_duplicates = find_duplicates(input_file_names)
    if input_file_names_duplicates:
      raise RuntimeError(
        "Found duplicate input files to produce output file %s: %s" % \
        (output_file_name, ", ".join(input_file_names_duplicates))
      )

    sbatch_template = """
from tthAnalysis.HiggsToTauTau.sbatchManager import sbatchManager
from tthAnalysis.HiggsToTauTau.ClusterHistogramAggregator import ClusterHistogramAggregator

m = sbatchManager(
  '{{pool_id}}', 
  verbose               = {{verbose}}, 
  dry_run               = {{dry_run}}, 
  use_home              = {{use_home}}, 
  min_file_size         = {{min_file_size}},
  max_num_submittedJobs = {{max_num_submittedJobs}},
)
m.setWorkingDir('{{working_dir}}')
m.log_completion = {{verbose}}

cluster_histogram_aggregator = ClusterHistogramAggregator(
  input_files             = {{input_file_names}},
  final_output_file       = '{{output_file_name}}',
  max_input_files_per_job = {{max_input_files_per_job}},
  waitForJobs             = {{waitForJobs}},
  sbatch_manager          = m,
  auxDirName              = '{{auxDirName}}',
  script_file_name        = '{{script_file_name}}',
  log_file_name           = '{{log_file_name}}',
)
cluster_histogram_aggregator.create_output_file()
"""
    sbatch_code = jinja2.Template(sbatch_template).render(**template_vars)
    num_jobs = 1

    return sbatch_code.splitlines(), num_jobs
