from tthAnalysis.HiggsToTauTau.configs.addMEMConfig import *
from tthAnalysis.HiggsToTauTau.analysisTools import create_cfg

import os

class addMEMConfig_3l_1tau(addMEMConfig):

  def __init__(self,
        treeName,
        outputDir,
        localDir,
        cfgDir,
        executable_addMEM,
        samples,
        era,
        check_output_files,
        leptonSelection,
        hadTauSelection,
        running_method,
        max_files_per_job,
        mem_integrations_per_job,
        max_mem_integrations,
        num_parallel_jobs,
        isDebug,
        jet_cleaning_by_index,
        central_or_shift,
        dry_run,
        use_nonnominal,
        use_home,
        submission_cmd = None,
        pool_id = '',
      ):
    addMEMConfig.__init__(self,
      treeName                 = treeName,
      outputDir                = outputDir,
      localDir                 = localDir,
      cfgDir                   = cfgDir,
      executable_addMEM        = executable_addMEM,
      samples                  = samples,
      era                      = era,
      check_output_files       = check_output_files,
      running_method           = running_method,
      max_files_per_job        = max_files_per_job,
      mem_integrations_per_job = mem_integrations_per_job,
      max_mem_integrations     = max_mem_integrations,
      num_parallel_jobs        = num_parallel_jobs,
      leptonSelection          = leptonSelection,
      hadTauSelection          = hadTauSelection,
      integration_choice       = '',
      jet_cleaning_by_index    = jet_cleaning_by_index,
      dry_run                  = dry_run,
      use_nonnominal           = use_nonnominal,
      use_home                 = use_home,
      channel                  = "3l_1tau",
      submission_cmd           = submission_cmd,
      pool_id                  = pool_id,
    )

    self.template_dir = os.path.join(
      os.getenv('CMSSW_BASE'), 'src', 'tthAnalysis', 'HiggsToTauTau', 'test', 'templates'
    )
    logging.info("Templates directory is: {templateDir}".format(templateDir = self.template_dir))
    self.cfgFile_addMEM_original = os.path.join(self.template_dir, "addMEM_3l_1tau_cfg.py")
    self.maxPermutations_branchName = "maxPermutations_addMEM_%s_lep%s_tau%s_%s" % (
      self.channel, self.leptonSelection, self.hadTauDefinition, self.hadTauWorkingPoint,
    )
    self.isDebug = isDebug
    self.central_or_shift = central_or_shift

  def createCfg_addMEM(self, inputFiles, startRange, endRange, outputFile, era, process, isMC, cfgFile_modified):
    """Create python configuration file for the addMEM_3l_1tau executable (MEM code)

    Args:
      inputFile: list of input files (Ntuples)
      outputFile: output file of the job
    """

    '''Let's assume that the following configuration options remain constant at all times and need not be modified
    - process.fwliteInput.outputEvery
    - process.addMEM_3l_1tau.branchName_electrons
    - process.addMEM_3l_1tau.branchName_muons
    - process.addMEM_3l_1tau.branchName_hadTaus
    - process.addMEM_3l_1tau.branchName_jets
    - process.addMEM_3l_1tau.branchName_met
    '''

    lines = []
    skipEvents = startRange
    maxEvents = endRange - startRange
    lines.append("process.fwliteInput.fileNames = cms.vstring(%s)" % inputFiles)
    lines.append("process.fwliteInput.skipEvents = cms.uint32(%s)" % skipEvents)
    lines.append("process.fwliteInput.maxEvents = cms.int32(%s)" % maxEvents)
    lines.append("process.fwliteOutput.fileName = cms.string('%s')" % os.path.basename(outputFile))
    lines.append("process.addMEM_3l_1tau.era = cms.string('%s')" % era)
    if skipEvents > 0:
      lines.append("process.addMEM_3l_1tau.copy_histograms = cms.vstring()")
    lines.append("process.addMEM_3l_1tau.leptonSelection = cms.string('%s')" % self.leptonSelection)
    lines.append("process.addMEM_3l_1tau.hadTauSelection = cms.string('%s')" % self.hadTauSelection)
    lines.append("process.addMEM_3l_1tau.isMC = cms.bool(%s)" % isMC)
    lines.append("process.addMEM_3l_1tau.isDEBUG = cms.bool(%s)" % self.isDebug)
    lines.append("process.addMEM_3l_1tau.central_or_shift = cms.vstring(%s)" % self.central_or_shift)
    lines.append("process.addMEM_3l_1tau.dryRun = cms.bool(%s)" % self.dry_run)
    lines.append("process.addMEM_3l_1tau.useNonNominal = cms.bool(%s)" % self.use_nonnominal)
    lines.append("process.addMEM_3l_1tau.jetCleaningByIndex = cms.bool(%s)" % self.jet_cleaning_by_index)

    create_cfg(self.cfgFile_addMEM_original, cfgFile_modified, lines)
