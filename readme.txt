To run everything we run: severe_paper_runner.R

This essentially runs steps 1-4 below

We run the following files (in this order) 
1) To initialise: init-decay.R 
   requires 
      WP4_modelImplementationFunctions.R
      severe_paper_helper_functions.R
      import_decay_data.R
      plot_decay_efficacy.R
      natmed_modelImplementationFunctions.R
   and data files:
      DECAY_Master_ForPaper.xlsx
      natmed_parameters.RData
      natmed_parameters_sympt_and_severe.RData
      BoostingTableSummary_20210726,csv
2) Generate imulations: simulateProgressionsUpdated.R 
   requires: 
      load_bootstrap_parameter_distributions.R (this is where the number of simulations, nruns, is set, currently nruns = 1e5)
      generateSimNatMedCurve.R
   in addition requires th following data file s(unless they are gettign regenerated:
      non_timecourse_paramseter_structures_n1e+05.RData
      simulated_natmed_efficacy_curve_n1e+05.RData
      time_course_summary_n1e+05.RData
3) Load and make the data required for the plots: loadDataForPlots.R 
  requires:
      plot_efficacies_revised.R  (which in turn requires ab_efficacy_estimates.R, natmed_parameters.RData and efficacies.RDS)
      load_bootstrap_parameter_distributions.R
4) Do regression (line 156) and make plots that are used in the manuscript: paperPlotsUpdated.R 
  requires:
      load_bootstrap_parameter_distributions.R
5) Generate the CI intervals and statistics for the manuscript: matchCIs_updated2.R