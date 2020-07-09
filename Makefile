
DATPATH=./DATA
FUNPATH=./PROGRAMS
RESPATH=./RESULTS

all: data figures simulations

data: $(DATPATH)/master_IPD.Rdata \
			$(DATPATH)/analysis_dset.Rdata 

simulations: $(DATPATH)/fit_pem_type.Rdata \
						 $(DATPATH)/sens_TC.Rdata \
						 $(DATPATH)/sens_IC.Rdata \
						 $(DATPATH)/sens_notNSCLC.Rdata \
						 $(DATPATH)/sens_NSCLC.Rdata \
						 $(DATPATH)/prior_sensitivity_2.Rdata \
						 $(DATPATH)/weighted_analysis.Rdata  \
						 $(DATPATH)/meta_designs_sims_2.Rdata \
						 $(DATPATH)/sens_not_placebo.Rdata

						 
figures: $(RESPATH)/Figure_1.pdf \
				 $(RESPATH)/Figure_2.pdf \
				 $(RESPATH)/Figure_3.pdf \
				 $(RESPATH)/forestplot_heterogeneity.pdf \
				 $(RESPATH)/flowchart_data.txt \
				 $(RESPATH)/Figure_prior_sensitivity_2.pdf \
				 $(RESPATH)/Figure_weighted_analysis.pdf \
				 $(RESPATH)/results_for_manuscript.txt \
				 $(RESPATH)/forestplot_het_rawests.pdf \
				 $(RESPATH)/meta_designs_power_2.txt \
				 $(RESPATH)/Figure_sens_placebo.pdf

				 
clean:
	rm -fv $(DATPATH)/*.Rdata
	rm -fv $(RESPATH)/*.pdf
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out

cleanlogs:
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out

### DATA

$(DATPATH)/master_IPD.Rdata: $(FUNPATH)/prepare_master_ipd_dataset.R \
                             $(DATPATH)/IPD/*.csv
	R CMD BATCH $(FUNPATH)/prepare_master_ipd_dataset.R

$(DATPATH)/analysis_dset.Rdata: $(FUNPATH)/reconstruct_pdl1_negative_groups.R \
                                $(DATPATH)/master_IPD.Rdata \
                                $(DATPATH)/Publication_level_variables/pubs_variables.csv
	R CMD BATCH $(FUNPATH)/reconstruct_pdl1_negative_groups.R

### SIMULTIONS

# Main analysis
$(DATPATH)/fit_pem_type.Rdata: $(FUNPATH)/functions_collapsed_gibbs_cpp_pem_types.R \
                     					 $(FUNPATH)/fit_collapsed_gibbs_cpp_pem_type.R \
                     					 $(FUNPATH)/cpp_code_pem_tumor_type.cpp \
                     					 $(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/fit_collapsed_gibbs_cpp_pem_type.R

# Stratified analysis
$(DATPATH)/sens_TC.Rdata \
$(DATPATH)/sens_IC.Rdata \
$(DATPATH)/sens_notNSCLC.Rdata \
$(DATPATH)/sens_NSCLC.Rdata: $(FUNPATH)/sensitivity_analysis.R \
														 $(FUNPATH)/functions_collapsed_gibbs_cpp_pem_types.R \
                     				 $(FUNPATH)/cpp_code_pem_tumor_type.cpp \
                     				 $(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/sensitivity_analysis.R														 

# Stratified analysis by assay type
$(DATPATH)/sens_ICH22C3.Rdata \
$(DATPATH)/sens_ICH288.Rdata \
$(DATPATH)/sens_ICH7310.Rdata \
$(DATPATH)/sens_ICHSP142.Rdata: $(FUNPATH)/suppl_strat_by_IHC.R \
																$(FUNPATH)/functions_collapsed_gibbs_cpp_pem_types.R \
                     				 		$(FUNPATH)/cpp_code_pem_tumor_type.cpp \
                     				 		$(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/suppl_strat_by_IHC.R	
 
# Simulation study of different designs
$(DATPATH)/meta_designs_sims_2.Rdata \
$(RESPATH)/meta_designs_power_2.txt: $(FUNPATH)/meta_analytic_designs_simulations_2.R \
																		 $(FUNPATH)/functions_permutation_test.R \
																		 $(DATPATH)/analysis_dset.Rdata \
																	   $(DATPATH)/fit_pem_type.Rdata
	R CMD BATCH $(FUNPATH)/meta_analytic_designs_simulations_2.R
												
# Sensitivity analysis for the prior distributions
$(DATPATH)/prior_sensitivity_2.Rdata: $(FUNPATH)/prior_sensitivity_2.R \
                     					      $(FUNPATH)/cpp_code_pem_tumor_type.cpp \
                     					      $(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/prior_sensitivity_2.R

# Weighted analysis
$(DATPATH)/weighted_analysis.Rdata: $(FUNPATH)/weighted_analysis.R \
																		$(DATPATH)/analysis_dset.Rdata \
																		$(FUNPATH)/functions_collapsed_gibbs_cpp_pem_types.R \
																		$(FUNPATH)/cpp_code_pem_tumor_type.cpp
	R CMD BATCH $(FUNPATH)/weighted_analysis.R

# Sensitivity analysis excluding placebo-controlled trials
$(DATPATH)/sens_not_placebo.Rdata: $(FUNPATH)/sensitivity_analysis_placebo.R \
                     					      $(FUNPATH)/cpp_code_pem_tumor_type.cpp \
                     					      $(FUNPATH)/functions_collapsed_gibbs_cpp_pem_types.R \
                     					      $(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/sensitivity_analysis_placebo.R
	
### FIGURES AND TABLES

$(RESPATH)/Figure_1.pdf: $(FUNPATH)/Figure1.R \
												 $(FUNPATH)/functions_piecewise_exponential.R \
												 $(DATPATH)/analysis_dset.Rdata \
												 $(DATPATH)/fit_pem_type.Rdata
	R CMD BATCH $(FUNPATH)/Figure1.R

$(RESPATH)/Figure_2.pdf: $(FUNPATH)/Figure2.R \
																	 $(DATPATH)/sens_TC.Rdata \
																	 $(DATPATH)/sens_IC.Rdata \
																	 $(DATPATH)/sens_notNSCLC.Rdata \
																	 $(DATPATH)/sens_NSCLC.Rdata \
																	 $(DATPATH)/sens_ICH22C3.Rdata \
																	 $(DATPATH)/sens_ICH288.Rdata \
																	 $(DATPATH)/sens_ICH7310.Rdata \
																	 $(DATPATH)/sens_ICHSP142.Rdata
	R CMD BATCH $(FUNPATH)/Figure2.R	

$(RESPATH)/Figure_3.pdf: $(FUNPATH)/Figure3.R \
												 $(FUNPATH)/functions_piecewise_exponential.R \
												 $(DATPATH)/analysis_dset.Rdata \
												 $(DATPATH)/fit_pem_type.Rdata
	R CMD BATCH $(FUNPATH)/Figure3.R

$(RESPATH)/forestplot_heterogeneity.pdf: $(FUNPATH)/heterogeneity.R \
																				 $(DATPATH)/analysis_dset.Rdata \
																				 $(DATPATH)/fit_pem_type.Rdata
	R CMD BATCH $(FUNPATH)/heterogeneity.R

$(RESPATH)/flowchart_data.txt: $(DATPATH)/Flowchart_data/IncludedTrials.csv \
                               $(FUNPATH)/flowchart.R
	R CMD BATCH $(FUNPATH)/flowchart.R

$(RESPATH)/Figure_prior_sensitivity_2.pdf: $(FUNPATH)/Figure_prior_sensitivity_2.R \
												 $(FUNPATH)/functions_piecewise_exponential.R \
												 $(DATPATH)/analysis_dset.Rdata \
												 $(DATPATH)/prior_sensitivity_2.Rdata
	R CMD BATCH $(FUNPATH)/Figure_prior_sensitivity_2.R


$(RESPATH)/Figure_weighted_analysis.pdf: $(FUNPATH)/Figure_weighted_analysis.R \
																				 $(FUNPATH)/functions_piecewise_exponential.R \
												 								 $(DATPATH)/analysis_dset.Rdata \
												 								 $(DATPATH)/weighted_analysis.Rdata
	R CMD BATCH $(FUNPATH)/Figure_weighted_analysis.R												 								 

$(RESPATH)/results_for_manuscript.txt: $(FUNPATH)/results_for_manuscript.R \
																			 $(DATPATH)/fit_pem_type.Rdata \
																			 $(DATPATH)/analysis_dset.Rdata \
																			 $(FUNPATH)/functions_piecewise_exponential.R
	R CMD BATCH $(FUNPATH)/results_for_manuscript.R 

$(RESPATH)/forestplot_het_rawests.pdf: $(FUNPATH)/meta_analytic_estimates.R \
																			 $(DATPATH)/analysis_dset.Rdata
	R CMD BATCH $(FUNPATH)/meta_analytic_estimates.R
	
$(RESPATH)/Figure_sens_placebo.pdf: $(FUNPATH)/Figure_sensitivity_analysis_placebo.R \
																		$(DATPATH)/sens_not_placebo.Rdata
	R CMD BATCH $(FUNPATH)/Figure_sensitivity_analysis_placebo.R
