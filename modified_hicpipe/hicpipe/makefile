
#####################################################################################################
# Utility gmake canned recipe
#####################################################################################################

# checks return code of last shell run
define ASSERT
if [ $$? -ne 0 ]; then exit 1; fi
endef

# user defined parameters are in the cfg file
include $(CFG)

######################################################################################################
# Parameters: mapping coordinates to fragment ends
######################################################################################################

# discard read of sum of distance to first cutter site is above this threshold
SEGMENT_LEN_THRESHOLD?=500

######################################################################################################
# Parameters: experimental bias modelling
######################################################################################################

# file describing normaliztion model
MFN?=models/map_len_gc.mdl

# estimate biases using (trans/far_cis/close_cis/both) contacts
BIAS_TYPE?=trans

# threshold between close and far cis definitions, relevant if BIAS_TYPE equals far_cis or close_cis
CIS_THRESHOLD?=1000000

######################################################################################################
# Parameters: output contact maps
######################################################################################################

# maximal simultaneous jobs running on sge cluster
CLUSTER_MAXJOBS?=200

# The following three vectors must be of the same length N>=1. They define:
# OUT_BINSIZE (bp): Squared binsize
# OUT_FILTER (both/trans/close_cis/far_cis): All vs. all or limit to bins close to cis diagonal
# OUT_DIAGONAL_DIST (bp): threshold between close and far cis, relevant for FILTER=close_cis/far_cis
# OUT_N: must be set to 1,2,..,N
OUT_BINSIZES?=500000 1000000 2000000
OUT_FILTERS?=close_cis close_cis both
OUT_DIAGONAL_DISTS?=10000000 20000000 0
OUT_N?=1 2 3

# window size used when smoothing matrix
SMOOTH_WIDTH?=10

######################################################################################################
# Parameters: R paths
######################################################################################################

# path to the Rscript binary
RSCRIPT?=Rscript

######################################################################################################
# Main pipeline: 
# - Check input files
# - Build regular and smooth contact matrices
######################################################################################################

all: check_user_params init
# print parameters
	@echo Input fragment ends table: $(FENDS_TABLE)
	@echo Input paired read files: $(READ_FILES)
	@echo Output directory: $(ODIR)
	@echo Dataset title: $(DATASET)
	@echo Star activity threshold \(bp\): $(SEGMENT_LEN_THRESHOLD)
	@echo Correction model file: $(MFN)
	@echo Correction bias type: $(BIAS_TYPE)
	@echo close/far cis threshold \(relevant only for cis biases\): $(CIS_THRESHOLD)
	@echo Binsizes: $(OUT_BINSIZES)
	@echo Smoothing window size \(in bins\): $(SMOOTH_WIDTH)
	@if [ ${CLUSTER} == 1 ]; then echo Distributing over Sun Grid Engine; else echo Running sequentially on local machine; fi
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	$(foreach I, $(OUT_N), make contact_map \
	                            CFG=$(CFG) \
	                            BINSIZE=$(word $(I), $(OUT_BINSIZES)) \
	                            FILTER=$(word $(I), $(OUT_FILTERS)) \
	                            DIAGONAL_DIST=$(word $(I), $(OUT_DIAGONAL_DISTS)); $(ASSERT); )

check_user_params:
# check user parameters
	@if [ -z ${CFG} ]; then echo The config file variable CFG is not defined, see README; exit 1; fi
	@if [ -z ${FENDS_TABLE} ]; then echo FENDS_TABLE not defined in config file ${CFG}, see README; exit 1; fi
	@if [ -z "${READ_FILES}" ]; then echo READ_FILES not definedin config file ${CFG} , see README; exit 1; fi
	@if [ -z ${DATASET} ]; then echo DATASET not defined in config file ${CFG}, see README; exit 1; fi
	@if [ -z ${ODIR} ]; then echo ODIR not defined in config file ${CFG}, see README; exit 1; fi
# check input files
	@if [ ! -e ${FENDS_TABLE} ]; then echo Input fends table ${FENDS_TABLE} does not exist; exit 1; fi
	@for IFILE in ${READ_FILES}; do { if [ ! -e $$IFILE ]; then echo Input paired read file $$IFILE does not exist; exit 1; fi }; done

# remove output directory
clean: check_user_params
	rm -rf $(ODIR)

#####################################################################################################
# Binaries
#####################################################################################################

bin/smooth_matrix: src/smooth_matrix.cpp
	@mkdir -p bin
	g++ $^ -O2 -o $@ -Wall
bin/model_integrate: src/model_integrate.cpp
	@mkdir -p bin
	g++ $^ -O2 -o $@ -Wall
init: bin/smooth_matrix bin/model_integrate

.DELETE_ON_ERROR:
.SECONDARY:

######################################################################################################
# PIPELINE
######################################################################################################

######################################################################################################
# Step 0: dataset initializtion
#
# create link from map_bin.f and map.bin_ranges into results dir
# extract limited fends table from full fends table
#
######################################################################################################

$(ODIR)/%.fends:
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 0: initialization, $*

	mkdir -p $(ODIR)
	mkdir -p tmp
	mkdir -p log
	ln -sf $(CURDIR)/const_correction/map_bin.f $(ODIR)/$*_map_bin.f
	ln -sf $(CURDIR)/const_correction/map.bin_ranges $(ODIR)/$*_map.bin_ranges
	head -n 1 ${FENDS_TABLE} > $@
	cat ${FENDS_TABLE} | awk '$$5 == 1' >> $@

######################################################################################################
# Step 1 : mapped paired reads to paired fends
# generate mat file (paired fends) from the paired reads
######################################################################################################

$(ODIR)/%.mat: $(ODIR)/%.fends
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 1: reads to fend pairs, $*
	./lscripts/coords2fends.pl ${FENDS_TABLE} $(ODIR)/$* 0 $(SEGMENT_LEN_THRESHOLD) $(ODIR)/$*.random $(READ_FILES)

######################################################################################################
# Step 2: prepare nm file
# create the nm file from the mat and fend tables
######################################################################################################

$(ODIR)/%.nm :$(ODIR)/%.fends $(ODIR)/%.mat
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 2: compute nm table, $*
	$(RSCRIPT) R/model_preprocess.r $(ODIR)/$* fends $(MFN) 3 $(ODIR)/$*.nm $(BIAS_TYPE) $(CIS_THRESHOLD) $(CLUSTER) $(CLUSTER_MAXJOBS) $(RSCRIPT)

######################################################################################################
# Step 3: learn model
# Create model biases (.f files) from the nm file
######################################################################################################

$(ODIR)/%.model: $(ODIR)/%.nm
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 3: learn model, $*
	$(RSCRIPT) R/learn_model.r $(ODIR)/$* $* $(MFN)
	touch $@

#####################################################################################################
# Step 4: Bin by coord
# create cbinned file from the fends table
#####################################################################################################

# add binsize to dataset name
$(ODIR)/%_$(BINSIZE).cbinned $(ODIR)/%_$(BINSIZE).cbins: $(ODIR)/%.model
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 4: spatial binning, $*
	ln -sf $*.mat $(ODIR)/$*_$(BINSIZE).mat
	ln -sf $*_frag_len_bin.f $(ODIR)/$*_$(BINSIZE)_frag_len_bin.f
	ln -sf $*_frag_gc_bin.f $(ODIR)/$*_$(BINSIZE)_frag_gc_bin.f
	ln -sf $*_map_bin.f $(ODIR)/$*_$(BINSIZE)_map_bin.f
	ln -sf $*.prior $(ODIR)/$*_$(BINSIZE).prior
	./lscripts/bin_coords.pl $(ODIR)/$*.binned $(ODIR)/$*_$(BINSIZE).cbinned $(ODIR)/$*_$(BINSIZE).cbins $(BINSIZE)

#####################################################################################################
# Step 5: compute observed and expected matrices, smooth them and unite results
# create various matrices (observed, expected, united)
#####################################################################################################

$(ODIR)/%.o_contact : $(ODIR)/%.cbinned
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 5a: observed matrix, $*
	./lscripts/observed_cbin_counts.pl $(ODIR)/$* $(DIAGONAL_DIST) $@

$(ODIR)/%.e_contact: $(ODIR)/%.cbinned
	@echo /////////////////////////////////////////////////////////////////////////////////////////
	@echo Step 5b: expected matrix, $*
	$(RSCRIPT) R/model_wrapper.r $(abspath $(ODIR)/$*) cbinned $(abspath $(MFN)) $(FILTER) $(DIAGONAL_DIST) \
                                 $(abspath $@) $(CLUSTER) $(CLUSTER_MAXJOBS) 1 coord_bin 1 coord_bin
	sed -i 's/value/expected_count/g' $@

$(ODIR)/%.o_smooth : $(ODIR)/%.o_contact
	./bin/smooth_matrix $(ODIR)/$*.cbins $(ODIR)/$*.o_contact coord_bin observed_count $(BINSIZE) $(SMOOTH_WIDTH) F $@
$(ODIR)/%.e_smooth : $(ODIR)/%.e_contact
	./bin/smooth_matrix $(ODIR)/$*.cbins $(ODIR)/$*.e_contact coord_bin expected_count $(BINSIZE) $(SMOOTH_WIDTH) F $@

# append observed field to expected file
$(ODIR)/%.n_smooth: $(ODIR)/%.o_smooth $(ODIR)/%.e_smooth
	./lscripts/append_table.pl $(ODIR)/$*.e_smooth $(ODIR)/$*.o_smooth $@ observed_count 0 T cbin1 cbin2
$(ODIR)/%.n_contact: $(ODIR)/%.o_contact $(ODIR)/%.e_contact
	./lscripts/append_table.pl $(ODIR)/$*.e_contact $(ODIR)/$*.o_contact $@ observed_count 0 T coord_bin1 coord_bin2
	sed -i 's/coord_bin/cbin/g' $@

contact_map: $(ODIR)/$(DATASET)_$(BINSIZE).n_contact
smooth_contact_map: $(ODIR)/$(DATASET)_$(BINSIZE).n_smooth