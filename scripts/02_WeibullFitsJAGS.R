# Load Libraries ----
# Load Data ----
plc_curves <- read.csv("data_clean/CavitronDataClean.csv",row.names = "X")# load clean cavitron data
plc_curves <- plc_curves[,c("Sample_ref_1","Sample_ref_2","Raw_conductance_kg_Mpa_s","Conductivity_SI_corrT","Pressure_Mpa","kmax","Ksmax","PLCclean")]# select columns of interest

# Extract Tree.No and sample Replicate ----
plc_curves$Tree.No <- trunc(plc_curves$Sample_ref_1)
plc_curves$Rep <- plc_curves$Sample_ref_1
plc_curves$Rep <- ifelse(test = (plc_curves$Rep - plc_curves$Tree.No)*10==0, yes = 1, no = (plc_curves$Rep - plc_curves$Tree.No)*10)

# merge data ----
sample_id <- read.csv("data_clean/dbg_cottonwood_plantID.csv")# load sample id information
head(sample_id)


