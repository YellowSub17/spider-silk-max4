



PROP_ID = 20251747
# DATA_PATH = f'/data/visitors/cosaxs/{PROP_ID}/2026052808'
DATA_PATH = f'/home/pat/max4/spider-silk-max4/data'

PHOTON_ENERGY = 12400.0 #keV
PHOTON_LAMBDA = 9.998725680096795e-11


EIGER_QBANDS = [ (0.0035,0.008), (0.009,0.06), (0.07,0.2) ]
PILATUS_QBANDS = [ (0.15,0.25), (0.3,0.6), (0.75,1.1), (1.2,1.4), (1.6,1.9) ]




Y2F_buffer8_runs = [i for i in range(102532, 102628)]
Y2F_aa_runs = [i for i in range(102628, 102726)]
buffer_aa_runs = [i for i in range(102726, 102809)]
capillary1_runs = [i for i in range(102809,102821)]
capillary2_runs = [i for i in range(102821, 102833)]
capillary3_runs = [i for i in range(102833, 102852)]
Y2F_h2o_h2o_runs = [i for i in range(102901, 103021)]
Y2F_h2o_kpi_1_runs = [i for i in range(103021, 103091)]
Y2F_h2o_kpi_2_runs = [i for i in range(103091, 103171)]
Y2F_h2o_kpi_nopump_runs = [i for i in range(103171, 103186)]
kpi_h2o_1_runs = [i for i in range(103186, 103231)]
kpi_h2o_2_runs = [i for i in range(103231, 103326)]
kpi_h2o_3_runs = [i for i in range(103326,103331)]
YR2A_h2o_kpi_runs = [i for i in range(103331, 103376)]
YR2A_h2o_h2o_runs = [i for i in range(103376, 103414)]
Y2F_kpi_h2o_runs = [i for i in range(103414, 103456)]
Y2F_kpi_AcA_1_runs = [i for i in range(103456, 103526)]
Y2F_kpi_AcA_2_runs = [i for i in range(103526, 103609)]
resin_runs = [i for i in range(103609, 103624)]
