import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_100_1_J9T.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_10_1_rI8.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_11_1_Y4D.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_12_1_7kH.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_13_1_TBr.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_14_1_BrI.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_15_1_mQ0.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_16_1_hKr.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_17_1_Wuh.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_18_1_Jea.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_19_1_LPm.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_1_1_D2A.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_20_1_q1e.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_21_1_bqk.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_22_1_QM8.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_23_1_lvn.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_24_1_wmh.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_25_1_A68.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_26_1_qbA.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_27_1_STT.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_28_1_DoF.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_29_1_SoA.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_2_1_rqd.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_30_1_XoG.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_31_1_wN7.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_32_1_cBZ.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_33_1_KE7.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_34_1_kEi.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_35_1_LWB.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_36_1_FjM.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_37_1_qqT.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_38_1_RYh.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_39_1_X5d.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_3_1_asX.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_40_1_H2b.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_41_1_UCs.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_42_1_JHh.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_43_1_QEv.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_44_1_WhW.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_45_1_Qn0.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_46_1_TBI.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_47_1_eAW.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_48_1_MvU.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_49_1_dSc.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_4_1_iMV.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_50_1_sQn.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_51_1_fJ5.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_52_1_Pyf.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_53_1_qJV.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_54_1_uc4.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_55_1_AsJ.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_56_1_MMJ.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_57_1_POd.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_58_1_HOe.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_59_1_vkH.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_5_1_qXR.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_60_1_1QB.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_61_1_T3q.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_62_1_C6y.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_63_1_Ksv.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_64_1_fFu.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_65_1_w2L.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_66_1_O4g.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_67_1_Pi4.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_68_1_BY6.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_69_1_8rG.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_6_1_0gn.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_70_1_qd2.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_71_1_E9a.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_72_1_lXG.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_73_1_kJx.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_74_1_kJs.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_75_1_6Ou.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_76_1_iLk.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_77_1_Mmd.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_78_1_gox.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_79_1_dby.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_7_1_0fu.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_80_1_puy.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_81_1_oJd.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_82_1_wPn.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_83_1_Ute.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_84_1_cxA.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_85_1_6ER.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_86_1_BBa.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_87_1_FIT.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_88_1_ulS.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_89_1_j0w.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_8_1_SB1.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_90_1_Tpb.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_91_1_a8L.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_92_1_y50.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_93_1_gO9.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_94_1_7aT.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_95_1_bki.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_96_1_hf9.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_97_1_1x4.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_98_1_EvH.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_99_1_Xre.root',
       '/store/user/zuranski/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/MH_120_MFF_50_CTau500_7TeVGEN_SIM_RAWDEBUG/7995b9c7916ba83563340a6bca2f4caf/GEN-SIM-RAWDEBUG_9_1_HFV.root' ] );


secFiles.extend( [
               ] )


