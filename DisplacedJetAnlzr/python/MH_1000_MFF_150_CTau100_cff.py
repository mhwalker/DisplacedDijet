import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_100_1_Qj2.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_101_1_16a.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_102_1_i8h.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_10_1_q0o.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_11_1_NEK.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_12_1_war.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_13_1_CYR.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_14_1_Mig.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_15_1_YT0.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_16_1_cNE.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_17_1_xl8.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_18_1_NoW.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_19_1_6r0.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_1_1_8n6.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_20_1_YRc.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_21_1_ZzS.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_22_1_Zpq.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_23_1_mQ6.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_24_1_RwU.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_25_1_4QI.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_26_1_gAt.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_27_1_nd7.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_28_1_i6K.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_29_1_jH1.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_2_1_fiU.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_30_1_5Am.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_31_1_wTK.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_32_1_bgX.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_33_1_zlv.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_34_1_HxN.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_35_1_zLw.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_36_1_TRb.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_37_1_EyZ.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_38_1_TTu.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_39_1_se7.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_3_1_weH.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_40_1_GNu.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_41_1_PxD.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_42_1_MXJ.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_43_1_QKu.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_44_1_10Q.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_45_1_f3l.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_46_1_nzT.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_47_1_sjl.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_48_1_6gg.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_49_1_rVs.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_4_1_5SA.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_50_1_3mu.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_51_1_f2T.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_52_1_wNH.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_53_1_kzf.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_54_1_eKx.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_55_1_7bI.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_56_1_6PV.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_57_1_7yv.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_58_1_SLn.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_59_1_x7c.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_5_1_RFk.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_60_1_mPH.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_61_1_lhF.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_62_1_Zuq.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_63_1_REZ.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_64_1_axK.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_65_1_n3z.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_66_1_Rdu.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_67_1_xSx.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_68_1_sQa.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_69_1_JTP.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_6_1_Wqm.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_70_1_Wnr.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_71_1_qGl.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_72_1_I3B.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_73_1_Qax.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_74_1_CzL.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_75_1_656.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_76_1_SYa.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_77_1_WYr.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_78_1_0Fq.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_79_1_Fxg.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_7_1_gJO.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_80_1_N78.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_81_1_KMP.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_82_1_GPe.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_83_1_qv6.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_84_1_tQ8.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_85_1_8Dz.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_86_1_Sma.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_87_1_hiY.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_88_1_F6H.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_89_1_tIn.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_8_1_j7d.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_90_1_DzI.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_91_1_80y.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_92_1_fei.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_93_1_f5K.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_94_1_rfg.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_95_1_0sA.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_96_1_mtY.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_97_1_lbd.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_98_1_QZP.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_99_1_P04.root',
       '/store/user/zuranski/MH_1000_MFF_150_CTau100_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_150_CTau100_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_9_1_AYO.root' ] );


secFiles.extend( [
               ] )


