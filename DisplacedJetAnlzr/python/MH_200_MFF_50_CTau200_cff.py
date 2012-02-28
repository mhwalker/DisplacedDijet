import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_100_1_v6y.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_101_1_CCW.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_10_1_ccx.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_11_1_cwB.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_12_1_GoF.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_13_1_nvd.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_14_1_A0q.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_15_1_mnM.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_16_1_VYA.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_17_1_6bw.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_18_1_Yq2.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_19_1_4p9.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_1_1_1JZ.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_20_1_Vdj.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_21_1_22E.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_22_1_Zj6.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_23_1_tJk.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_24_1_ZRd.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_25_1_bSO.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_26_1_liA.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_27_1_HzA.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_28_1_5x4.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_29_1_4Ch.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_2_1_9RL.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_30_1_7Ha.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_31_1_qnn.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_32_1_lJj.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_33_1_SGv.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_34_1_J10.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_35_1_AVT.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_36_1_3Nq.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_37_1_aLV.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_38_1_nkW.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_39_1_end.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_3_1_2Xs.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_40_1_sXx.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_41_1_lja.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_42_1_YSG.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_43_1_SQl.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_44_1_ZIl.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_45_1_rF0.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_46_1_9Rp.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_47_1_ffr.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_48_1_HfL.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_49_1_Jyz.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_4_1_iC7.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_50_1_k79.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_51_1_vIs.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_52_1_upl.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_53_1_BxM.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_54_1_s1i.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_55_1_ZFG.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_56_1_GuW.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_57_1_dX5.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_58_1_IEs.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_59_1_qF0.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_5_1_ht8.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_60_1_Gw6.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_61_1_0yG.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_62_1_Y9N.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_63_1_EAo.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_64_1_PDv.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_65_1_CMp.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_66_1_nMd.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_67_1_sz9.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_68_1_07c.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_69_1_fnz.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_6_1_c6r.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_70_1_Rkv.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_71_1_9uh.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_72_1_usQ.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_73_1_Fog.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_74_1_07R.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_75_1_gqY.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_76_1_yvB.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_77_1_ze3.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_78_1_rEG.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_79_1_Qt1.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_7_1_Sqb.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_80_1_bog.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_81_1_vi0.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_82_1_50H.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_83_1_b21.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_84_1_mAH.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_85_1_RPx.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_86_1_kGU.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_87_1_LZy.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_88_1_8XB.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_89_1_zbb.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_8_1_dDC.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_90_1_8Z6.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_91_1_eIw.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_92_1_Nf4.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_93_1_dh9.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_94_1_tFQ.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_95_1_Gmn.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_96_1_gQ4.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_97_1_lWF.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_98_1_fWG.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_99_1_QMd.root',
       '/store/user/zuranski/MH_200_MFF_50_CTau200_7TeVGEN_SIM_RAWDEBUG/MH_200_MFF_50_CTau200_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_9_1_0IW.root' ] );


secFiles.extend( [
               ] )

