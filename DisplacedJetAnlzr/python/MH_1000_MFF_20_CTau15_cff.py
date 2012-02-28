import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_100_1_0M3.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_101_1_KWB.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_102_1_PmN.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_10_1_Idd.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_11_1_eZn.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_13_1_sd6.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_14_1_KiD.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_15_1_Enx.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_16_1_DRA.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_17_1_XI9.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_18_1_yrx.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_19_1_zGk.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_1_1_CEJ.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_20_1_Zfz.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_21_1_eWP.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_22_1_bqH.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_23_1_0C6.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_24_1_8Wu.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_25_1_7ah.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_26_1_Ps5.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_27_1_Fzm.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_28_1_VuO.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_29_1_AJr.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_2_1_Vrv.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_30_1_zrR.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_31_1_ZeR.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_32_1_Kah.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_33_1_xoJ.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_34_1_X8w.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_35_1_6Yr.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_36_1_iIF.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_37_1_xsW.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_38_1_gT1.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_39_1_PmN.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_3_1_iYg.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_40_1_jCN.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_41_1_6NB.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_42_1_u1E.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_43_1_ir2.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_44_1_xnG.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_45_1_YWo.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_46_1_iC2.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_47_1_nfg.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_48_1_9UY.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_49_1_snQ.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_4_1_J5O.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_50_1_CZ9.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_51_1_rkb.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_52_1_bgV.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_53_1_3T0.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_54_1_w2h.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_55_1_1Y8.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_56_1_mT6.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_57_1_0hG.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_58_1_pn4.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_59_1_uOt.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_5_1_HSx.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_60_1_cBH.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_61_1_IOW.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_62_1_yJT.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_63_1_acs.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_64_1_8Lo.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_65_1_UD5.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_66_1_HmU.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_67_1_UUO.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_68_1_Xiq.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_69_1_SyF.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_6_1_vLO.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_70_1_ydY.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_71_1_A6g.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_72_1_IhT.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_73_1_O6Y.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_74_1_DWc.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_75_1_0zF.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_76_1_R1t.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_77_1_GkA.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_78_1_L4r.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_79_1_dsM.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_7_1_Pkm.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_80_1_dwy.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_81_1_WOQ.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_82_1_Hdx.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_83_1_fFM.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_84_1_rLg.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_85_1_gsC.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_86_1_Y3D.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_87_1_agE.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_88_1_J1w.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_89_1_U0a.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_8_1_XUB.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_90_1_dgm.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_91_1_Ojn.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_92_1_UHQ.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_93_1_R70.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_94_1_v3q.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_95_1_lsb.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_96_1_cTg.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_97_1_dpb.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_98_1_lZP.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_99_1_AqC.root',
       '/store/user/zuranski/MH_1000_MFF_20_CTau15_7TeVGEN_SIM_RAWDEBUG/MH_1000_MFF_20_CTau15_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_9_1_cKS.root' ] );


secFiles.extend( [
               ] )

