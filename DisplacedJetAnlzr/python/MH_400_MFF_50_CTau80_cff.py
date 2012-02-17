import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_100_1_s3C.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_101_1_Ip3.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_10_1_XVj.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_11_1_5p9.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_12_1_nQu.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_13_1_aYC.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_14_1_0xa.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_15_1_hGL.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_16_1_E3P.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_17_1_79m.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_18_1_lGz.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_19_1_9zc.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_1_1_XEm.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_20_1_RXu.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_21_1_bWQ.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_22_1_Q4U.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_23_1_HVy.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_24_1_nWS.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_25_1_EAX.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_26_1_uMJ.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_27_1_eUy.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_28_1_D0V.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_29_1_KOw.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_2_1_RKu.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_30_1_TAY.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_31_1_7GA.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_32_1_uVh.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_33_1_dej.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_34_1_rbm.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_35_1_eLN.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_36_1_zOo.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_37_1_vrl.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_38_1_tj6.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_39_1_IUb.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_3_1_pXi.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_40_1_fuO.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_41_1_dxH.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_42_1_22R.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_43_1_ata.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_44_1_arB.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_45_1_XRs.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_46_1_fjb.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_47_1_o4x.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_48_1_kBc.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_49_1_BaP.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_4_1_26c.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_50_1_R2r.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_51_1_BcI.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_52_1_NXC.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_53_1_h96.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_54_1_8o6.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_55_1_YEA.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_56_1_TyQ.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_57_1_5GF.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_58_1_3sI.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_59_1_H1v.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_5_1_uLn.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_60_1_SIM.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_61_1_GOX.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_62_1_l1w.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_63_1_N7r.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_64_1_2TT.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_65_1_R8P.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_66_1_DIf.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_67_1_liJ.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_68_1_RWn.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_69_1_hQ9.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_6_1_oj0.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_70_1_y6x.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_71_1_idP.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_72_1_F9l.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_73_1_UKs.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_74_1_rA8.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_75_1_NiV.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_76_1_nld.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_77_1_Xc8.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_78_1_UCK.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_79_1_l1U.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_7_1_ngT.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_80_1_il0.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_81_1_bvf.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_82_1_Tqr.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_83_1_TvV.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_84_1_dp4.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_85_1_a5U.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_86_1_sxO.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_87_1_o7J.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_88_1_yGO.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_89_1_zyK.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_8_1_XRz.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_90_1_iwz.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_91_1_UQB.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_92_1_eJ5.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_93_1_Mwj.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_94_1_HRj.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_95_1_kTB.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_96_1_LQE.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_97_1_Pnr.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_98_1_QlI.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_99_1_3gt.root',
       '/store/user/zuranski/MH_400_MFF_50_CTau80_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_50_CTau80_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_9_1_d9g.root' ] );


secFiles.extend( [
               ] )


