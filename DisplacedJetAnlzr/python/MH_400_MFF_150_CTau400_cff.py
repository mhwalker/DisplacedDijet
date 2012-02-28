import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_100_1_ny9.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_101_1_pfY.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_102_1_nEv.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_10_1_BNT.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_11_1_YvB.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_12_1_DgN.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_13_1_es6.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_14_1_Pq5.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_15_1_uJV.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_16_1_NTK.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_17_1_yAK.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_18_1_L4z.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_19_1_4NC.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_1_1_Iua.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_20_1_Nrk.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_21_1_NqK.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_22_1_H38.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_23_1_ZJv.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_24_1_Kxf.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_25_1_0Co.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_26_1_f5s.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_27_1_QMZ.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_28_1_Y6F.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_29_1_IVd.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_2_1_qHe.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_30_1_EpG.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_31_1_wgp.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_32_1_iCa.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_33_1_b57.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_34_1_fpI.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_35_1_HgB.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_36_1_Kab.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_37_1_jd9.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_38_1_iPg.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_39_1_5pu.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_3_1_Xzt.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_40_1_pfs.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_41_1_yPE.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_42_1_Yry.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_43_1_Q46.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_44_1_V8g.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_45_1_hva.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_46_1_Nkl.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_47_1_2Gv.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_48_1_54E.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_49_1_Q3D.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_4_1_bhP.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_50_1_hUq.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_51_1_W7R.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_52_1_1Ly.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_53_1_26B.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_54_1_4Mp.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_55_1_4tP.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_56_1_D4O.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_57_1_aXT.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_58_1_D99.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_59_1_5WT.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_5_1_lk9.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_60_1_XzZ.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_61_1_XgX.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_62_1_Qo1.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_63_1_WLR.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_64_1_m5s.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_65_1_8Mm.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_66_1_jFT.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_67_1_Dxj.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_68_1_IlU.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_69_1_KQM.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_6_1_RpQ.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_70_1_Qz9.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_71_1_XFq.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_72_1_aV7.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_73_1_n2r.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_74_1_7Pk.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_75_1_ORV.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_76_1_AA8.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_77_1_I0k.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_78_1_Io3.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_79_1_Yj0.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_7_1_oJk.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_80_1_VIr.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_81_1_h6g.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_82_1_9Dw.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_83_1_pLq.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_84_1_tWH.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_85_1_M7R.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_86_1_ihM.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_87_1_YWl.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_88_1_xuN.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_89_1_Sqs.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_8_1_p1B.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_90_1_YNo.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_91_1_ESV.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_92_1_q0A.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_93_1_jbZ.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_94_1_J6X.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_95_1_RG0.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_96_1_GUk.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_97_1_ggm.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_98_1_aOi.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_99_1_CXe.root',
       '/store/user/zuranski/MH_400_MFF_150_CTau400_7TeVGEN_SIM_RAWDEBUG/MH_400_MFF_150_CTau400_7TeV_GEN_SIM_RECODEBUG/793492cd38bb7dd018eb13f14958f6ca/reco_9_1_MBT.root' ] );


secFiles.extend( [
               ] )

