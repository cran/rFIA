utils::globalVariables(c('ADJ_FACTOR_MACR', 'ADJ_FACTOR_MICR', 'ADJ_FACTOR_SUBP', 'ANN_BA_GROWTH', 'x',
                         'ANN_DIA_GROWTH', 'ANN_HT_GROWTH', 'ANN_NET_GROWTH', 'AREA AREA_SE', 'AREA_TOTAL',
                         'AREA_TOTAL_SE', 'AREA_USED', 'BAA', 'BAA_GROW', 'BAA_GROW_SE', 'BAA_PERC', 'BAA_PERC_SE',
                         'BAA_SE', 'BAA_TOTAL', 'BA_GROW', 'BA_GROW_SE', 'BA_SE', 'BA_TOTAL', 'BA_TOTAL_full', 'BIO',
                         'BIO_1000HR', 'BIO_1000HR_ACRE', 'BIO_100HR', 'BIO_100HR_ACRE', 'BIO_10HR',
                         'BIO_10HR_ACRE', 'BIO_1HR', 'BIO_1HR_ACRE', 'BIO_1HR_ACRE_SE', 'BIO_1HR_SE',
                         'BIO_ACRE', 'BIO_ACRE_SE', 'BIO_AG_ACRE', 'BIO_AG_TOTAL', 'BIO_BG_ACRE',
                         'BIO_BG_TOTAL', 'BIO_PILE', 'BIO_PILE_ACRE', 'BIO_PILE_ACRE_SE', 'BIO_PILE_SE',
                         'BIO_TOTAL', 'BIO_TOTAL_SE', 'CARB CARBON_AG', 'CARBON_BG CARB_1000HR',
                         'CARB_1000HR_ACRE', 'CARB_100HR', 'CARB_100HR_ACRE', 'CARB_10HR', 'CARB_10HR_ACRE',
                         'CARB_1HR', 'CARB_1HR_ACRE', 'CARB_1HR_ACRE_SE', 'CARB_1HR_SE', 'CARB_ACRE',
                         'CARB_ACRE_SE', 'CARB_AG_ACRE', 'CARB_AG_TOTAL', 'CARB_BG_ACRE','CARB_BG_TOTAL',
                         'CARB_PILE', 'CARB_PILE_ACRE', 'CARB_PILE_ACRE_SE', 'CARB_PILE_SE', 'CARB_TOTAL',
                         'CARB_TOTAL_SE', 'CCLCD', 'CN', 'CND_CN COMPONENT', 'CONDID', 'CONDPROP_UNADJ',
                         'COND_CHANGE_CD', 'COND_STATUS_CD', 'COND_STATUS_CD.CURR', 'COND_STATUS_CD.PREV',
                         'COUNTYCD', 'COVER_AREA', 'COVER_AREA_SE', 'COVER_PCT', 'COVER_PCT_SE',
                         'CWD_CARBON_ADJ', 'CWD_DRYBIO_ADJ', 'CWD_VOLCF_ADJ', 'DESIGNCD', 'DIA DIA_GROW',
                         'DIA_GROW_SE', 'DIA_TOTAL', 'DRYBIO_AG', 'DRYBIO_BG', 'END_INVYR', 'ESTN_METHOD',
                         'ESTN_UNIT_CN', 'EVALID', 'EVAL_CN', 'EVAL_TYP', 'EXPNS', 'EhCond', 'EhPlot', 'EhStrat', 'Eh_a',
                         'Eh_a_SE', 'Eh_b', 'Eh_g', 'Ehv', 'FUEL_TYPE', 'FWD_LG_CARBON_ADJ', 'FWD_LG_DRYBIO_ADJ',
                         'FWD_LG_VOLCF_ADJ', 'FWD_MD_CARBON_ADJ', 'FWD_MD_DRYBIO_ADJ', 'FWD_MD_VOLCF_ADJ',
                         'FWD_SM_CARBON_ADJ', 'FWD_SM_DRYBIO_ADJ', 'FWD_SM_VOLCF_ADJ', 'GENUS HT_GROW',
                         'HT_GROW_SE', 'HT_TOTAL', 'H_a', 'H_a_SE', 'H_b', 'H_g', 'INVASIVE_SAMPLING_STATUS_CD',
                         'LATE_AREA', 'LATE_AREA_SE', 'LATE_PERC', 'LATE_PERC_SE', 'MACRO_BREAKPOINT_DIA',
                         'MATURE_AREA', 'MATURE_AREA_SE', 'MATURE_PERC', 'MATURE_PERC_SE', 'MORT_PERC',
                         'MORT_PERC_SE', 'MORT_TOTAL', 'MORT_TOTAL_SE', 'MORT_TPA', 'MORT_TPA_SE',
                         'MOSAIC_AREA', 'MOSAIC_AREA_SE', 'MOSAIC_PERC', 'MOSAIC_PERC_SE', 'NETVOL_ACRE',
                         'NETVOL_AC_TOT', 'NETVOL_GROW', 'NETVOL_GROW_AC', 'NETVOL_GROW_AC_SE',
                         'NETVOL_GROW_SE', 'NETVOL_TOTAL', 'P1PNTCNT_EU', 'P1POINTCNT', 'P2POINTCNT',
                         'PERC_AREA', 'PERC_AREA_SE', 'PERC_SE', 'PILE_CARBON_ADJ', 'PILE_DRYBIO_ADJ',
                         'PILE_VOLCF_ADJ', 'PLOT PLT_CN', 'POLE_AREA', 'POLE_AREA_SE', 'POLE_PERC',
                         'POLE_PERC_SE', 'POP_EVAL', 'PROP_BASIS', 'RECR_PERC', 'RECR_PERC_SE', 'RECR_TOTAL',
                         'RECR_TOTAL_SE', 'RECR_TPA', 'RECR_TPA_SE', 'REMPER', 'REMV_PERC', 'REMV_PERC_SE',
                         'REMV_TOTAL', 'REMV_TOTAL_SE', 'REMV_TPA', 'REMV_TPA_SE', 'REPORT_YEAR_NM',
                         'SAWVOL_ACRE', 'SAWVOL_TOTAL', 'SPCD', 'SPECIES', 'STAGE', 'STATECD', 'STRATUM_CN', 'SUBP',
                         'SUBPTYP', 'SUBPTYP_GRM', 'SYMBOL', 'S_a', 'S_a_SE', 'S_b', 'S_g', 'TOTAL_AREA',
                         'TOTAL_AREA_SE', 'TOTAL_SE', 'TOTAL_TPA', 'TPA', 'TPAGROW_UNADJ', 'TPAMORT_UNADJ',
                         'TPAREMV_UNADJ', 'TPA_PERC', 'TPA_PERC_SE', 'TPA_SE', 'TPA_UNADJ', 'TREE', 'TREE_SE',
                         'TREE_TOTAL', 'TREE_TOTAL_SE', 'TREE_TOTAL_full', 'TRE_CN', 'VEG_SPCD', 'VOL', 'VOLCFNET',
                         'VOLCSNET', 'VOL_1000HR', 'VOL_1000HR_ACRE', 'VOL_100HR', 'VOL_100HR_ACRE', 'VOL_10HR',
                         'VOL_10HR_ACRE', 'VOL_1HR', 'VOL_1HR_ACRE', 'VOL_1HR_ACRE_SE', 'VOL_1HR_SE',
                         'VOL_ACRE', 'VOL_ACRE_SE', 'VOL_PILE', 'VOL_PILE_ACRE', 'VOL_PILE_ACRE_SE',
                         'VOL_PILE_SE', 'VOL_TOTAL', 'VOL_TOTAL_SE', 'YEAR', 'a', 'aAdj', 'aDI', 'aEst', 'aStrat', 'aVar',
                         'a_b', 'a_i', 'a_t', 'areaVar', 'av', 'bCV', 'bEst', 'bOut', 'bPlot', 'bStrat', 'bTEst', 'bTPlot',
                         'bTStrat', 'bTVar', 'bTv', 'bV', 'bVar', 'baVar', 'baaEst', 'baaPlot', 'baaStrat', 'baaVar',
                         'baagVar', 'baav', 'bagCV', 'bagEst', 'bagPlot', 'bagStrat', 'bagVar', 'bagaVar', 'bagv', 'bbgCV',
                         'bbgEst', 'bbgPlot', 'bbgStrat', 'bbgVar', 'bbgaVar', 'bbgv', 'bcCV', 'bcEst', 'bcPlot', 'bcStrat',
                         'bcV', 'bcVar', 'bgVar', 'blgCV', 'blgEst', 'blgPlot', 'blgStrat', 'blgV', 'blgVar', 'bmdCV',
                         'bmdEst', 'bmdPlot', 'bmdStrat', 'bmdV', 'bmdVar', 'bpCV', 'bpEst', 'bpPlot', 'bpStrat', 'bpV',
                         'bpVar', 'bsmCV', 'bsmEst', 'bsmPlot', 'bsmStrat', 'bsmV', 'bsmVar', 'btCV', 'btEst', 'btPlot',
                         'btStrat', 'btVar', 'btaVar', 'btv', 'bv', 'cCV', 'cEst', 'cPlot', 'cStrat', 'cV', 'cVar', 'cagCV',
                         'cagEst', 'cagPlot', 'cagStrat', 'cagVar', 'cagaVar', 'cagv', 'cbEst', 'cbgCV', 'cbgEst',
                         'cbgPlot', 'cbgStrat', 'cbgVar', 'cbgaVar', 'cbgv', 'ccCV', 'ccEst', 'ccPlot', 'ccStrat', 'ccV',
                         'ccVar', 'clgCV', 'clgEst', 'clgPlot', 'clgStrat', 'clgV', 'clgVar', 'cmdCV', 'cmdEst', 'cmdPlot',
                         'cmdStrat', 'cmdV', 'cmdVar', 'condArea', 'coverVar', 'cpCV', 'cpEst', 'cpPlot', 'cpStrat', 'cpV',
                         'cpVar', 'csmCV', 'csmEst', 'csmPlot', 'csmStrat', 'csmV', 'csmVar', 'ctCV', 'ctEst', 'ctPlot',
                         'ctStrat', 'ctVar', 'ctaVar', 'ctv', 'cv', 'cvB', 'cvBAA', 'cvBT', 'cvD', 'cvEst', 'cvEst_b', 'cvEst_bT',
                         'cvEst_baa', 'cvEst_bag', 'cvEst_bbg', 'cvEst_bt', 'cvEst_cag', 'cvEst_cbg', 'cvEst_ct',
                         'cvEst_d', 'cvEst_eh', 'cvEst_g', 'cvEst_ga', 'cvEst_h', 'cvEst_hT', 'cvEst_ht', 'cvEst_l',
                         'cvEst_m', 'cvEst_mT', 'cvEst_ma', 'cvEst_mo', 'cvEst_nv', 'cvEst_p', 'cvEst_r', 'cvEst_rT',
                         'cvEst_s', 'cvEst_sv', 'cvEst_t', 'cvEst_tT', 'cvG', 'cvGA', 'cvH', 'cvHT', 'cvL', 'cvM', 'cvMT', 'cvMa',
                         'cvMo', 'cvP', 'cvR', 'cvRT', 'cvS', 'cvStrat', 'cvStrat_Eh', 'cvStrat_b', 'cvStrat_baa',
                         'cvStrat_bag', 'cvStrat_bbg', 'cvStrat_bt', 'cvStrat_cag', 'cvStrat_cbg', 'cvStrat_ct',
                         'cvStrat_d', 'cvStrat_g', 'cvStrat_ga', 'cvStrat_h', 'cvStrat_hT', 'cvStrat_ht',
                         'cvStrat_l', 'cvStrat_m', 'cvStrat_mT', 'cvStrat_ma', 'cvStrat_mo', 'cvStrat_nv',
                         'cvStrat_p', 'cvStrat_r', 'cvStrat_rT', 'cvStrat_s', 'cvStrat_sv', 'cvStrat_t', 'cvT',
                         'cvTT', 'cveH', 'dEst', 'dPlot', 'dStrat', 'dVar', 'db', 'dgVar', 'dv', 'eh', 'ehVar', 'expInv', 'fVar', 'fa',
                         'faFull', 'faFullmo', 'fullDI', 'fullEst', 'fullStrat', 'fullVar', 'fullv', 'gEst', 'gPlot',
                         'gStrat', 'gVar', 'gaEst', 'gaPlot', 'gaStrat', 'gaVar', 'gagVar', 'gav', 'ggVar', 'gv', 'h', 'hCond',
                         'hEst', 'hPlot', 'hStrat', 'hVar', 'haVar', 'htEst', 'htPlot', 'htStrat', 'htVar', 'htgVar', 'htv', 'hv',
                         'iEst', 'iPlot', 'iPlots', 'iStrat', 'iVar', 'iv', 'l', 'lEst', 'lStrat', 'lVar', 'lpVar', 'lv', 'mEst',
                         'mPlot', 'mStrat', 'mVar', 'ma', 'maEst', 'maStrat', 'maVar', 'mapVar', 'mav', 'mo', 'moEst', 'moStrat',
                         'moVar', 'mopVar', 'mov', 'mtVar', 'mv', 'nPlots', 'nPlots_AREA', 'nPlots_INV', 'nPlots_TREE',
                         'nPlots_VOL', 'nStands', 'nh', 'nvCV', 'nvEst', 'nvPlot', 'nvStrat', 'nvVar', 'nvaVar', 'nvv', 'p', 'p1',
                         'p1EU', 'p1EU_a', 'p1_b', 'p1_i', 'p1_t', 'p2', 'p2_b', 'p2_i', 'p2_t', 'pDI', 'pEst', 'pStrat', 'pVar',
                         'plotIn', 'plotIn_AREA', 'plotIn_TREE', 'plotIn_a', 'plotIn_t', 'plotSF', 'plotsInv',
                         'polyID', 'ppVar', 'propVar', 'pv', 'rEst', 'rPlot', 'rStrat', 'rVar', 'raVar', 'rtVar', 'rv', 's', 'sCond',
                         'sPlot', 'sStrat', 'sVar', 'stage', 'sv', 'svCV', 'svEst', 'svPlot', 'svStrat', 'svVar', 'svaVar', 'svv',
                         'tAdj', 'tDI', 'tEst', 'tPlot', 'tStrat', 'tTEst', 'tTPlot', 'tTStrat', 'tTVar', 'tTv', 'tVar', 'totals',
                         'tpVar', 'tpaVar', 'treeVar', 'tv', 'vCV', 'vEst', 'vPlot', 'vStrat', 'vV', 'vVar', 'vcCV', 'vcEst',
                         'vcPlot', 'vcStrat', 'vcV', 'vcVar', 'vlgCV', 'vlgEst', 'vlgPlot', 'vlgStrat', 'vlgV', 'vlgVar',
                         'vmdCV', 'vmdEst', 'vmdPlot', 'vmdStrat', 'vmdV', 'vmdVar', 'vpCV', 'vpEst', 'vpPlot', 'vpStrat',
                         'vpV', 'vpVar', 'vsmCV', 'vsmEst', 'vsmPlot', 'vsmStrat', 'vsmV', 'vsmVar', 'w',
                         'AREA', 'AREA_SE', 'CARB', 'CARBON_AG', 'CARBON_BG', 'CARB_1000HR', 'CND_CN', 'COMPONENT',
                         'DIA', 'DIA_GROW', 'GENUS', 'HT_GROW', 'PLOT', 'PLT_CN', 'maxYear', 'aD', 'LOCATION_NM',
                         'grpVar', 'm', 'xVar', 'yVar', 'LAT', 'LON', 'PREVCOND', 'PREV_PLT_CN', 'STATE',
                         'BIO_DUFF', 'BIO_DUFF_ACRE', 'BIO_DUFF_ACRE_SE', 'BIO_DUFF_SE', 'BIO_LITTER',
                         'BIO_LITTER_ACRE', 'CARB_DUFF', 'CARB_DUFF_ACRE', 'CARB_DUFF_ACRE_SE',
                         'CARB_DUFF_SE', 'CARB_LITTER', 'CARB_LITTER_ACRE', 'DUFF_BIOMASS', 'DUFF_CARBON',
                         'LITTER_BIOMASS', 'LITTER_CARBON', 'SUBP_COMPONENT_AL_FOREST',
                         'SUBP_COMPONENT_AL_TIMBER', 'SUBP_COMPONENT_GS_FOREST',
                         'SUBP_COMPONENT_GS_TIMBER', 'SUBP_SUBPTYP_GRM_AL_FOREST',
                         'SUBP_SUBPTYP_GRM_AL_TIMBER', 'SUBP_SUBPTYP_GRM_GS_FOREST',
                         'SUBP_SUBPTYP_GRM_GS_TIMBER', 'SUBP_TPAGROW_UNADJ_AL_FOREST',
                         'SUBP_TPAGROW_UNADJ_AL_TIMBER', 'SUBP_TPAGROW_UNADJ_GS_FOREST',
                         'SUBP_TPAGROW_UNADJ_GS_TIMBER', 'SUBP_TPAMORT_UNADJ_AL_FOREST',
                         'SUBP_TPAMORT_UNADJ_AL_TIMBER', 'SUBP_TPAMORT_UNADJ_GS_FOREST',
                         'SUBP_TPAMORT_UNADJ_GS_TIMBER', 'SUBP_TPAREMV_UNADJ_AL_FOREST',
                         'SUBP_TPAREMV_UNADJ_AL_TIMBER', 'SUBP_TPAREMV_UNADJ_GS_FOREST',
                         'SUBP_TPAREMV_UNADJ_GS_TIMBER', 'VOL_DUFF', 'VOL_DUFF_ACRE', 'VOL_DUFF_ACRE_SE',
                         'VOL_DUFF_SE', 'bdCV', 'bdEst', 'bdPlot', 'bdStrat', 'bdV', 'bdVar', 'blCV', 'blEst', 'blPlot',
                         'blStrat', 'blV', 'blVar', 'cdCV', 'cdEst', 'cdPlot', 'cdStrat', 'cdV', 'cdVar', 'clCV', 'clEst',
                         'clPlot', 'clStrat', 'clV', 'clVar', 'dbConnect', 'dbListTables', 'nPlots_MORT',
                         'nPlots_RECR', 'nPlots_REMV', 'p1_eu', 'plotIn_g', 'plotIn_h', 'plotIn_m', 'plotIn_r', 'progress', 'totalPlots',
                         'COND_STATUS_CD.chng', 'COND_STATUS_CD.prev', 'GROWTH_ACCT', 'geometry', 'lazy_dt',
                         'BIO_AC_TOT', 'BIO_GROW', 'BIO_GROW_AC', 'BIO_GROW_AC_SE', 'BIO_GROW_SE', 'DIA.mid',
                         'DIA.prev', 'DRYBIO_AG.mid', 'DRYBIO_AG.prev', 'INVYR', 'VOLCFNET.mid',
                         'VOLCFNET.prev', 'bioAPlot', 'bioEst', 'bioPlot', 'bioStrat', 'bioVar', 'bioaEst',
                         'bioaStrat', 'bioaVar', 'bioagVar', 'bioav', 'biogVar', 'biov', 'cvBIO', 'cvBIOA', 'cvEst_bio',
                         'cvEst_bioa', 'cvStrat_bio', 'cvStrat_bioa',
                         'BIO_AG_ACRE_SE', 'BIO_BG_ACRE_SE', 'CARB_AG_ACRE_SE', 'CARB_BG_ACRE_SE',
                         'NETVOL_ACRE_SE', 'SAWVOL_ACRE_SE', 'SUBPTYP_PROP_CHNG',
                         'BA', 'BA1', 'BA2', 'DIA.beg', 'DIA1', 'DIA2', 'DRYBIO_AG.beg', 'DRYBIO_AG1', 'DRYBIO_AG2',
                         'ONEORTWO', 'VOLCFNET.beg', 'VOLCFNET1', 'VOLCFNET2', 'VOLCSNET.beg', 'VOLCSNET.mid',
                         'VOLCSNET.prev', 'ga', 'se', 'ba', 'baa', 'd', 'vol', 'volA', 'place',
                         'ANN_INVENTORY', 'MANUAL', 'RSCD', 'SEEDLING', 'SEEDLING_SE', 'SEEDLING_TOTAL',
                         'SEEDLING_TOTAL_full', 'STATENM', 'nPlots_SEEDLING', 'plotIn_SEEDLING',
                         'BIO_SE', 'CARB_SE','MEASYEAR', 'PERC', 'PLOT_BASIS', 'TREECOUNT_CALC', 'UNITCD', 'VOL_SE',
                         'aD_c', 'bbaVar', 'caaVar', 'cbaVar', 'cvEst_bc', 'cvEst_bd', 'cvEst_bl', 'cvEst_blg',
                         'cvEst_bmd', 'cvEst_bp', 'cvEst_bsm', 'cvEst_c', 'cvEst_cc', 'cvEst_cd', 'cvEst_cl',
                         'cvEst_clg', 'cvEst_cmd', 'cvEst_cp', 'cvEst_csm', 'cvEst_v', 'cvEst_vc', 'cvEst_vlg',
                         'cvEst_vmd', 'cvEst_vp', 'cvEst_vsm', 'cvStrat_eh', 'ehEst', 'ehStrat', 'ehaVar', 'ehv', 'grp',
                         'landD', 'nPlots_DWM', 'ndif', 'p2eu', 'pltID', 'prev', 'r_t', 'sEst', 'saVar', 'sp', 'state', 'typeD',
                         'wgt', 'yrPrev', 'INV_AREA_TOTAL', 'SUBPCOND_PROP', 'cover', 'cvEst_i', 'cvStrat_i', 'iTEst', 'iTStrat',
                         'plotIn_INV', 'minyr', "TPARECR_UNADJ", 'aDI_ga', 'cvEst_hp', 'cvEst_mp', 'cvEst_rp', 'cvStrat_hp', 'cvStrat_mp',
                         'cvStrat_rp', 'fa_ga', 'hPlot_ga', 'hpVar', 'mPlot_ga', 'mpVar', 'plotIn_ga', 'rPlot_ga',
                         'rpVar', 'tDI_ga', 'tPlot_ga', 'BAA_TOTAL_SE', 'BA_GROW_AC', 'BA_GROW_AC_SE', 'BA_TOTAL_SE', 'BIOA_TOTAL',
                         'BIOA_TOTAL_SE', 'DIA_TOTAL_SE', 'NETVOLA_TOTAL', 'NETVOLA_TOTAL_SE',
                         'NETVOL_TOTAL_SE', 'PLOT_STATUS_CD', 'bioAEst', 'bioAStrat', 'bioAVar', 'bioAgVar',
                         'bioAv', 'cvEst_bioA', 'cvStrat_bioA', 'tDI_ga_r', 'tDI_r', 'polys', 'diffYear', 'state_recr',
                         'P2POINTCNT_INVYR', 'sumwgt', 'seVar',
                         'BAA_RATE', 'BAA_RATE_SE', 'BAA_RATE_ga', 'BAA_ga', 'DIA_BEGIN', 'DIA_END', 'HARV_BAA',
                         'HARV_RATE', 'HARV_RATE_SE', 'HARV_RATE_ga', 'LAMBDA', 'LAMBDA_SE', 'LAMBDA_ga', 'M',
                         'MORT_BAA', 'MORT_RATE', 'MORT_RATE_SE', 'MORT_RATE_ga', 'MORT_TPA_ga', 'M_ga', 'P1',
                         'PREV_BAA', 'PREV_BAA_ga', 'PREV_TPA', 'PREV_TPA_ga', 'Q1', 'RECR_RATE', 'RECR_RATE_SE',
                         'RECR_RATE_ga', 'RECR_TPA_ga', 'REMV_TPA_ga', 'S1', 'STATUSCD', 'STATUSCD.prev',
                         'SUST_INDEX', 'SUST_INDEX_SE', 'TPA_UNADJ.prev', 'cvEst_q', 'cvStrat_q', 'lPlot', 'nLive',
                         'nStems', 'qEst', 'qStrat', 'qVar', 'qv', 'x_ga', 'y', 'y_ga',
                         'BAA1', 'BAA2', 'BAA_RATE_INT', 'BAA_STATUS', 'CHNG_BAA','CHNG_BAA_SE', 'CHNG_TPA',
                         'CHNG_TPA_SE', 'DIA_MIDPT', 'M_TPA', 'PLOT_STATUS_CD1', 'PLOT_STATUS_CD2', 'PREVDIA',
                         'PREV_BAA_SE', 'PREV_HARV_BAA', 'PREV_MORT_BAA', 'PREV_TPA_SE', 'PREV_TRE_CN',
                         'RECR_BAA', 'REMV_BAA', 'R_TPA', 'SI_STATUS', 'TPA_RATE', 'TPA_RATE_INT', 'TPA_RATE_SE',
                         'TPA_STATUS', 'TPA_UNADJ1', 'TPA_UNADJ2', 'T_TPA', 'cb', 'cbPlot', 'cbStrat', 'cbVar', 'cbv', 'ct',
                         'cvEst_cb', 'cvStrat_cb', 'pb', 'pbEst', 'pbPlot', 'pbStrat', 'pbVar', 'pbv', 'pt', 'ptEst', 'ptPlot',
                         'ptStrat', 'ptVar', 'ptv', 'tDI1', 'tDI2', 'AGENTCD', 'ANIMAL', 'ANIMAL_RATE', 'ANIMAL_RATE_SE',
                         'BUG', 'BUG_RATE', 'DISEASE',
                         'DISEASE_RATE', 'DISEASE_RATE_SE', 'FIRE', 'FIRE_RATE', 'FIRE_RATE_SE', 'INSECT_RATE',
                         'INSECT_RATE_SE', 'MORTYR', 'SILV', 'SILV_RATE', 'SILV_RATE_SE', 'UNKNOWN',
                         'UNKNOWN_RATE', 'UNKNOWN_RATE_SE', 'VEG', 'VEG_RATE', 'VEG_RATE_SE', 'WEATHER',
                         'WEATHER_RATE', 'WEATHER_RATE_SE', 'animal', 'animalEst', 'animalPlot', 'animalStrat',
                         'animalVar', 'animalv', 'bug', 'bugEst', 'bugPlot', 'bugStrat', 'bugVar', 'bugv',
                         'cvEst_animal', 'cvEst_bug', 'cvEst_disease', 'cvEst_fire', 'cvEst_silv', 'cvEst_un',
                         'cvEst_veg', 'cvEst_weather', 'cvStrat_animal', 'cvStrat_bug', 'cvStrat_disease',
                         'cvStrat_fire', 'cvStrat_silv', 'cvStrat_un', 'cvStrat_veg', 'cvStrat_weather',
                         'disease', 'diseaseEst', 'diseasePlot', 'diseaseStrat', 'diseaseVar', 'diseasev', 'fire',
                         'fireEst', 'firePlot', 'fireStrat', 'fireVar', 'firev', 'silv', 'silvEst', 'silvPlot',
                         'silvStrat', 'silvVar', 'silvv', 'unEst', 'unPlot', 'unStrat', 'unVar', 'unknown', 'unv', 'veg',
                         'vegEst', 'vegPlot', 'vegStrat', 'vegVar', 'vegv', 'weather', 'weatherEst', 'weatherPlot',
                         'weatherStrat', 'weatherVar', 'weatherv',   'ALL_SEV', 'ALL_SEV_SE', 'DROUGHT_SEV', 'DROUGHT_SEV_SE', 'GROW_ALL_SEV',
                         'GROW_ALL_SEV_SE', 'GROW_DROUGHT_SEV', 'GROW_DROUGHT_SEV_SE', 'GROW_WET_SEV',
                         'GROW_WET_SEV_SE', 'WET_SEV', 'WET_SEV_SE', 'all_sev', 'cvEst_a', 'cvEst_gd', 'cvEst_gw',
                         'cvEst_w', 'cvStrat_a', 'cvStrat_gd', 'cvStrat_gw', 'cvStrat_w', 'drought_sev', 'faEst',
                         'faStrat', 'faVar', 'fav', 'gdEst', 'gdStrat', 'gdVar', 'gdv', 'grow_all_sev',
                         'grow_drought_sev', 'grow_wet_sev', 'gwEst', 'gwStrat', 'gwVar', 'gwv', 'wEst', 'wStrat',
                         'wVar', 'wet_sev', 'wv', 'MORT', 'STATUSCD1', 'STATUSCD2', 'cvEst_mort', 'cvStrat_mort', 'mort', 'mortEst',
                         'mortPlot', 'mortStrat', 'mortVar', 'mortv',
                         'GROW_TMAX_ANOM', 'GROW_TMAX_ANOM_SE', 'GROW_TMEAN_ANOM', 'GROW_TMEAN_ANOM_SE',
                         'GROW_VPD_ANOM', 'GROW_VPD_ANOM_SE', 'TMAX_ANOM', 'TMAX_ANOM_SE', 'TMEAN_ANOM',
                         'TMEAN_ANOM_SE', 'VPD_ANOM', 'VPD_ANOM_SE', 'cvEst_gtmax', 'cvEst_gtmean',
                         'cvEst_gvpd', 'cvEst_tmax', 'cvEst_tmean', 'cvEst_vpd', 'cvStrat_gtmax',
                         'cvStrat_gtmean', 'cvStrat_gvpd', 'cvStrat_tmax', 'cvStrat_tmean', 'cvStrat_vpd',
                         'grow_tmax_anom', 'grow_tmean_anom', 'grow_vpd_anom', 'gtmaxEst', 'gtmaxStrat',
                         'gtmaxVar', 'gtmaxv', 'gtmeanEst', 'gtmeanStrat', 'gtmeanVar', 'gtmeanv', 'gvpdEst',
                         'gvpdStrat', 'gvpdVar', 'gvpdv', 'tmaxEst', 'tmaxStrat', 'tmaxVar', 'tmax_anom', 'tmaxv',
                         'tmeanEst', 'tmeanStrat', 'tmeanVar', 'tmean_anom', 'tmeanv', 'vpdEst', 'vpdStrat', 'vpdVar',
                         'vpd_anom', 'vpdv', 'CURR_BAA', 'CURR_TPA', 'SUST_INDEX_INT', 'cvEst_si', 'cvStrat_si', 'siEst', 'siPlot',
                         'siStrat', 'siVar', 'siv', 'BAA_CHNG_PERC', 'BAA_MORT', 'BAA_MORT_SE', 'BAA_RECR', 'BAA_RECR_SE', 'ELEV ELEV_SE',
                         'SI', 'SI_INT', 'SI_SE', 'SURV', 'TPA_CHNG_PERC', 'TPA_MORT', 'TPA_MORT_SE', 'TPA_RECR',
                         'TPA_RECR_SE', 'bgrow', 'bgrowEst', 'bgrowStrat', 'bgrowVar', 'bgrowv', 'bmort', 'bmortEst',
                         'bmortStrat', 'bmortVar', 'bmortv', 'brecr', 'brecrEst', 'brecrStrat', 'brecrVar', 'brecrv',
                         'cvEst_bgrow', 'cvEst_bmort', 'cvEst_brecr', 'cvEst_elev', 'cvEst_tmort',
                         'cvEst_trecr', 'cvStrat_bgrow', 'cvStrat_bmort', 'cvStrat_brecr', 'cvStrat_elev',
                         'cvStrat_tmort', 'cvStrat_trecr', 'elevEst', 'elevStrat', 'elevVar', 'elevv', 'tmort',
                         'tmortEst', 'tmortStrat', 'tmortVar', 'tmortv', 'trecr', 'trecrEst', 'trecrStrat',
                         'ALSTK', 'BAA_RATE_VAR', 'ELEV', 'ELEV_SE', 'Eh_a_species', 'Eh_a_species_SE',
                         'Eh_a_struct', 'Eh_a_struct_SE', 'Eh_b_species', 'Eh_b_struct', 'Eh_g_species',
                         'Eh_g_struct', 'Eh_species', 'Eh_struct', 'H_a_species', 'H_a_species_SE',
                         'H_a_struct', 'H_a_struct_SE', 'H_b_species', 'H_b_struct', 'H_g_species',
                         'H_g_struct', 'H_species', 'H_struct', 'N', 'SI_VAR', 'SPCD1', 'SPCD2', 'STOCKING',
                         'STOCKING_SE', 'S_a_species', 'S_a_species_SE', 'S_a_struct', 'S_a_struct_SE',
                         'S_b_species', 'S_b_struct', 'S_g_species', 'S_g_struct', 'S_species', 'S_struct',
                         'TPA_RATE_VAR', 'baaMean', 'baaSD', 'cvEst_esp', 'cvEst_est', 'cvEst_hsp', 'cvEst_hst',
                         'cvEst_ssp', 'cvEst_sst', 'cvEst_stk', 'cvStrat_esp', 'cvStrat_est', 'cvStrat_hsp',
                         'cvStrat_hst', 'cvStrat_ssp', 'cvStrat_sst', 'cvStrat_stk', 'espEst', 'espStrat',
                         'espVar', 'espv', 'estEst', 'estStrat', 'estVar', 'estv', 'hspEst', 'hspStrat', 'hspVar', 'hspv',
                         'hstEst', 'hstStrat', 'hstVar', 'hstv', 'htClass', 'minLive', 'nTotal', 'qt', 'sd', 'sspEst',
                         'sspStrat', 'sspVar', 'sspv', 'sstEst', 'sstStrat', 'sstVar', 'sstv', 'stkEst', 'stkStrat',
                         'stkVar', 'stkv', 'tpaMean', 'tpaSD', 'trecrVar', 'trecrv', 'ymax', 'ymin', 'EVAL_DESCR',
                         'AG_LIVE', 'AG_OVER_DEAD', 'AG_OVER_LIVE', 'AG_UNDER_LIVE', 'BG_LIVE', 'BG_OVER_DEAD',
                         'BG_OVER_LIVE', 'BG_UNDER_LIVE', 'CARBON_DOWN_DEAD', 'CARBON_LITTER',
                         'CARBON_SOIL_ORG', 'CARBON_STANDING_DEAD', 'CARBON_UNDERSTORY_AG',
                         'CARBON_UNDERSTORY_BG', 'DEAD_WOOD', 'DOWN_DEAD', 'LITTER', 'SOIL_ORG', 'STAND_DEAD',
                         'STAND_DEAD_MOD', 'all_of', 'caVar', 'cvStrat_c', 'dead', 'live',
                         'MEASDAY', 'MEASMON', 'b_rate', 'coef', 'df', 'lm', 'meas', 't_rate', 'treID',
                         'aD_p', 'amax', 'amin', 'maxDate', 'minDate', 'treID1', 'treID2', 'STATEAB', 'pops',
                         'grpC', 'grpP', 'grpT', 'GROWTH_HABIT', 'LAYER', 'P2VEG_SAMPLING_STATUS_CD', 'TOTAL_COVER_AREA',
                         'cTEst', 'cTStrat', 'plotIn_VEG', 'BA_RATE', 'BA_RATE_SE', 'BA_RATE_VAR', 'CHNG_BA',
                         'CHNG_QMD', 'CURR_BA', 'CURR_QMD',
                         'DESIGNCD.prev', 'DESIGNCD1', 'FSI', 'FSI_INT', 'FSI_SE', 'FSI_STATUS', 'FSI_VAR',
                         'PERC_FSI', 'PERC_FSI_INT', 'PERC_FSI_SE', 'PERC_FSI_VAR', 'PLOT_STATUS_CD.prev',
                         'PREV_BA', 'PREV_QMD', 'b1', 'b2', 'baSD', 'brate', 'cvEst_psi', 'cvStrat_psi', 'dt',
                         'psiVar', 'si',
                         'si1', 'si1Est', 'si1Strat', 'si1Var', 'si1v', 't1', 't2', 'trate',
                         '.', 'FSI1', 'FSI2', 'TPA1', 'b', 'grps', 'int', 'pDI1', 'pDI2', 'ranef',
                         'rate', 'slope',
                         'bSD', 'tSD', 'DSTRBCD1', 'DSTRBCD2', 'DSTRBCD3', 'TRTCD1', 'TRTCD2', 'TRTCD3',
                         'na.omit', 'ra', 'sds', 'teEst','teStrat', 'teVar', 'tev', 'texpect', 'tmax1', 'tmax2',
                         'meanBA', 'percentile', 'CURR_RD', 'CURR_RD_VAR', 'PREV_RD', 'PREV_RD_VAR', 'cvEst_ra1',
                         'cvEst_ra2', 'cvStrat_ra1', 'cvStrat_ra2', 'ra1', 'ra1Est', 'ra1Strat', 'ra1Var', 'ra1v', 'ra2', 'ra2Est',
                         'ra2Strat', 'ra2Var', 'ra2v', 'TPA2', 'disturb', 'int_fixed', 'rate_fixed', 'rd', 'tmax',
                         'estimate', 'fe_int', 'fe_rate', 'grp_index', 'lower', 'mean_fe_int',  'mean_fe_rate',
                         'mean_int', 'mean_rate', 'quantile', 'term', 'skew2',
                         'fsiHelper1_lm',  'fullRemp',  'series',
                         'AREA_TOTAL_VAR', 'BAA_PERC_VAR', 'BAA_TOTAL_VAR', 'BAA_VAR', 'BA_GROW_AC_VAR',
                         'BA_GROW_VAR', 'BA_TOTAL_VAR', 'BA_VAR', 'BIOA_TOTAL_VAR', 'BIO_GROW_AC_VAR',
                         'BIO_GROW_VAR', 'BIO_TOTAL_VAR', 'DIA_GROW_VAR', 'DIA_TOTAL_VAR', 'Eh_a_VAR',
                         'H_a_VAR', 'MORT_PERC_VAR', 'MORT_TOTAL_VAR', 'MORT_TPA_VAR', 'NETVOLA_TOTAL_VAR',
                         'NETVOL_GROW_AC_VAR', 'NETVOL_GROW_VAR', 'NETVOL_TOTAL_VAR', 'RECR_PERC_VAR',
                         'RECR_TOTAL_VAR', 'RECR_TPA_VAR', 'REMV_PERC_VAR', 'REMV_TOTAL_VAR', 'REMV_TPA_VAR',
                         'S_a_VAR', 'TPA_PERC_VAR', 'TPA_VAR', 'TREE_TOTAL_VAR', 'TREE_VAR',
                         'mr', 'adj', 'cvStrat_bT', 'cvStrat_tT', 'tpaHelper1_test', 'tpaHelper2_test',
                         'P1POINTCNT_INVYR', 'PLOT_STATUS', 'STRATUM_DESCR', 'adist', 'aes', 'anim_save', 'bT',
                         'buff', 'cvStrat_bc', 'cvStrat_bd', 'cvStrat_bl', 'cvStrat_blg',
                         'cvStrat_bmd', 'cvStrat_bp', 'cvStrat_bsm', 'cvStrat_cc', 'cvStrat_cd', 'cvStrat_cl',
                         'cvStrat_clg', 'cvStrat_cmd', 'cvStrat_cp', 'cvStrat_csm', 'cvStrat_v', 'cvStrat_vc',
                         'cvStrat_vlg', 'cvStrat_vmd', 'cvStrat_vp', 'cvStrat_vsm', 'element_rect',
                         'element_text', 'facet_wrap', 'geom_errorbar', 'geom_line', 'geom_sf', 'ggplot', 'ggsave',
                         'ggtitle', 'installed.packages', 'labs', 'nsum', 'plotIn_', 'proj4string<-',
                           'propSampled', 'scale_colour_viridis_c', 'scale_colour_viridis_d',
                         'scale_fill_viridis_c', 'stratID', 'stratWgt', 'stratWgt_INVYR', 'tD', 'tT', 'theme',
                         'theme_bw', 'theme_minimal', 'transition_manual', 'unit', 'waiver', 'xlab', 'ylab',
                         'lower_fe_int', 'lower_fe_rate', 'lower_int', 'lower_rate', 'upper_fe_int', 'fixed_rate',
                         'upper_fe_rate', 'upper_int', 'upper_rate', 'keep', 'nplts', 'p2eu_INVYR', 'fixed_alpha',
                         'PERC_AREA_VAR', 'aCV', 'acv', 'at', 'atEst', 'atStrat', 'atVar', 'atv', 'good',
                         'REMPER_VAR', 'cvEst_remp', 'cvStrat_remp', 'rempEst', 'rempStrat', 'rempVar', 'rempv',
                         'AREA_CHNG', 'AREA_CHNG_SE', 'AREA_CHNG_VAR', 'CONDID1', 'CONDID2', 'CONDPROP_CHNG',
                         'PERC_CHNG', 'PERC_CHNG_SE', 'PERC_CHNG_VAR', 'PREV_AREA', 'PREV_AREA_SE',
                         'PREV_AREA_VAR', 'PREV_CONDPROP', 'aD_c1', 'aD_c2', 'aD_p1', 'aD_p2', 'ac', 'cov_cv',
                         'dropThese', 'landD1', 'landD2', 'sp1', 'sp2', 'tD1', 'tD2',
                         'bio', 'MORT_TREE_TOTAL', 'MORT_TREE_TOTAL_SE', 'MORT_TREE_TOTAL_VAR',
                         'RECR_TREE_TOTAL', 'RECR_TREE_TOTAL_SE', 'RECR_TREE_TOTAL_VAR',
                         'REMV_TREE_TOTAL', 'REMV_TREE_TOTAL_SE', 'REMV_TREE_TOTAL_VAR', 'aZero', 'macro',
                         'AREA_DOMAIN1', 'AREA_DOMAIN2', 'TREE_DOMAIN1', 'TREE_DOMAIN2',
                         'COMMON_NAME', 'DRYBIO', 'DRYBIO_', 'DRYBIO_BOLE', 'DRYBIO_FOLIAGE', 'DRYBIO_SAPLING',
                         'DRYBIO_STUMP', 'DRYBIO_TOP', 'DRYBIO_WDLD_SPP', 'JENKINS_FOLIAGE_RATIO_B1',
                         'JENKINS_FOLIAGE_RATIO_B2', 'JENKINS_STEM_BARK_RATIO_B1',
                         'JENKINS_STEM_BARK_RATIO_B2', 'JENKINS_STEM_WOOD_RATIO_B1',
                         'JENKINS_STEM_WOOD_RATIO_B2', 'JENKINS_TOTAL_B1', 'JENKINS_TOTAL_B2',
                         'barkRatio', 'jTotal', 'leafRatio', 'stemRatio',
                         'BOLE_CF_ACRE', 'BOLE_CF_ACRE_SE', 'BOLE_CF_ACRE_VAR', 'BOLE_CF_TOTAL',
                         'BOLE_CF_TOTAL_SE', 'BOLE_CF_TOTAL_VAR', 'SAW_CF_ACRE', 'SAW_CF_ACRE_SE',
                         'SAW_CF_ACRE_VAR', 'SAW_CF_TOTAL', 'SAW_CF_TOTAL_SE', 'SAW_CF_TOTAL_VAR',
                         'SAW_MBF_ACRE', 'SAW_MBF_ACRE_SE', 'SAW_MBF_ACRE_VAR', 'SAW_MBF_TOTAL',
                         'SAW_MBF_TOTAL_SE', 'SAW_MBF_TOTAL_VAR', 'bcf', 'bcfEst', 'bcfPlot', 'bcfStrat', 'bcfVar',
                         'bcfVar_ratio', 'bcfv', 'cvEst_bcf', 'cvEst_sbf', 'cvEst_scf', 'cvStrat_bcf',
                         'cvStrat_sbf', 'cvStrat_scf', 'sbf', 'sbfEst', 'sbfPlot', 'sbfStrat', 'sbfVar',
                         'sbfVar_ratio', 'sbfv', 'scf', 'scfEst', 'scfPlot', 'scfStrat', 'scfVar', 'scfVar_ratio',
                         'scfv', 'SAWVOL_GROW', 'SAWVOL_GROW_AC', 'SAWVOL_GROW_AC_SE', 'SAWVOL_GROW_AC_VAR',
                         'SAWVOL_GROW_SE', 'SAWVOL_GROW_VAR', 'SAWVOL_TOTAL_SE', 'SAWVOL_TOTAL_VAR',
                         'VOLBFNET', 'VOLBFNET.beg', 'VOLBFNET.mid', 'VOLBFNET.prev', 'VOLBFNET1', 'VOLBFNET2',
                         'cvEst_sa', 'cvStrat_sa', 'sagVar', 'sgVar', 'svol', 'jBoleBio', 'nStrata', 'nStrata_INVYR',
                         'SUBP_COMPONENT_SL_FOREST', 'SUBP_COMPONENT_SL_TIMBER',
                         'SUBP_SUBPTYP_GRM_SL_FOREST', 'SUBP_SUBPTYP_GRM_SL_TIMBER',
                         'SUBP_TPAGROW_UNADJ_SL_FOREST', 'SUBP_TPAGROW_UNADJ_SL_TIMBER',
                         'SUBP_TPAMORT_UNADJ_SL_FOREST', 'SUBP_TPAMORT_UNADJ_SL_TIMBER',
                         'SUBP_TPAREMV_UNADJ_SL_FOREST', 'SUBP_TPAREMV_UNADJ_SL_TIMBER',
                         'STRATUM_WGT', 'A', 'cut.these',
                         'AREA_BASIS', 'BIO_ACRE_VAR', 'BIO_cv', 'BIO_mean', 'BIO_var', 'CARB_ACRE_VAR',
                         'CARB_TOTAL_VAR', 'CARB_cv', 'CARB_mean', 'CARB_var', 'COVER_AREA_TOTAL',
                         'COVER_AREA_TOTAL_SE', 'COVER_AREA_TOTAL_VAR', 'COVER_PCT_VAR',
                         'ESTN_UNIT_DESCR', 'Eh', 'Eh_cv', 'Eh_mean', 'Eh_var', 'H', 'H_cv', 'H_mean', 'H_var',
                         'INV_AREA_TOTAL_SE', 'INV_AREA_TOTAL_VAR', 'P2PNTCNT_EU', 'P2PNTCNT_EU_INVYR',
                         'PROP_FOREST', 'SCIENTIFIC_NAME', 'STAGE_AREA_TOTAL', 'STAGE_AREA_TOTAL_SE',
                         'STAGE_AREA_TOTAL_VAR', 'STATUS1', 'STATUS2', 'S_cv', 'S_mean', 'S_var', 'TREE_BASIS',
                         'VOL_ACRE_VAR', 'VOL_TOTAL_VAR', 'VOL_cv', 'VOL_mean', 'VOL_var', 'aChng', 'aD.prev', 'aD1',
                         'aD2', 'aGrpBy', 'ac_cv', 'ac_mean', 'ac_var', 'bPlot_cv', 'bPlot_cv_t', 'bPlot_mean',
                         'bPlot_var', 'bcf_cv', 'bcf_mean', 'bcf_var', 'bioPlot_cv', 'bioPlot_cv_t',
                         'bioPlot_mean', 'bioPlot_var', 'cPlot_cv', 'cPlot_mean', 'cPlot_var', 'cover_cv',
                         'cover_mean', 'cover_var', 'dPlot_cv', 'dPlot_cv_t', 'dPlot_mean', 'dPlot_var', 'fa_cv',
                         'fa_mean', 'fa_var', 'fad', 'fad_mean', 'fad_var', 'full.area', 'gPlot_cv', 'gPlot_cv_t',
                         'gPlot_mean', 'gPlot_var', 'grpByOrig', 'hPlot_cv', 'hPlot_cv_t', 'hPlot_mean',
                         'hPlot_var', 'landD.prev', 'mPlot_cv', 'mPlot_cv_t', 'mPlot_mean', 'mPlot_var', 'method',
                         'nPlots.x', 'nPlots.y', 'nPlots_VEG', 'prev_mean', 'prev_var', 'rPlot_cv', 'rPlot_cv_t',
                         'rPlot_mean', 'rPlot_var', 'sPlot_cv', 'sPlot_cv_t', 'sPlot_mean', 'sPlot_var', 'sbf_cv',
                         'sbf_mean', 'sbf_var', 'scf_cv', 'scf_mean', 'scf_var', 'sp.prev', 'tChng', 'tD.prev',
                         'tPlot_cv', 'tPlot_mean', 'tPlot_var', 'typeD.prev',
                         'CHNG_PERC', 'CHNG_PERC_SE', 'CHNG_PERC_VAR', 'CHNG_TOTAL', 'CHNG_TOTAL_SE',
                         'CHNG_TOTAL_VAR', 'CURR_TOTAL', 'CURR_TOTAL_SE', 'CURR_TOTAL_VAR', 'GROW_PERC',
                         'GROW_PERC_VAR', 'GROW_TOTAL', 'GROW_TOTAL_VAR', 'GROW_TPA', 'GROW_TPA_VAR',
                         'PREV_TOTAL', 'PREV_TOTAL_SE', 'PREV_TOTAL_VAR', 'cPlot_cv_t', 'pPlot', 'pPlot_mean',
                         'pPlot_var', 'state.prev'

))
