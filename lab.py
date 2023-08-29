
import const

tt = {'folder':'D:/data/charged_defect/'}

mos2_tt = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2_tt/',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tt/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     3.186*8],
    'sc': [6,6,1],
    'fft_g':[181,181,241]
}

mos2_12x12 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2_12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     3.186*8],
    'sc': [12,12,1],
    'fft_g':[180,180,120]
}

bn_3x3 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/dft/3x3/',
    'folder_dft': tt['folder']+'bn/dft/3x3/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [3,3,1],
    'fft_g':[55,55,181]
}

bn_6x6 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/6x6/',
    'folder_dft': tt['folder']+'bn/dft/6x6/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'folder_out': tt['folder']+'bn/out/6x6',
    'rho_tot': tt['folder']+'bn/dft/6x6/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [6,6,1],
   # 'fft_g':[109,109,181]
    'fft_g':[108,108,180]
}

bn_7x7 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/7x7/',
    'folder_dft': tt['folder']+'bn/dft/7x7/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [7,7,1],
    'fft_g':[128,128,180]
}

bn_8x8 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/8x8/',
    'folder_dft': tt['folder']+'bn/dft/8x8/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [8,8,1],
    'fft_g':[144,144,180]
}

bn_9x9 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/9x9/',
    'folder_dft': tt['folder']+'bn/dft/9x9/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [9,9,1],
    'fft_g':[162,162,180]
}

bn_10x10 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/10x10/',
    'folder_dft': tt['folder']+'bn/dft/10x10/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [10,10,1],
    'fft_g':[180,180,180]
}

bn_11x11 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/11x11/',
    'folder_dft': tt['folder']+'bn/dft/11x11/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [11,11,1],
    'fft_g':[200,200,180]
}

bn_12x12 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/12x12/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [12,12,1],
    'fft_g':[216,216,180],
    'folder_G_info': tt['folder']+'bn_2/G_info/',
    'folder_model': tt['folder']+'bn_2/model_interp/'
}
bn_13x13 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/13x13/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [13,13,1],
    'fft_g':[240,240,180]
}

bn_14x14 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/14x14/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [14,14,1],
    'fft_g':[256,256,180]
}

bn_15x15 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn_2/15x15/',
    'folder_dft': tt['folder']+'bn/dft/15x15/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [15,15,1],
    'fft_g':[270,270,180]
}
