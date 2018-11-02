from PYLIQ.PYLIQ_functions import Youd_2002, single_plot, SUMMARY_1, SUMMARY_2, SUMMARY_3
import pandas as pd
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import warnings
import os
name = "PYLIQ"

def analyze(input_file, file_dir):
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    output_dir = os.path.join(file_dir, 'PYLIQ Output/')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    seed_87_fc = [0, 10, 25, 50, 75]
    seed_87_delta_N160 = [0, 1, 2, 4, 5]
    plt.rc('font', family='arial', size=13)
    plt.rc('legend', fontsize=10)
    plt.rc('xtick', labelsize=11.5)
    plt.rc('xtick.major', size=5)
    plt.rc('ytick', labelsize=11.5)
    plt.rc('ytick.major', size=5)

    depth_all, subdivided_wt, ws_all, LPI_all, LSN_all, LDI_all, LD_Youd_all, Sv1D_all, H1_all, H2_all, \
    fines_susceptibility_all, wc_LL_PI_all = ([] for i in range(12))

    # General input parameters
    Water_unit_weight = 9.81
    Pa = 101.325
    soil_types = ['Sand', 'Gravel', 'Clay', 'Silt', 'Peat', 'Organic', 'Other', np.NaN]

    file = pd.ExcelFile(os.path.join(file_dir, input_file))
    input_df = pd.read_excel(file, 'SUMMARY')
    input_df.set_index('Sheet', inplace=True)
    sheet_list = input_df.index.values.tolist()
    options_df = pd.read_excel(file, 'OPTIONS', header=None, index_col=0)
    anal_depth = options_df.loc['Liquefaction max. depth (m)'].values[0]
    LD_min_thick = options_df.loc['Lateral displacement minimum thickness (m)'].values[0]
    calculate_Youd = options_df.loc['Calculate LD from Youd (2002)'].values[0]
    calculate_fines = options_df.loc['Calculate fines liq. Susceptibility from Bray & Sancio (2006)'].values[0]
    calculate_Sr = options_df.loc['Calculate post-liquefaction residual strength, Sr'].values[0]
    OS_wf = options_df.loc['Weighting factor for Olson & Stark (2002) model'].values[0]
    IB_wf = options_df.loc['Weighting factor for Idriss & Boulanger (2008) model'].values[0]
    W_wf = options_df.loc['Weighting factor for Weber (2015) model'].values[0]
    export_pdf = options_df.loc['Export figures as .PDF'].values[0]
    export_png = options_df.loc['Export figures as .PNG'].values[0]
    export_xlsx = options_df.loc['Export results as .XLSX'].values[0]
    writer = pd.ExcelWriter(output_dir + 'Output.xlsx', engine='openpyxl')

    for ws in sheet_list:
        if input_df.loc[ws, 'run'] != 'YES': continue
        ID = ws
        water_table, Mw, PGA, ER, borehole_diameter, rod_extension, site_condition, slope, H_free_face, L_free_face = \
            input_df.loc[
                ws, ['Water table (m)', 'Mw', 'PGA (g)', 'ER (%)', 'Borehole diameter (mm)', 'Rod extension (m)',
                     'Site condition', 'Ground slope (%)', 'H free face', 'L free face']]
        gen_parameters = 'I&B (2014): Mw = ' + str(Mw) + '; PGA (g) = ' + str(PGA)
        df = pd.read_excel(file, ws)
        rows = df.shape[0]
        water_table_idx = np.searchsorted(df['Depth (m)'], water_table).item()
        soil_type = df['Soil type (USCS)'].str.split('-').str[0].replace(
            {'SM': "Sand", 'SC': "Sand", 'SP': "Sand", 'SW': "Sand", 'GM': "Gravel", 'GC': "Gravel", 'GP': "Gravel",
             'GW': "Gravel", 'CH': "Clay", 'CL': "Clay", 'MH': "Silt", 'ML': "Silt", 'PT': "Peat", 'OH': "Organic",
             'OL': "Organic"})

        soil_type[~soil_type.isin(soil_types)] = 'Other'
        if water_table_idx != 0 and soil_type[water_table_idx] == 'Sand' and df['Depth (m)'][
            water_table_idx] != water_table:
            subdivided_wt.append('YES')
            repeat_array = ([1] * rows)
            repeat_array[water_table_idx] = 2
            soil_type = np.repeat(soil_type, repeat_array).reset_index(drop=True)
            df = pd.DataFrame(np.repeat(df.values, repeat_array, axis=0), columns=df.columns)
            df['Depth (m)'][water_table_idx] = water_table
            rows = df.shape[0]
        else:
            subdivided_wt.append('NO')

        uscs, depth, NSPT, FC_fill, FC_no_fill, unit_weight, D50 = df['Soil type (USCS)'].fillna(''), df[
            'Depth (m)'].astype(float), df.NSPT.astype(float), df['FC (%)'].fillna(method='ffill').fillna(
            method='bfill').astype(float), df['FC (%)'].astype(float), df['Unit weight (kN/m3)'].astype(float), df[
                                                                       'D50 (mm)']
        FC = np.where(soil_type == 'Sand', FC_fill, FC_no_fill)

        depth_max = np.max(depth)
        thickness = np.append(depth[0], np.diff(depth))
        rd = np.clip(np.exp(
            -1.012 - 1.126 * np.sin(depth / 11.73 + 5.133) + Mw * (0.106 + 0.118 * np.sin(depth / 11.28 + 5.142))),
            None, 1)
        total_stress = np.cumsum(thickness * unit_weight)
        effective_stress = total_stress - (np.clip((depth - water_table), 0, None) * Water_unit_weight)
        CSR = 0.65 * PGA * (total_stress / effective_stress) * rd

        Ce = [ER / 60] * rows

        if borehole_diameter < 115.1:
            Cb = 1
        elif borehole_diameter < 150.1:
            Cb = 1.05
        else:
            Cb = 1.15

        Cb = [Cb] * rows

        Cr_bins = np.array([0.75, 0.8, 0.85, 0.95, 1])
        Cr = Cr_bins[np.digitize(depth + rod_extension, [3, 4, 6, 10])]

        Cs = [1] * rows

        flag = (soil_type == 'Sand') & (pd.to_numeric(NSPT, errors='coerce').notnull())
        delta_N160 = np.where(flag, np.exp(1.63 + 9.7 / (FC + 0.01) - (15.7 / (FC + 0.01)) ** 2), np.NaN)

        N60 = NSPT * Ce * Cb * Cr * Cs

        N160CS = []
        for i in range(rows):
            if flag[i]:
                def f(N160CS):
                    return (np.clip(
                        (Pa / effective_stress[i]) ** (0.784 - 0.0768 * np.sqrt(np.clip(N160CS, 1, 46))), None,
                        1.7) - (
                                    N160CS - delta_N160[i]) / N60[i])

                N160CS.append(brentq(f, 1, 200))
            else:
                N160CS.append(np.NaN)
        N160CS = np.asarray(N160CS)
        N160 = N160CS - delta_N160
        Cn = N160 / N60

        MSF = 1 + ((np.clip(1.09 + (N160CS / 31.5) ** 2, None, 2.2) - 1) * (8.64 * np.exp(-Mw / 4) - 1.325))
        K_sigma = np.clip(
            1 - (1 / (18.9 - 2.55 * np.sqrt(np.clip(N160CS, None, 37))) * np.log(effective_stress / Pa)),
            None,
            1.1)
        K_alpha = 1
        flag_FS = (flag & (depth > water_table) & (depth <= anal_depth) & (N160CS < 37.5))
        CRR_75 = np.where(flag_FS,
                          np.exp(N160CS / 14.1 + (N160CS / 126) ** 2 - (N160CS / 23.6) ** 3 + (
                                  N160CS / 25.4) ** 4 - 2.8),
                          2)
        CRR = np.clip(CRR_75 * K_sigma * K_alpha * MSF, None, 2)
        CRR = np.where(np.isnan(CRR), 2, CRR)

        FS_liq = np.clip(CRR / CSR, None, 2)
        mid_depth = depth - thickness / 2
        gamma_lim = np.clip(1.859 * (1.1 - np.sqrt(N160CS / 46)) ** 3, 0, 0.5)
        gamma_lim = np.where(np.isnan(gamma_lim), 0, gamma_lim)
        F_alpha = 0.032 + 0.69 * np.sqrt(np.clip(N160CS, 7, None)) - 0.13 * np.clip(N160CS, 7, None)
        F_alpha = np.where(np.isnan(F_alpha), 0, F_alpha)
        gamma_max = np.where(FS_liq > 2, 0,
                             np.where(FS_liq < F_alpha, gamma_lim,
                                      np.clip(gamma_lim, None,
                                              0.035 * (1 - F_alpha) * (2 - FS_liq) / (FS_liq - F_alpha))))
        xi_v = 1.5 * np.exp(-0.369 * np.sqrt(N160CS)) * np.clip(gamma_max, None, 0.08)
        xi_v = np.where(np.isnan(xi_v), 0, xi_v)
        contiguous_soil_type_thickness = pd.DataFrame(thickness).groupby(
            (soil_type != soil_type.shift()).cumsum()).transform('sum')[0]
        LDIi = gamma_max * thickness
        LDIi[contiguous_soil_type_thickness < LD_min_thick] = 0
        LDI = np.cumsum(LDIi[::-1])[::-1]
        Sv1Di = xi_v * thickness
        Sv1D = np.cumsum(Sv1Di[::-1])[::-1]
        LSNi = 1000 * (xi_v / mid_depth) * thickness
        LSNi = np.where(np.isnan(LSNi), 0, LSNi)
        LSN = np.cumsum(LSNi[::-1])[::-1]
        LPIi = np.where(FS_liq < 1, 1 - FS_liq, 0) * np.where(mid_depth <= 20, 10 - 0.5 * mid_depth, 0) * thickness
        LPI = np.cumsum(LPIi[::-1])[::-1]
        a = next((i for i, v in enumerate(FS_liq) if v < 1), rows - 1)
        H1 = depth[a - 1] if a > 0 else 0
        H2 = np.clip(thickness[FS_liq < 1].sum(), None, 9)

        if calculate_Youd == 'YES':
            R = options_df.loc['R (km) [for Youd (2002) only]'].values[0]
            flag_T15 = (flag & (depth > water_table) & (depth <= anal_depth) & (N160 < 15) & (
                    contiguous_soil_type_thickness > LD_min_thick))
            DH = Youd_2002(site_condition, thickness, flag_T15, FC, D50, slope, Mw, R, H_free_face, L_free_face)
            LD_Youd = LDI * (DH / np.max(LDI))
        else:
            DH = np.NaN
            LD_Youd = ([np.NaN] * rows)
        if export_pdf in ['Borings+SUMMARY', 'Borings'] or export_png in ['Borings+SUMMARY', 'Borings']:
            lyr_subdivision = options_df.loc['Plot soil column layer subdivisions'].values[0]
            single_plot(NSPT, N160, FC, FC_no_fill, CSR, CRR, FS_liq, LPI, LDI, Sv1D, LSN, gamma_max, xi_v, depth,
                        mid_depth, thickness, soil_type, soil_types, uscs, water_table, lyr_subdivision, ID, output_dir,
                        export_pdf, export_png, gen_parameters)
        if calculate_Sr == 'YES':
            Sr_OS_2002 = np.where(FS_liq < 1, effective_stress * (0.03 + 0.0075 * np.clip(N160, None, 20)), np.NaN)
            Sr_Weber_2015 = np.where(FS_liq < 1, np.clip(
                6 / 125 * (np.exp(0.1407 * N160CS + 4.2399 * (effective_stress / Pa) ** 0.120) - 0.43991 * (
                        N160CS ** 1.45 + 0.2 * N160CS * (effective_stress / Pa) ** 2.48 + 41.13)), None,
                0.29 * effective_stress), np.NaN)
            delta_seed = [seed_87_delta_N160[i] for i in np.searchsorted(seed_87_fc, FC_fill + 0.01) - 1]
            N160CS_seed_87 = N160 + delta_seed
            Sr_IB_2008_void_red = np.where(FS_liq < 1, effective_stress * np.exp(
                N160CS_seed_87 / 16 + ((N160CS_seed_87 - 16) / 21.2) ** 3 - 3), np.NaN)
            Sr_weighted = OS_wf * Sr_OS_2002 + IB_wf * Sr_IB_2008_void_red + W_wf * Sr_Weber_2015
        else:
            Sr_OS_2002 = pd.Series([np.NaN] * rows)
            Sr_Weber_2015 = pd.Series([np.NaN] * rows)
            Sr_IB_2008_void_red = pd.Series([np.NaN] * rows)
            Sr_weighted = pd.Series([np.NaN] * rows)
        if calculate_fines == 'YES':
            PI = df['PI']
            wc_LL = df['wc'] / df['LL']
            fines_susceptibility = np.where((PI <= 12) & (wc_LL > 0.85), 'Susceptible',
                                            np.where((PI <= 18) & (wc_LL > 0.8), 'Moderate', None))
            wc_LL_PI_all.append([wc_LL.tolist(), PI.tolist()])
            if 'Moderate' in fines_susceptibility and 'Susceptible' not in fines_susceptibility:
                fines_susceptibility_all.append('Moderate')
            elif 'Susceptible' in fines_susceptibility:
                fines_susceptibility_all.append('Susceptible')
            else:
                fines_susceptibility_all.append('Non')
        else:
            fines_susceptibility = pd.Series(None)
            fines_susceptibility_all.append(None)
        print(ID + ' successful')
        depth_all.append(depth_max)
        ws_all.append(ws)
        LPI_all.append(np.max(LPI))
        LSN_all.append(np.max(LSN))
        LDI_all.append(np.max(LDI))
        LD_Youd_all.append(DH)
        Sv1D_all.append(np.max(Sv1D))
        H1_all.append(H1)
        H2_all.append(H2)

        if export_xlsx == 'YES':
            df2 = pd.DataFrame(
                data=(
                    [soil_type, total_stress, effective_stress, Cn, Cr, pd.Series(Cb), pd.Series(Cs), pd.Series(Ce),
                     N60,
                     N160,
                     delta_N160,
                     N160CS, rd, MSF, CSR, CRR,
                     FS_liq, LPI, gamma_lim, F_alpha, gamma_max, xi_v, LDI, LD_Youd, Sv1D, LSN, fines_susceptibility,
                     Sr_OS_2002, Sr_IB_2008_void_red, Sr_Weber_2015, Sr_weighted])).T
            df2.columns = ['Soil type', 'Total stress (kPa)', 'Effective stress (kPa)', 'Cn', 'Cr', 'Cb', 'Cs', 'Ce',
                           'N60',
                           'N160',
                           'Δ_N160', 'N160CS', 'rd', 'MSF', 'CSR', 'CRR', 'FS_liq', 'LPI', 'γ_lim', 'F_α',
                           'γ_max', 'ξ_v', 'LDI (m)', 'LD (m) Youd (2002)', 'Sv-1D (m)', 'LSN', 'Fines susceptibility',
                           'Sr (kPa) Olson & Stark (2002)', 'Sr (kPa) Idriss & Boulanger (2008)',
                           'Sr (kPa) Weber (2015)', 'Sr (kPa)']
            df = df.join([df2])
            df.to_excel(writer, index=False, sheet_name=ws, float_format='%0.3f')

    if export_pdf in ['Borings+SUMMARY', 'SUMMARY'] or export_png in ['Borings+SUMMARY', 'SUMMARY']:
        SUMMARY_1(LPI_all, LSN_all, LDI_all, LD_Youd_all, Sv1D_all, ws_all, calculate_Youd, LD_min_thick, output_dir,
                  export_pdf, export_png, gen_parameters)
        SUMMARY_2(H1_all, H2_all, ws_all, output_dir, export_pdf, export_png, gen_parameters)
        if calculate_fines == 'YES':
            SUMMARY_3(wc_LL_PI_all, ws_all, output_dir, export_pdf, export_png)

    if export_xlsx == 'YES':
        LPI_ranges = [0, 5, 15]
        LPI_labels = ['Low', 'Moderate', 'High']
        LSN_ranges = [0, 10, 20, 30, 40, 50]
        LSN_labels = ['Little to non', 'Minor', 'Moderate', 'Moderate to severe', 'Mayor', 'Severe damage']
        LPI_results_labels = [LPI_labels[i] for i in np.searchsorted(LPI_ranges, LPI_all) - 1]
        LSN_results_labels = [LSN_labels[i] for i in np.searchsorted(LSN_ranges, LSN_all) - 1]
        df2 = pd.DataFrame(data=(
            [depth_all, subdivided_wt, LPI_all, LPI_results_labels, LSN_all, LSN_results_labels, LDI_all, LD_Youd_all,
             Sv1D_all,
             H1_all, H2_all, fines_susceptibility_all])).T
        df2.columns = ['Total depth (m)', 'Subdivided @water_table', 'LPI', 'LPI label', 'LSN', 'LSN label', 'LDI (m)',
                       'LD (m) Youd (2002)',
                       'Sv-1D (m)', 'H1 (m)', 'H2 (m)', 'Fines susceptibility']
        df = input_df[input_df['run'] == 'YES'].reset_index().join([df2])
        df.to_excel(writer, index=False, sheet_name='SUMMARY', float_format='%0.3f')
        writer.save()

    input('Press ENTER to exit')
