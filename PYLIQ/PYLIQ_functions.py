import matplotlib.pyplot as plt
import numpy as np
import itertools
import pandas as pd
from matplotlib.ticker import AutoMinorLocator
from math import log10


def single_plot(NSPT, N160, FC, FC_no_fill, CSR, CRR, FS_liq, LPI, LDI, Sv1D, LSN, gamma_max, xi_v, depth,
                mid_depth, thickness, soil_type, soil_types, uscs, water_table, lyr_subdivision, ID, output_dir,
                export_pdf, export_png, gen_parameters):
    depth_max_rounded = np.ceil(np.max(depth))
    hatches = ['...', '', '//', '|', '|', '///', '', '']
    colors = ('white', 'orange', 'white', 'white', 'white', 'white', 'white', 'white')
    edgecolors = ('grey', 'black', 'green', 'orange', 'black', 'black', 'black', 'black')
    primary_x_vars = [[NSPT, N160], [FC, FC_no_fill], [CSR, CRR], [FS_liq], [LPI], [LDI], [Sv1D]]
    primary_x_vars_labels = [['NSPT', 'N160'], ['FC filled', 'FC'], ['CSR', 'CRR'], ['FS_liq'], ['LPI'], ['LDI'],
                             ['Sv1D']]
    primary_x_vars_styles = [['ko-', 'bo'], ['bo', 'ko'], ['k-', 'ko-'], ['ko-'], ['ko-'], ['ko-'], ['ko-']]
    primary_x_labels = ['Blow count', 'Fines content (%)', 'CSR & CRR', 'FS Liq', 'LPI', 'LDI (m)', 'Sv-1D (m)']
    primary_y_vars = [depth, depth, depth, depth, mid_depth, mid_depth, mid_depth]
    secondary_x_vars = [None, None, None, None, LSN, gamma_max, xi_v, None]
    secondary_x_vars_labels = [None, None, None, None, 'LSN', 'γ max', 'ξv', None]
    contiguous_soil_type = soil_type.groupby(
        (uscs != uscs.shift()).cumsum()).first()
    contiguous_uscs = uscs.groupby(
        (uscs != uscs.shift()).cumsum()).first()
    contiguous_uscs_thickness = pd.DataFrame(thickness).groupby(
        (uscs != uscs.shift()).cumsum())[0].sum()
    fig, axes = plt.subplots(1, 8, sharey=True)
    for i, ax in enumerate(axes.flat):
        if i < 7:
            for j, var in enumerate(primary_x_vars[i]):
                ax.plot(var, primary_y_vars[i], primary_x_vars_styles[i][j], label=primary_x_vars_labels[i][j])
            ax.legend() if i < 3 else None
            ax.set_xlim(left=0)
            ax.set(xlabel=primary_x_labels[i])
            ax.xaxis.set_minor_locator(AutoMinorLocator())
        else:
            soil_type_indices = [soil_types.index(i) for i in contiguous_soil_type]
            hatch = [hatches[i] for i in soil_type_indices]
            color = [colors[i] for i in soil_type_indices]
            edgecolor = [edgecolors[i] for i in soil_type_indices]
            bottom = np.append(0, np.cumsum(contiguous_uscs_thickness)[:-1])
            patches = ax.bar(1, contiguous_uscs_thickness, bottom=bottom, color=color,
                             lw=(0 if lyr_subdivision == 'NO' else 1),
                             edgecolor=edgecolor)
            for patch_set, hatch in zip(patches, hatch):
                plt.setp(patch_set, hatch=hatch)
            cut_z = 0
            for ri, r in enumerate(patches):
                cut_z = cut_z + contiguous_uscs_thickness[ri + 1]
                cut_z_mid = cut_z - contiguous_uscs_thickness[ri + 1] / 2
                ax.text(r.get_x() + r.get_width() / 2, cut_z_mid, str(contiguous_uscs[ri + 1]), ha="center",
                        va="center", color="black", fontsize=10, fontweight='bold')
        if secondary_x_vars[i] is not None:
            ax_secondary = ax.twiny()
            ax_secondary.plot(secondary_x_vars[i], primary_y_vars[i], 'o-', color='grey',
                              label=secondary_x_vars_labels[i])
            ax_secondary.set_xlabel(secondary_x_vars_labels[i], color='grey')
            ax_secondary.tick_params(axis='x', which='both', colors='grey')
            ax_secondary.xaxis.set(label_position='bottom', ticks_position='bottom')
            ax_secondary.set_xlim(left=0)
            ax_secondary.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
    axes[0].yaxis.set_minor_locator(AutoMinorLocator())
    axes[0].set(ylabel='Depth (m)')
    axes[3].axvline(x=1, color='r')
    axes[7].set(title='Soil column', xticks=[])
    axes[7].axhline(y=water_table, color='b')
    axes[7].spines['top'].set_visible(False)
    axes[7].spines['right'].set_visible(True)
    axes[7].spines['bottom'].set_visible(False)
    axes[7].tick_params(axis='y', which='both', labelright=True, right=True)
    plt.ylim(bottom=0, top=depth_max_rounded)
    axes[0].invert_yaxis()
    fig.set_size_inches(16, 10)
    plt.suptitle('Soil liquefaction analysis from SPT ' + ID, fontsize=14)
    fig.text(0.1, 0.965, gen_parameters, color='blue', fontsize=12)
    fig.text(1, 0, 'PYLIQ', color='grey', fontsize=12)
    plt.subplots_adjust(left=None, bottom=0.03, right=1, top=0.9, wspace=0.1, hspace=None)
    if export_pdf in ['Borings+SUMMARY', 'Borings']: plt.savefig(output_dir + ID + '.pdf',
                                                                 format='pdf', dpi=900,
                                                                 bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'Borings']: plt.savefig(output_dir + ID + '.png',
                                                                 format='png', bbox_inches='tight')
    plt.close()


def Youd_2002(site_condition, thickness, flag_T15, FC, D50, slope, Mw, R, H_free_face, L_free_face):
    if site_condition == 'Level ground':
        DH = 0
    else:
        T15 = thickness[flag_T15].sum()
        if T15 != 0:
            F15 = FC[flag_T15].mean()
            D50_filled = D50.fillna(method='ffill').fillna(method='bfill')
            D5015 = D50_filled[flag_T15].mean()
            if site_condition == 'Sloping ground':
                if 0.1 <= slope <= 6:
                    DH = 10 ** (-16.213 + 1.532 * Mw - 1.406 * log10(
                        10 ** (0.89 * Mw - 5.64) + R) - 0.012 * R + 0.338 * log10(
                        slope) + 0.54 * log10(T15) + 3.413 * log10(100 - F15) - 0.795 * log10(D5015 + 0.1))
                else:
                    print(ID + ' Youd (2002): Ground slope must be between 0.1 to 6 %')
                    DH = np.NaN
            if site_condition == 'Free face':
                if 1 <= 100 * (H_free_face / L_free_face) <= 20:
                    DH = 10 ** (-16.713 + 1.532 * Mw - 1.406 * log10(
                        10 ** (0.89 * Mw - 5.64) + R) - 0.012 * R + 0.592 * log10(
                        100 * (H_free_face / L_free_face)) + 0.54 * log10(T15) + 3.413 * log10(
                        100 - F15) - 0.795 * log10(
                        D5015 + 0.1))
                else:
                    print(ID + ' Youd (2002): Free face ratio (H/L) must be between 1 to 20 %')
                    DH = np.NaN
        else:
            DH = 0
    return DH


def SUMMARY_1(LPI_all, LSN_all, LDI_all, LD_Youd_all, Sv1D_all, ws_all, calculate_Youd, LD_min_thick, output_dir,
              export_pdf, export_png, gen_parameters):
    LPI_ranges = [0, 5, 15]
    LPI_colors = ['yellow', 'orange', 'red']
    LPI_labels = ['Low', 'Moderate', 'High']
    LSN_ranges = [0, 10, 20, 30, 40, 50]
    LSN_colors = ['slateblue', 'blue', 'gold', 'yellow', 'orange', 'red']
    LSN_labels = ['Little to non', 'Minor', 'Moderate', 'Moderate to severe', 'Mayor', 'Severe damage']
    fig, ax = plt.subplots(2, 2)
    fig.set_size_inches(16, 10)
    fig.text(0.97, 0.96, 'PYLIQ', color='grey', fontsize=12)
    fig.text(0.1, 0.965, gen_parameters, color='blue', fontsize=12)
    fig.suptitle('Soil liquefaction analysis from SPT SUMMARY')
    primary_y_vars = [LPI_all, LSN_all]
    ranges = [LPI_ranges, LSN_ranges]
    labels = [LPI_labels, LSN_labels]
    colors = [LPI_colors, LSN_colors]
    titles = ['LPI from Iwasaki et al. (1978)', 'LSN from van Ballegooy et al. (2014)']
    for i, var in enumerate(primary_y_vars):
        ax[0, i].bar(ws_all, var, color='black', zorder=2)
        ax[0, i].set(title=titles[i])
        bottom, top = ax[0, i].get_ylim()
        ranges[i].append(top)
        for j, label in enumerate(labels[i]):
            ax[0, i].axhspan(ranges[i][j], ranges[i][j + 1], color=colors[i][j], zorder=1, label=label)
        ax[0, i].set(ylim=[0, top])
        ax[0, i].legend()
    bar_width = 0.35
    ws_all_indices = np.arange(len(ws_all))
    if calculate_Youd == 'YES':
        ax[1, 0].bar(ws_all_indices - bar_width / 2, LDI_all, bar_width, color='black', label='LDI I&B (2014)')
        ax[1, 0].bar(ws_all_indices + bar_width / 2, LD_Youd_all, bar_width, color='red', label='LD Youd (2002)')
    else:
        ax[1, 0].bar(ws_all_indices, LDI_all, color='black', label='LDI I&B (2014)')
    ax[1, 0].set_xticks(ws_all_indices)
    ax[1, 0].legend()
    ax[1, 0].set(title='LD (m), minimum thickness considered (m): ' + str(LD_min_thick))
    ax[1, 1].bar(ws_all, Sv1D_all, color='black')
    ax[1, 1].set(title='Sv-1D (m)')
    plt.subplots_adjust(left=None, bottom=0.11, right=1, top=0.92, wspace=0.1, hspace=0.25)
    for i, axi in enumerate(ax.flat):
        axi.set_xticklabels(ws_all, rotation=30, ha='right')
    if export_pdf in ['Borings+SUMMARY', 'SUMMARY']: plt.savefig(output_dir + 'SUMMARY_1.pdf',
                                                                 format='pdf', dpi=900,
                                                                 bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'SUMMARY']: plt.savefig(output_dir + 'SUMMARY_1.png',
                                                                 format='png', bbox_inches='tight')
    plt.close()


def SUMMARY_2(H1_all, H2_all, ws_all, output_dir, export_pdf, export_png, gen_parameters):
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16, 10)
    fig.text(0.87, 0.065, 'PYLIQ', color='grey', fontsize=12)
    fig.text(0.1, 0.89, gen_parameters, color='blue', fontsize=12)
    plt.title('Site identification of liquefaction-induced damage from Ishihara (1985)')
    H1_200gal = [0.833, 1.077, 1.313, 1.6, 1.881, 2.213, 2.449, 2.59, 2.738, 2.908, 2.983, 3., 3.]
    H2_200gal = [0.781, 1.022, 1.263, 1.572, 1.871, 2.276, 2.605, 2.856, 3.146, 3.563, 3.951, 4.212, 8.009]
    H1_300gal = [1.245, 1.776, 2.218, 2.814, 3.537, 4.111, 4.789, 5.158, 5.306, 5.498, 5.72, 5.883, 5.994,
                 6.091,
                 6.159,
                 6.219, 6.25, 6.287, 6.318]
    H2_300gal = [0.826, 1.172, 1.461, 1.845, 2.325, 2.72, 3.259, 3.615, 3.789, 4.079, 4.485, 4.911, 5.338,
                 5.832,
                 6.269,
                 6.715, 7.084, 7.433, 7.928]
    H1_400_500gal = [1.09, 0.994, 1.237, 1.436, 1.737, 2.179, 2.68, 3.003, 3.563, 4.078, 4.91, 5.536, 6.169,
                     6.949,
                     7.244,
                     7.583, 7.834, 8.085, 8.292, 8.477, 8.662, 8.825, 8.989, 9.145, 9.272]
    H2_400_500gal = [0.42, 0.382, 0.477, 0.553, 0.658, 0.83, 1.021, 1.135, 1.354, 1.545, 1.879, 2.079, 2.269,
                     2.574,
                     2.727,
                     3.007, 3.306, 3.654, 4.002, 4.428, 4.912, 5.377, 5.919, 6.52, 7.053]
    ax.plot(H1_200gal, H2_200gal, 'k', label='amax≈0.2g')
    ax.plot(H1_300gal, H2_300gal, 'k--', label='amax≈0.3g')
    ax.plot(H1_400_500gal, H2_400_500gal, 'k-.', label='amax≈0.4g-0.5g')
    ax.scatter(H1_all, H2_all, label='borings')
    for i, txt in enumerate(ws_all):
        ann = ax.annotate(' ' + txt, (H1_all[i], H2_all[i]))
        ann.set_clip_box(ax.bbox)
    ax.set(xlabel='Thickness of surface layer, H1 (m)', ylabel='Thickness of liquefiable sand layer, H2 (m)')
    ax.legend()
    plt.xlim(plt.xlim()[0], 10)
    if export_pdf in ['Borings+SUMMARY', 'SUMMARY']: fig.savefig(output_dir + 'SUMMARY_2.pdf',
                                                                 format='pdf', dpi=900,
                                                                 bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'SUMMARY']: fig.savefig(output_dir + 'SUMMARY_2.png',
                                                                 format='png', bbox_inches='tight')
    plt.close()


def SUMMARY_3(wc_LL_PI_all, ws_all, output_dir, export_pdf, export_png):
    markers = itertools.cycle(('o', 's', '^', 'D'))
    fig, ax = plt.subplots(1, 1)
    fig.set_size_inches(16, 10)
    fig.text(0.87, 0.065, 'PYLIQ', color='grey', fontsize=12)
    plt.title('Assessment of the Liquefaction Susceptibility of Fine-Grained Soils from Bray & Sancio (2006)')
    for i, wc_LL_PI in enumerate(wc_LL_PI_all):
        ax.plot(*wc_LL_PI, marker=next(itertools.cycle(markers)), linestyle='None', label=ws_all[i])
    plt.xlim(left=0.4, right=1.6)
    plt.ylim(bottom=0, top=60)
    plt.plot([0.8, 0.8, 1.6], [0, 18, 18], 'k', lw=1.0)
    plt.text(1.2, 9, 'Susceptible', ha='center')
    plt.plot([0.85, 0.85, 1.6], [0, 12, 12], 'k', lw=1.0)
    plt.text(1.2, 15, 'Mod. susceptible', ha='center')
    plt.text(0.5, 45, 'Not susceptible', ha='center')
    plt.ylabel('Plasticity index, PI')
    plt.xlabel('water content / Liquid Limit, wc/LL')
    ax.legend()
    if export_pdf in ['Borings+SUMMARY', 'SUMMARY']: fig.savefig(output_dir + 'SUMMARY_3.pdf',
                                                                 format='pdf', dpi=900,
                                                                 bbox_inches='tight')
    if export_png in ['Borings+SUMMARY', 'SUMMARY']: fig.savefig(output_dir + 'SUMMARY_3.png',
                                                                 format='png', bbox_inches='tight')
    plt.close()
