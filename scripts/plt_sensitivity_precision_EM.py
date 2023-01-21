import os
import pathlib
import sys

import seaborn
import pandas
import seaborn
from matplotlib import pyplot as plt

#%% Input and output
# minread0_tsv_path = "min_read_count_0_10_20_50/control_counts_minread_count0.tsv"
# minread10_tsv_path = "min_read_count_0_10_20_50/control_counts_minread_count10.tsv"
# minread20_tsv_path = "min_read_count_0_10_20_50/control_counts_minread_count20.tsv"
# minread60_tsv_path = "min_read_count_0_10_20_50/control_counts_minread_count60.tsv"

precision_bat_png_path = "precision_bat.png"
precision_fish_png_path = "precision_fish.png"
precision_shark_png_path = "precision_shark.png"

minread0_tsv_path = sys.argv[1]
minread10_tsv_path = sys.argv[2]
minread40_tsv_path = sys.argv[3]
minread60_tsv_path = sys.argv[4]
outdir = sys.argv[5]
pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

#%%
df0 = pandas.read_csv(minread0_tsv_path, sep="\t")
df0['Min. read count'] = 0
df10 = pandas.read_csv(minread10_tsv_path, sep="\t")
df10['Min. read count'] = 10
df40 = pandas.read_csv(minread40_tsv_path, sep="\t")
df40['Min. read count'] = 40
df60 = pandas.read_csv(minread60_tsv_path, sep="\t")
df60['Min. read count'] = 60

#%%
cat_df = pandas.concat([df0, df10], axis=0)
cat_df = pandas.concat([cat_df, df40], axis=0)
cat_df = pandas.concat([cat_df, df60], axis=0)
cat_df.reset_index(inplace=True)

cat_df.replace('bat', 'Bat', inplace=True)
cat_df.replace('fish', 'Fish', inplace=True)
cat_df.replace('shark', 'Shark', inplace=True)

#%%
cat_df.rename({'precision (TP/(FP+TP)': 'Precision (TP/(FP+TP)'}, axis=1, inplace=True)
cat_df.rename({'sensitivity (TP/(TP+FN)': 'Sensitivity (TP/(TP+FN)'}, axis=1, inplace=True)
# cat_df.replace('vtam', 'VTAM', inplace=True)
cat_df.replace('dalu', 'DALU', inplace=True)
cat_df.replace('obibar', 'OBIbar', inplace=True)

#%%
vtam_precision_dic={}
vtam_sensitivity_dic={}

vtam_precision_dic['Bat'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Bat") & (cat_df['Min. read count'] == 60), 'Precision (TP/(FP+TP)'].values[0]
vtam_sensitivity_dic['Bat'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Bat") & (cat_df['Min. read count'] == 60), 'Sensitivity (TP/(TP+FN)'].values[0]
vtam_precision_dic['Fish'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Fish") & (cat_df['Min. read count'] == 10), 'Precision (TP/(FP+TP)'].values[0]
vtam_sensitivity_dic['Fish'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Fish") & (cat_df['Min. read count'] == 10), 'Sensitivity (TP/(TP+FN)'].values[0]
vtam_precision_dic['Shark'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Shark") & (cat_df['Min. read count'] == 40), 'Precision (TP/(FP+TP)'].values[0]
vtam_sensitivity_dic['Shark'] = cat_df.loc[(cat_df['Pipeline'] == "vtam") & (cat_df['Dataset'] == "Shark") & (cat_df['Min. read count'] == 40), 'Sensitivity (TP/(TP+FN)'].values[0]

cat_df = cat_df.loc[cat_df['Pipeline'] != "vtam"]

#%%
seaborn.set_theme(style="darkgrid")

datasets = ["Bat", "Fish", "Shark"]

for ds in datasets:
    if ds == "Bat":
        text_label = "VTAM read count cutoff=60"
    elif ds == "Fish":
        text_label = "VTAM read count cutoff=10"
    elif ds == "Shark":
        text_label = "VTAM read count cutoff=40"
    else:
        sys.exit(1)
    plot_df = cat_df.loc[cat_df['Dataset'] == ds, ]
    g = seaborn.catplot(
        data=plot_df, kind="bar", x="Pipeline", y="Precision (TP/(FP+TP)", hue="Min. read count", palette="rocket_r", order=["DALU", "OBIbar"], legend=False)
    g.despine(left=True)
    plt.axhline(vtam_precision_dic[ds], color='blue')
    plt.text(1, vtam_precision_dic[ds]-0.05, text_label, horizontalalignment='right', size='medium', color='blue', weight='semibold')
    plt.title("{} data set".format(ds))
    plt.ylim([0, 1.1])
    plt.legend(loc='upper center',  title='Read count cutoff', framealpha=1, facecolor='white')
    plt.tight_layout()
    preci_png_path = os.path.join(outdir, "precision_{}.png".format(ds))
    plt.savefig(preci_png_path)

    g = seaborn.catplot(
        data=plot_df, kind="bar", x="Pipeline", y="Sensitivity (TP/(TP+FN)", hue="Min. read count", palette="rocket_r", order=["DALU", "OBIbar"], legend=False)
    g.despine(left=True)
    plt.axhline(vtam_sensitivity_dic[ds], color='blue')
    plt.text(1, vtam_sensitivity_dic[ds]+0.05, text_label, horizontalalignment='right', size='medium', color='blue', weight='semibold')
    plt.ylim([0, 1.1])
    plt.title("{} data set".format(ds))
    plt.legend(loc='lower center',  title='Read count cutoff', framealpha=1, facecolor='white')
    plt.tight_layout()
    sensi_png_path = os.path.join(outdir, "sensitivity_{}.png".format(ds))
    plt.savefig(sensi_png_path)
