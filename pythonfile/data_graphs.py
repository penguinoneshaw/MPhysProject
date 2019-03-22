import matplotlib.dates as mdates
from pathlib import Path

import statsmodels.api as sm
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import pprint
matplotlib.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,
    "errorbar.capsize": 0.5,
    "pgf.preamble": [
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{mathpazo}",
        r"\usepackage[version-1-compatibility]{siunitx}"
    ]
}
matplotlib.rcParams.update(pgf_with_pdflatex)

import matplotlib.pyplot as plt

ANALYSES = [
    {
        "folder": "UNESCO-10",
        "t": 10,
        "title": "SOFAR minima using UNESCO algorithm",
        "full": True,
        "temp": "Temperature Difference between Surface and Calculated SOFAR Axis",
        "ps": "Normalised Power Spectrum using UNESCO algorithm",
        "sst": "Sea Surface Temperature"
    },
    {
        "folder": "LEROY-10",
        "t": 10,
        "title": "SOFAR minima using Leroy et al. equation",
        "full": True,
        "temp": False,
        "ps": False,
        "sst": False
    },
    {
        "folder": "UNESCO-10-REJECT",
        "title": "UNESCO algorithm rejecting shallower profiles than 1km",
        "t": 10,
        "temp": False,
        "full": True,
        "ps": False,
        "sst": False
    },
    {
        "folder": "UNESCO-30",
        "t": 30,
        "title": "UNESCO algorithm with a 30 day averaging",
        "temp": False,
        "full": True,
        "ps": False,
        "sst": False
    },
    {
        "folder": "UNESCO-IDEAL",
        "t": 30,
        "title": "UNESCO algorithm using the 'ideal' fitting function.",
        "temp": False,
        "full": True,
        "ps": False,
        "sst": False
    }
]

regions = list(range(0, 8))

locations = pd.read_csv(Path('locations.csv'), dtype={'mesh': 'str'})
"""
for a in ANALYSES:
    fig, plots = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                              'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))

    fig.text(
        0.02, 0.5, r"Distance below sea level (m)", va='center', rotation='vertical')
    fig.text(
        0.5, 0.02, r"Year", ha='center')
    for i in reversed(regions):
        region = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.sos-results.csv"), parse_dates=[
            'date'],  index_col='date')
        region['speed of sound minimum depth (m)'].plot(ax=plots.T.flat[i], yerr=region['err'], ylim=[
            1500, 0], ecolor='red', marker='+', barsabove=True, elinewidth=0.2, linestyle='None')
        region[region.columns[1]].rolling(
            12, win_type='hamming').mean().plot(ax=plots.T.flat[i], style="k")

        plots.T.flat[i].set_title(f"Mesh Region {i}")
        plots.T.flat[i].set_xlabel('')
        plots.T.flat[i].set_xlim("1990-01-01", "2019-01-01")
        plots.T.flat[i].xaxis.set_minor_locator(mdates.YearLocator())
        fig.suptitle(a['title'])
        plots.T.flat[i].grid()
    fig.savefig(Path(f"final_output/figures/{a['folder']}-toplevel.pdf"))
    plt.close(fig)

    if a['full']:
        fig3, plots3 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
            'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
        fig6, plots6 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
            'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
        fig3.text(
            0.02, 0.5, r"Distance below sea level (m)", va='center', rotation='vertical')
        fig6.text(
            0.02, 0.5, r"Distance below sea level (m)", va='center', rotation='vertical')
        fig3.text(
            0.5, 0.02, r"Year", ha='center')
        fig6.text(
            0.5, 0.02, r"Year", ha='center')

        for j in locations.index:
            i = locations['mesh'][j]
            three = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.3.sos-results.csv"), parse_dates=[
                'date'],  index_col='date')
            three['speed of sound minimum depth (m)'].plot(ax=plots3.flat[j], yerr=three['err'], ylim=[
                1500, 0], ecolor='red', marker='+', elinewidth=0.2, barsabove=True, linestyle='None')
            three[three.columns[1]].rolling(
                12, win_type='hamming').mean().plot(ax=plots3.flat[j], style="k")

            plots3.flat[j].set_title(f"{i[:3]}: {locations['name'][j]}")
            plots3.flat[j].set_xlabel('')
            fig3.suptitle(a['title'])
            plots3.flat[j].grid()

            try:
                six = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.6.sos-results.csv"), parse_dates=[
                    'date'],  index_col='date')
                six['speed of sound minimum depth (m)'].plot(ax=plots6.flat[j], yerr=six['err'], ylim=[
                    1500, 0], ecolor='red', marker='+', elinewidth=0.2, barsabove=True, linestyle='None')
                six[six.columns[1]].rolling(
                    12, win_type='hamming').mean().plot(ax=plots6.flat[j], style="k")
                plots6.flat[j].set_title(f"{i}: {locations['name'][j]}")
                plots6.flat[j].set_xlabel('')
                fig6.suptitle(a['title'])
                plots6.flat[j].grid()
            except FileNotFoundError as f:
                # fig6.delaxes(plots6[j])
                continue
        fig3.savefig(
            Path(f"final_output/figures/{a['folder']}-threelevel.pdf"))
        fig6.savefig(Path(f"final_output/figures/{a['folder']}-sixlevel.pdf"))

        plt.close(fig3)
        plt.close(fig6)

    if a['temp']:
        fig, plots = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
            'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
        fig.text(
            0.02, 0.5, r"Temperature Difference (${}^\circ$C)", va='center', rotation='vertical')
        fig.text(
            0.5, 0.02, r"Year", ha='center')
        for i in reversed(regions):
            temp = pd.read_csv(Path(f"final_output/{ a['folder'] }/temp/{i}.temp-results.csv"), parse_dates=[
                'date'],  index_col='date')
            temp['average temperature across profiles'].plot(
                ax=plots.T.flat[i], marker='+', linestyle='None')
            temp['average temperature across profiles'].rolling(
                12, win_type='hamming').mean().plot(ax=plots.T.flat[i], style="k")

            plots.T.flat[i].set_title(f"Mesh Region {i}")
            plots.T.flat[i].set_xlabel('')
            plots.T.flat[i].xaxis.set_minor_locator(mdates.YearLocator())
            plots.T.flat[i].set_xlim("1990-01-01", "2019-01-01")
            fig.suptitle(a['temp'])
            plots.T.flat[i].grid()
        fig.savefig(
            Path(f"final_output/figures/{a['folder']}-temp-toplevel.pdf"))
        if a['full']:
            fig3, plots3 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig6, plots6 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig3.text(
                0.02, 0.5, r"Temperature Difference (${}^\circ$C)", va='center', rotation='vertical')
            fig6.text(
                0.02, 0.5, r"Temperature Difference (${}^\circ$C)", va='center', rotation='vertical')
            fig3.text(
                0.5, 0.02, r"Year", ha='center')
            fig6.text(
                0.5, 0.02, r"Year", ha='center')
            fig6.suptitle(a['temp'])
            fig3.suptitle(a['temp'])

            for j in locations.index:
                i = locations['mesh'][j]
                three = pd.read_csv(Path(f"final_output/{ a['folder'] }/temp/{i}.3.temp-results.csv"), parse_dates=[
                    'date'],  index_col='date')
                three['average temperature across profiles'].plot(
                    ax=plots3.flat[j],  marker='+',  linestyle='None')
                three['average temperature across profiles'].rolling(
                    12, win_type='hamming').mean().plot(ax=plots3.flat[j], style="k")

                plots3.flat[j].set_title(f"{i[:3]}: {locations['name'][j]}")
                plots3.flat[j].set_xlabel('')
                plots3.flat[j].grid()

                try:
                    six = pd.read_csv(Path(f"final_output/{ a['folder'] }/temp/{i}.6.temp-results.csv"), parse_dates=[
                        'date'],  index_col='date')
                    six['average temperature across profiles'].plot(
                        ax=plots6.flat[j], marker='+', linestyle='None')
                    six['average temperature across profiles'].rolling(
                        12, win_type='hamming').mean().plot(ax=plots6.flat[j], style="k")
                    plots6.flat[j].set_title(f"{i}: {locations['name'][j]}")
                    plots6.flat[j].set_xlabel('')
                    plots6.flat[j].grid()
                except FileNotFoundError as f:
                    fig6.delaxes(plots6[j])
                    continue
            fig3.savefig(
                Path(f"final_output/figures/{a['folder']}-temp-threelevel.pdf"))
            fig6.savefig(
                Path(f"final_output/figures/{a['folder']}-temp-sixlevel.pdf"))
            plt.close(fig3)
            plt.close(fig6)
    if a['ps']:
        fig, plots = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
            'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
        fig.text(
            0.02, 0.5, r"Temperature Difference (${}^\circ$C)", va='center', rotation='vertical')
        fig.text(
            0.5, 0.02, r"Frequency (10 days${}^{-1})$", ha='center')
        for i in reversed(regions):
            temp = pd.read_csv(
                Path(f"final_output/{ a['folder'] }/power_spectra/{i}.ps-results.csv"))
            temp.plot(x='frequency', y='Normalised Power', ax=plots.T.flat[i])

            plots.T.flat[i].set_title(f"Mesh Region {i}")
            plots.T.flat[i].set_xlabel('')
            fig.suptitle(a['ps'])
            plots.T.flat[i].grid()
        fig.savefig(
            Path(f"final_output/figures/{a['folder']}-ps-toplevel.pdf"))
        if a['full']:
            fig3, plots3 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig6, plots6 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig3.text(
                0.02, 0.5, r"Normalised Power", va='center', rotation='vertical')
            fig6.text(
                0.02, 0.5, r"Normalised Power", va='center', rotation='vertical')
            fig3.text(
                0.5, 0.02, r"Frequency (10 days${}^{-1})$", ha='center')
            fig6.text(
                0.5, 0.02, r"Frequency (10 days${}^{-1})$", ha='center')
            fig6.suptitle(a['ps'])
            fig3.suptitle(a['ps'])

            for j in locations.index:
                i = locations['mesh'][j]
                three = pd.read_csv(
                    Path(f"final_output/{ a['folder'] }/power_spectra/{i}.3.ps-results.csv"))
                three.plot(
                    x='frequency', y='Normalised Power', ax=plots3.flat[j])
                plots3.flat[j].set_title(f"{i[:3]}: {locations['name'][j]}")
                plots3.flat[j].set_xlabel('')
                plots3.flat[j].grid()

                try:
                    six = pd.read_csv(
                        Path(f"final_output/{ a['folder'] }/power_spectra/{i}.6.ps-results.csv"))
                    six.plot(
                        x='frequency', y='Normalised Power', ax=plots6.flat[j])
                    plots6.flat[j].set_title(f"{i}: {locations['name'][j]}")
                    plots6.flat[j].set_xlabel('')
                    plots6.flat[j].grid()
                except FileNotFoundError as f:
                    fig6.delaxes(plots6[j])
                    continue
            fig3.savefig(
                Path(f"final_output/figures/{a['folder']}-ps-threelevel.pdf"))
            fig6.savefig(
                Path(f"final_output/figures/{a['folder']}-ps-sixlevel.pdf"))
    if a['sst']:
        fig, plots = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
            'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
        fig.text(
            0.02, 0.5, r"Sea Surface Temperature (${}^\circ$C)", va='center', rotation='vertical')
        fig.text(
            0.5, 0.02, r"Year", ha='center')
        for i in reversed(regions):
            temp = pd.read_csv(Path(f"final_output/{ a['folder'] }/sst/{i}.temp-results.csv"), parse_dates=[
                'date'],  index_col='date')
            temp['sea surface temperatures'].plot(
                ax=plots.T.flat[i], marker='+', linestyle='None')
            temp['sea surface temperatures'].rolling(
                12, win_type='hamming').mean().plot(ax=plots.T.flat[i], style="k")

            plots.T.flat[i].set_title(f"Mesh Region {i}")
            plots.T.flat[i].set_xlabel('')
            plots.T.flat[i].xaxis.set_minor_locator(mdates.YearLocator())
            plots.T.flat[i].set_xlim("1990-01-01", "2019-01-01")
            fig.suptitle(a['sst'])
            plots.T.flat[i].grid()
        fig.savefig(
            Path(f"final_output/figures/{a['folder']}-sst-toplevel.pdf"))
        if a['full']:
            fig3, plots3 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig6, plots6 = plt.subplots(4, 2, 'col', 'row',  True, gridspec_kw={
                'hspace': 0.3, 'bottom': 0.08, 'top': 0.92}, figsize=(6, 7))
            fig3.text(
                0.02, 0.5, r"Sea Surface Temperature (${}^\circ$C)", va='center', rotation='vertical')
            fig6.text(
                0.02, 0.5, r"Sea Surface Temperature (${}^\circ$C)", va='center', rotation='vertical')
            fig3.text(
                0.5, 0.02, r"Year", ha='center')
            fig6.text(
                0.5, 0.02, r"Year", ha='center')
            fig6.suptitle(a['sst'])
            fig3.suptitle(a['sst'])

            for j in locations.index:
                i = locations['mesh'][j]
                three = pd.read_csv(Path(f"final_output/{ a['folder'] }/sst/{i}.3.temp-results.csv"), parse_dates=[
                    'date'],  index_col='date')
                three['sea surface temperatures'].plot(
                    ax=plots3.flat[j],  marker='+',  linestyle='None')
                three['sea surface temperatures'].rolling(
                    12, win_type='hamming').mean().plot(ax=plots3.flat[j], style="k")

                plots3.flat[j].set_title(f"{i[:3]}: {locations['name'][j]}")
                plots3.flat[j].set_xlabel('')
                plots3.flat[j].grid()

                try:
                    six = pd.read_csv(Path(f"final_output/{ a['folder'] }/sst/{i}.6.temp-results.csv"), parse_dates=[
                        'date'],  index_col='date')
                    six['sea surface temperatures'].plot(
                        ax=plots6.flat[j], marker='+', linestyle='None')
                    six['sea surface temperatures'].rolling(
                        12, win_type='hamming').mean().plot(ax=plots6.flat[j], style="k")
                    plots6.flat[j].set_title(f"{i}: {locations['name'][j]}")
                    plots6.flat[j].set_xlabel('')
                    plots6.flat[j].grid()
                except FileNotFoundError as f:
                    fig6.delaxes(plots6[j])
                    continue
            fig3.savefig(
                Path(f"final_output/figures/{a['folder']}-sst-threelevel.pdf"))
            fig6.savefig(
                Path(f"final_output/figures/{a['folder']}-sst-sixlevel.pdf"))
            plt.close(fig3)
            plt.close(fig6)
# plt.show()
"""
a = ANALYSES[0]
gradients = {}
for i in regions:
    region = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.sos-results.csv"), parse_dates=[
        'date'],  index_col='date')
    X = sm.add_constant(region[region.columns[0]]["2007-01-01":])
    y = region[region.columns[1]]["2007-01-01":]
    eps = region[region.columns[2]]["2007-01-01":]
    mod = sm.WLS(y, X, sigma=eps)
    fmod = mod.fit()
    gradients[i] = (fmod.params[region.columns[0]], fmod.bse[region.columns[0]])

for j in locations.index:
    i = locations['mesh'][j]
    three = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.3.sos-results.csv"), parse_dates=[
        'date'],  index_col='date')

    X = sm.add_constant(three[three.columns[0]]["2007-01-01":])
    y = three[three.columns[1]]["2007-01-01":]
    eps = three[three.columns[2]]["2007-01-01":]
    mod = sm.WLS(y, X, sigma=eps)
    fmod = mod.fit()
    gradients[i[:3]] = (fmod.params[region.columns[0]], fmod.bse[region.columns[0]])

    six = pd.read_csv(Path(f"final_output/{ a['folder'] }/speed_of_sound/{i}.6.sos-results.csv"), parse_dates=[
        'date'],  index_col='date')
    X = sm.add_constant(six[six.columns[0]]["2007-01-01":])
    y = six[six.columns[1]]["2007-01-01":]
    eps = six[six.columns[2]]["2007-01-01":]
    mod = sm.WLS(y, X, sigma=eps)
    fmod = mod.fit()
    gradients[i] = (fmod.params[region.columns[0]], fmod.bse[region.columns[0]])
pprint.pprint(gradients)