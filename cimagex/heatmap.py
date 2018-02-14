from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib
import pathlib
import csv


def heatmap(ratios, x_labels, y_labels, output_path):
    # matplotlib.rcParams.update({
    #     'font.family': 'sans-serif',
    #     'font.sans-serif': 'Univers LT Std 55',
    #     'svg.fonttype': 'svgfont'
    # })

    fig = plt.gcf()
    ax = plt.subplot(131)

    fig.set_size_inches(9, 9, forward=True)
    im = ax.imshow(ratios, interpolation='nearest', cmap=cm.Blues, vmin=0, vmax=20, aspect='auto')

    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels, size=10)

    ax.get_xaxis().tick_top()
    plt.xticks(rotation=90)
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)

    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.02, ticks=[0, 5, 10, 15, 20])
    cbar.ax.set_ylabel('R', rotation='horizontal', va='center')
    cbar.ax.get_yaxis().set_label_coords(-0.1, 0.35)

    plt.tight_layout()

    plt.savefig(str(pathlib.Path(output_path).with_suffix('.svg')), format='svg')

    with pathlib.Path(output_path).with_suffix('.csv').open('w', encoding='UTF-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow([''] + x_labels)

        for label, data in zip(y_labels, ratios):
            writer.writerow([label] + data)
