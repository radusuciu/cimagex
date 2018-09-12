from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import pathlib
import csv





def heatmap(ratios, x_labels, y_labels, output_path):
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': 'Arial',
        # 'svg.fonttype': 'svgfont'
    })

    fig = plt.gcf()
    # ax = plt.subplot(131)
    ax = plt.subplot()
    fig.set_size_inches(6, 76, forward=True)

    cmap = cm.get_cmap('Blues', lut=4)
    cmap.set_under('lightgray')
    bounds = [0, 0.1, 2, 4, 10, 20]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    im = ax.imshow(ratios, interpolation='nearest', cmap=cmap, norm=norm, vmin=0.1, vmax=20, aspect='auto')

    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels, size=10)

    ax.get_xaxis().tick_top()
    plt.xticks(rotation=60, ha='left')
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)

    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.02, boundaries=[0, 0.1, 2, 4, 10, 20])
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_ylabel('R', rotation='horizontal', va='center')
    cbar.ax.get_yaxis().set_label_coords(-0.1, 0.35)
    cbar.set_ticks([0, 1, 3, 7, 15])
    cbar.set_ticklabels(['ND', '<2', '2-4', '4-10', '>10'])

    plt.tight_layout()

    plt.savefig(str(pathlib.Path(output_path).with_suffix('.svg')), format='svg')
    plt.close()

    with pathlib.Path(output_path).with_suffix('.csv').open('w') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow([''] + x_labels)

        for label, data in zip(y_labels, ratios):
            writer.writerow([label] + data)
