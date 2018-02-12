from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import pathlib
import csv


def heatmap(ratios, x_labels, y_labels, output_path):
    fig, ax = plt.subplots()
    fig.set_size_inches(6, 15, forward=True)
    im = ax.imshow(ratios, interpolation='nearest',  cmap=cm.Blues)

    ax.set_yticks(range(len(y_labels)))
    ax.set_yticklabels(y_labels, family='sans-serif', size=10)

    ax.get_xaxis().tick_top()
    plt.xticks(rotation=90)
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)

    plt.tight_layout()
    plt.savefig(str(pathlib.Path(output_path).with_suffix('.svg')), format='svg')

    with pathlib.Path(output_path).with_suffix('.csv').open('w', encoding='UTF-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow([''] + x_labels)

        for label, data in zip(y_labels, ratios):
            writer.writerow([label] + data)
