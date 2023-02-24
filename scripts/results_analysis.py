import csv
from sys import argv
from os import path, pardir
from itertools import groupby, chain, repeat, count
import subprocess

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def get_working_dir():
    return path.join(path.dirname(path.realpath(__file__)), pardir, "out")


def factorials():
    f = 1
    for i in count():
        yield f
        f *= (i + 1)


def plot_with_derivatives(ax, t, h, ders, style, alpha, linewidth=1):
    p = [v / f for v, f in zip(ders, factorials())]
    p.reverse()
    hs = np.linspace(0, h, 10)
    ys = np.polyval(p, hs)
    ts = t + hs
    ax.plot(ts, ys, style, linewidth=linewidth, alpha=alpha)

def plot_piecewise_with_derivatives(ax, steps, values, style, alpha, linewidth=1):
    t = 0
    for h, ders in zip(steps, values):
        plot_with_derivatives(ax, t, h, ders, style, alpha, linewidth)
        t += h

def process_signals(ax, csvFileName, port, style, alpha, linewidth=1):
    with open(csvFileName, 'r') as csvFile:
        csvReader = csv.reader(csvFile, dialect='excel')
        for rowMaxDer in csvReader:
            rowSteps = next(csvReader) # moving reader inside the loop
            rowValues = next(csvReader)

            instanceName, portReference = rowMaxDer[1:3]
            maxDer = int(rowMaxDer[3])
            if port[0] != instanceName or str(port[1]) != portReference or rowValues[0] == 'defect':
                continue # if the port doesn't match go to the next two rows
            steps = list(map(float, rowSteps[3:]))
            vs = list(map(float, rowValues[3:]))
            sampleKeys = chain.from_iterable(repeat(x, maxDer + 1) for x in count())
            samples = zip(sampleKeys, vs)
            values = [tuple(v[1] for v in sample) for _, sample in groupby(samples, lambda s: s[0])]
            plot_piecewise_with_derivatives(ax, steps, values, style, alpha, linewidth)
            ax.set_xlim([0, 50])


def produce_results(figureNum):
    subprocess.call([path.join(get_working_dir(), "Experiments.exe"), str(figureNum)], cwd=get_working_dir())


def experiment_1():
    print('Running Experiment 1.')
    produce_results(1)

    ports = [('tau2omega', 1), ('omega2tau', 1)]
    networks = [
        ('euler_subsystem', 'k'),
    ]
    _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12,4))

    for axIdx, port in enumerate(ports):
        for network, style in networks:
            csvFile = path.join(get_working_dir(), 'experiment_1', f'{network}_0.csv')
            ax = axs[axIdx]
            process_signals(ax, csvFile, port, style, 1., 0.75)

    axs[0].set_ylabel('$\widetilde{y}_{11}(t)$')
    axs[1].set_ylabel('$\widetilde{y}_{21}(t)$')

    axs[0].set_xlabel('$t$')
    axs[1].set_xlabel('$t$')

    for i in range(2):
        axs[i].yaxis.label.set_fontsize(16)
        axs[i].xaxis.label.set_fontsize(16)

    axs[1].legend(handles=[
        Line2D([0], [0], label='$Co-simulation$', color='k', alpha=1., linewidth=0.75),
        Line2D([0], [0], label='$Monolithic$', color='r', linestyle='dashed', alpha=0.9, linewidth=0.5),
    ])


    ports = [('twomass', 3), ('twomass', 16)]
    for axIdx, port in enumerate(ports):
        csvFile = path.join(get_working_dir(), 'experiment_1', 'cvode_monolithic_2.csv')
        process_signals(axs[axIdx], csvFile, port, 'r--', 0.9, 0.5)

    plt.savefig(path.join(get_working_dir(), 'fixedstepresponse.pdf'), bbox_inches='tight')


def experiment_2():
    print('Running Experiment 2.')
    produce_results(2)

    ports = [('tau2omega', 1), ('omega2tau', 1)]
    networks = [
        ('euler_subsystem', 'k'),   
    ]
    _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12,4))

    for axIdx, port in enumerate(ports):
        for network, style in networks:
            csvFile = path.join(get_working_dir(), 'experiment_2', f'{network}_0.csv')
            ax = axs[axIdx]
            process_signals(ax, csvFile, port, style, 1., 0.75)

    axs[0].set_ylabel('$\widetilde{y}_{11}(t)$')
    axs[1].set_ylabel('$\widetilde{y}_{21}(t)$')

    axs[0].set_xlabel('$t$')
    axs[1].set_xlabel('$t$')

    for i in range(2):
        axs[i].yaxis.label.set_fontsize(16)
        axs[i].xaxis.label.set_fontsize(16)

    axs[1].legend(handles=[
        Line2D([0], [0], label='$Co-simulation$', color='k', alpha=1., linewidth=0.75),
        Line2D([0], [0], label='$Monolithic$', color='r', linestyle='dashed', alpha=0.9, linewidth=0.5),
    ])


    ports = [('twomass', 3), ('twomass', 16)]
    for axIdx, port in enumerate(ports):
        csvFile = path.join(get_working_dir(), 'experiment_2', 'cvode_monolithic_2.csv')
        process_signals(axs[axIdx], csvFile, port, 'r--', 0.9, 0.5)

    plt.savefig(path.join(get_working_dir(), 'variablestepresponse.pdf'), bbox_inches='tight')


def process_errors(ax, csvFileName, style, alpha=1, linewidth=1):
    with open(csvFileName, 'r') as csvFile:
        csvReader = csv.reader(csvFile, dialect='excel')
        rowSteps = [float(step) for step in next(csvReader) if step] # moving reader inside the loop
        rowValues = [float(value) for value in next(csvReader) if value]
        ax.plot(rowSteps, rowValues, style, alpha=alpha, linewidth=linewidth)


def get_ylabels(errType, port):
    if errType == 'errors' and port == 1:
        return '$RMS(\Delta \widetilde{y}_{11})$', '$RMS(\Delta \widetilde{y}_{21})$'
    elif errType == 'defects' and port == 1:
        return '$RMS(\delta \widetilde{y}_{11})$', '$RMS(\delta \widetilde{y}_{21})$'
    elif errType == 'defects' and port == 0:
        return '$RMS(\delta \widetilde{u}_{21})$', '$RMS(\delta \widetilde{u}_{11})$'


def add_legend(axs, euler0Condition, showM=True, showN=True):
    euler0 = r", RMS = 0 " if euler0Condition else ""
    mLabel = 'm_i = ' if showM else ''
    nLabel = 'n_i = ' if showN else ''
    axs[0].legend(handles=[
        Line2D([0], [0], label='$' + mLabel + nLabel + '0$', color='k', linestyle='--', alpha=0.5, linewidth=1),
        Line2D([0], [0], label='$' + mLabel + nLabel + '1'  + euler0  + '$', color='k', linestyle='-.', alpha=0.75, linewidth=2),
        Line2D([0], [0], label='$' + mLabel + nLabel + '2'  + euler0  + '$', color='k', linestyle=':', alpha=1, linewidth=3),
    ])

    axs[1].legend(handles=[
        Line2D([0], [0], label='$' + mLabel + nLabel + '0$', color='k', linestyle='--', alpha=0.5, linewidth=1),
        Line2D([0], [0], label=f'$' + mLabel + nLabel + '1'  + euler0  + '$', color='k', linestyle='-.', alpha=0.75, linewidth=2),
        Line2D([0], [0], label=f'$' + mLabel + nLabel + '2' + euler0  + '$', color='k', linestyle=':', alpha=1, linewidth=3),
    ])

def experiment_3():
    print('Running Experiment 3.')
    produce_results(3)

    instances = ['tau2omega', 'omega2tau']
    networks = [
        ('euler_subsystem', 'k'),
    ]
    xStyles = ['--', '-.', ':']

    for errType, port, figure in [('errors', 1, 'fixedsteperrors'), ('defects', 1, 'fixedstepoutputdefect'), ('defects', 0, 'fixedstepconnectiondefect')]:
        if errType == 'errors':
            _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))
            axs[0].set_xlabel('$h$')
        else:
            _, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(5, 10))
        for extrapolationOrder in range(3):
            for axIdx, instance in enumerate(instances):
                axs[axIdx].set_yscale('log')
                axs[axIdx].set_xscale('log')
                for network, color in networks:
                    csvFile = path.join(get_working_dir(), 'experiment_3', f'{errType}_{network}_{instance}_{port}_{extrapolationOrder}_{extrapolationOrder}.csv')
                    print(csvFile)
                    if not (figure == 'fixedstepoutputdefect' and network == 'euler_subsystem' and extrapolationOrder > 0):
                        process_errors(axs[axIdx], csvFile, color + xStyles[extrapolationOrder], 0.5 + extrapolationOrder * 0.25, 1 + extrapolationOrder)
                    else:
                        print(f'{figure} {network} {extrapolationOrder} {axIdx}')

        label1, label2 = get_ylabels(errType, port)
        axs[0].set_ylabel(label1)
        axs[1].set_ylabel(label2)

        axs[0].set_xlim([1e-3, 1.])
        axs[1].set_xlim([1e-3, 1.])
        
        axs[1].set_xlabel('$h$')

        add_legend(axs, figure == "fixedstepoutputdefect")

        plt.savefig(path.join(get_working_dir(), figure + '.pdf'), bbox_inches='tight')


def experiment_4():
    print('Running Experiment 4.')
    produce_results(4)

    instances = ['tau2omega', 'omega2tau']
    networks = [
        ('euler_subsystem', 'k'),
    ]
    xStyles = ['--', '-.', ':']

    for errType, port, figure in [('errors', 1, 'fixedsteperrorsu0'), ('defects', 1, 'fixedstepoutputdefectu0'), ('defects', 0, 'fixedstepconnectiondefectu0')]:
        _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))

        print(figure)
        for extrapolationOrder in range(3):
            for axIdx, instance in enumerate(instances):
                axs[axIdx].set_yscale('log')
                axs[axIdx].set_xscale('log')
                for network, style in networks:
                    csvFile = path.join(get_working_dir(), 'experiment_4', f'{errType}_{network}_{instance}_{port}_0_{extrapolationOrder}.csv')
                    print(csvFile)
                    if not (figure == 'fixedstepoutputdefectu0' and network == 'euler_subsystem' and extrapolationOrder > 0):
                        process_errors(axs[axIdx], csvFile, style + xStyles[extrapolationOrder], 0.5 + extrapolationOrder * 0.25, 1 + extrapolationOrder)

        label1, label2 = get_ylabels(errType, port)
        axs[0].set_ylabel(label1)
        axs[1].set_ylabel(label2)

        axs[0].set_xlim([1e-3, 1.])
        axs[1].set_xlim([1e-3, 1.])

        axs[0].set_xlabel('$h$')
        axs[1].set_xlabel('$h$')

        axs[0].xaxis.label.set_fontsize(16)
        axs[1].xaxis.label.set_fontsize(16)
        axs[0].yaxis.label.set_fontsize(16)
        axs[1].yaxis.label.set_fontsize(16)

        add_legend(axs, figure == 'fixedstepoutputdefectu0', showM=False)

        plt.savefig(path.join(get_working_dir(), figure + '.pdf'), bbox_inches='tight')


def experiment_5():
    print('Running Experiment 5.')
    produce_results(5)

    instances = ['tau2omega', 'omega2tau']
    networks = [
        ('euler_subsystem', 'k'),
    ]
    xStyles = ['--', '-.', ':']

    for errType, port, figure in [('errors', 1, 'fixedsteperrorsy0'), ('defects', 1, 'fixedstepoutputdefecty0'), ('defects', 0, 'fixedstepconnectiondefecty0')]:
        _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))

        for extrapolationOrder in range(3):
            for axIdx, instance in enumerate(instances):
                axs[axIdx].set_yscale('log')
                axs[axIdx].set_xscale('log')
                for network, style in networks:
                    csvFile = path.join(get_working_dir(), 'experiment_5', f'{errType}_{network}_{instance}_{port}_{extrapolationOrder}_0.csv')
                    print(csvFile)
                    process_errors(axs[axIdx], csvFile, style + xStyles[extrapolationOrder], 0.5 + extrapolationOrder * 0.25, 1 + extrapolationOrder)

        label1, label2 = get_ylabels(errType, port)
        axs[0].set_ylabel(label1)
        axs[1].set_ylabel(label2)

        axs[0].set_xlim([1e-3, 1.])
        axs[1].set_xlim([1e-3, 1.])

        axs[0].set_xlabel('$h$')
        axs[1].set_xlabel('$h$')

        axs[0].xaxis.label.set_fontsize(16)
        axs[1].xaxis.label.set_fontsize(16)
        axs[0].yaxis.label.set_fontsize(16)
        axs[1].yaxis.label.set_fontsize(16)

        add_legend(axs, False, showN=False)

        plt.savefig(path.join(get_working_dir(), figure + '.pdf'), bbox_inches='tight')


def _experiment_6_errors():
    instances = ['tau2omega', 'omega2tau']
    networks = [
        ('euler_subsystem', 'k'),
    ]
    xStyles = ['--', '-.', ':']

    for errType, port, figure in [('errors', 1, 'variablesteperrors'), ('defects', 1, 'variablestepoutputdefect'), ('defects', 0, 'variablestepconnectiondefect')]:
        _, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 6))


        for extrapolationOrder in range(3):
            for axIdx, instance in enumerate(instances):
                axs[axIdx].set_yscale('log')
                axs[axIdx].set_xscale('log')
                for network, style in networks:
                    csvFile = path.join(get_working_dir(), 'experiment_6', f'{errType}_{network}_{instance}_{port}_{extrapolationOrder}_{extrapolationOrder}.csv')
                    print(csvFile)
                    if not (figure == 'variablestepoutputdefect' and network == 'euler_subsystem' and extrapolationOrder > 0):
                        process_errors(axs[axIdx], csvFile, style + xStyles[extrapolationOrder], 0.5 + extrapolationOrder * 0.25, 1 + extrapolationOrder)

        tols = np.logspace(-3, 0)
        for i in range(2):
            if errType != 'errors':
                axs[i].plot(tols, tols, 'r--', alpha=0.2)

        label1, label2 = get_ylabels(errType, port)
        axs[0].set_ylabel(label1)
        axs[1].set_ylabel(label2)

        axs[0].set_xlim([1e-3, 1.])
        axs[1].set_xlim([1e-3, 1.])

        axs[0].set_xlabel('$tol$')
        axs[1].set_xlabel('$tol$')

        add_legend(axs, figure == 'variablestepoutputdefect')

        plt.savefig(path.join(get_working_dir(), figure + '.pdf'), bbox_inches='tight')


def _experiment_6_steps():
    networks = [
        ('euler_subsystem', 0, '$m_1 = n_i = 0$', 'k--', 1, .5),
        ('euler_subsystem', 1, '$m_1 = n_i = 1$', 'k-.', 2, .75),
        ('euler_subsystem', 2, '$m_1 = n_i = 2$', 'k:', 3, 1.),
    ]
    plt.figure()
    plt.plot()

    for network, extrapolationOrder, networkLabel, style, width, alpha in networks:
        csvFileName = path.join(get_working_dir(), 'experiment_6', f'step_size_{network}_{extrapolationOrder}_{extrapolationOrder}.csv')
        print(csvFileName)
        with open(csvFileName, 'r') as csvFile:
            csvReader = csv.reader(csvFile, dialect='excel')
            rowTols = [float(tol) for tol in next(csvReader) if tol] # moving reader inside the loop
            rowStepMins = [float(value) for value in next(csvReader) if value]                
            rowStepAvgs = [float(value) for value in next(csvReader) if value]
            rowStepMaxs = [float(value) for value in next(csvReader) if value]
            xtrapLabel = ', $m_i = n_i = {0}$'.format(extrapolationOrder) if networkLabel == '$Euler$' else ''
            plt.loglog(
                rowTols, rowStepAvgs, style,
                linewidth=width, alpha=alpha,
                label=networkLabel + xtrapLabel
            )

    plt.xlabel('$tol$')
    plt.ylabel('$\overline{H}$')
    plt.xlim([1e-3, 1.])
    plt.legend()
    plt.savefig(path.join(get_working_dir(), 'variablestepsteps.pdf'), bbox_inches='tight')


def experiment_6():
    print('Running Experiment 6.')
    produce_results(6)
    _experiment_6_errors()
    _experiment_6_steps()


def article_experiments():
    return [
        experiment_1,
        experiment_2,
        experiment_3,
        experiment_4,
        experiment_5,
        experiment_6,
    ]


def produce_all_figures():
    for experiment in article_experiments():
        experiment()


if __name__ == '__main__':
    args = argv[1:]
    if not args:
        produce_all_figures()
    else:
        figIdx = int(args[0]) - 1
        article_experiments()[figIdx]()
