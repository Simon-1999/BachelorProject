import matplotlib.pyplot as plt
import matplotlib as mpl

def default_style():

    plt.style.use('default')
    plt.rc('axes', unicode_minus=False)
    mpl.rcParams['font.family'] = "TW Cen MT"
    # plt.rc('font', family='serif')
    # mpl.rcParams['font.size'] = 14
    # mpl.rcParams['patch.linewidth'] = .8
    # plt.rcParams['xtick.major.size'] = 5.0
    # plt.rcParams['xtick.minor.size'] = 3.0
    plt.rcParams['xtick.top'] = True
    # plt.rcParams['ytick.major.size'] = 5.0
    # plt.rcParams['ytick.minor.size'] = 3.0
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')    
    # plt.rc('legend', fontsize='x-small')
    # plt.minorticks_on()


def science_style():

    plt.style.use('default')
    # plt.rcParams["axes.labelweight"] = "light"
    # plt.rcParams["font.weight"] = "light"
    # plt.rc('font', family='serif')
    mpl.rcParams['font.family'] = "Times New Roman"

    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.major.size'] = 5.0
    plt.rcParams['xtick.minor.size'] = 3.0
    plt.rcParams['xtick.direction'] = 'in'

    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.major.size'] = 5.0
    plt.rcParams['ytick.minor.size'] = 3.0
    plt.rcParams['ytick.direction'] = 'in'