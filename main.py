import mne
import numpy as np


#legacy
def get_freq(info):
    raw = info
    power, freq = raw.compute_psd().get_data(return_freqs=True)
    return freq



# функция обработки edf > mne.eeg
def input_eeg(filename, type):
    ten_twenty_montage = mne.channels.make_standard_montage('standard_1020')

    if(type =="edf") :
        data = mne.io.read_raw_edf(filename)
        raw_data = data.get_data()
        channels = data.ch_names
        print(channels)
        return data

    elif(type =="bdf"):
        data = mne.io.read_raw_bdf(filename)
        return data


# функция mna.eeg -> ассиметрия альфа ритма
# добавить больше лобных электродов, оформить лучше код
def get_alpha_assimetry(raw):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)
    alpha_frontal = None
    for index, band in enumerate(bands):
        if (band_names[index] == 'Alpha'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)
            a_fp1 = data[:, raw.ch_names.index('Fp1')]
            a_fp2 = data[:, raw.ch_names.index('Fp2')]
            alpha_frontal = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
    return alpha_frontal

def get_beta_assimetry(raw, cols):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)
    beta_assimetry = 0
    count = 0
    hight = 0
    lower = 100
    for index, band in enumerate(bands):
        if (band_names[index] == 'Beta'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            if 'Fp1' in cols and 'Fp2' in cols:
                a_fp1 = data[:, raw.ch_names.index('Fp1')]
                a_fp2 = data[:, raw.ch_names.index('Fp2')]
                beta_fp12 = ((a_fp2 - a_fp1) / (a_fp2 + a_fp1))[0]
                beta_assimetry += beta_fp12
                count+=1
                if beta_fp12 > hight:
                    hight = beta_fp12
                if beta_fp12 < lower:
                    lower = beta_fp12

            if 'F7' in cols and 'F8' in cols:
                a_f7 = data[:, raw.ch_names.index('F7')]
                a_f8 = data[:, raw.ch_names.index('F8')]
                beta_f78 = ((a_f7 - a_f8) / (a_f7 + a_f8))[0]
                beta_assimetry += beta_f78
                count+=1
                if beta_f78 > hight:
                    hight = beta_f78
                if beta_f78 < lower:
                    lower = beta_f78

            if 'F3' in cols and 'F4' in cols:
                a_f3 = data[:, raw.ch_names.index('F3')]
                a_f4 = data[:, raw.ch_names.index('F4')]
                beta_f34 = ((a_f3 - a_f4) / (a_f3 + a_f4))[0]
                beta_assimetry += beta_f34
                count+=1
                if beta_f34 > hight:
                    hight = beta_f34
                if beta_f34 < lower:
                    lower = beta_f34

            if 'T3' in cols and 'T4' in cols:
                a_t3 = data[:, raw.ch_names.index('T3')]
                a_t4 = data[:, raw.ch_names.index('T4')]
                beta_t34 = ((a_t3 - a_t4) / (a_t3 + a_t4))[0]
                beta_assimetry += beta_t34
                count+=1
                if beta_t34 > hight:
                    hight = beta_t34
                if beta_t34 < lower:
                    lower = beta_t34

            if 'C3' in cols and 'C4' in cols:
                a_c3 = data[:, raw.ch_names.index('C3')]
                a_c4 = data[:, raw.ch_names.index('C4')]
                beta_c34 = ((a_c3 - a_c4) / (a_c3 + a_c4))[0]
                beta_assimetry += beta_c34
                count+=1
                if beta_c34 > hight:
                    hight = beta_c34
                if beta_c34 < lower:
                    lower = beta_c34

            if 'P3' in cols and 'P4' in cols:
                a_p3 = data[:, raw.ch_names.index('P3')]
                a_p4 = data[:, raw.ch_names.index('P4')]
                beta_p34 = ((a_p3 - a_p4) / (a_p3 + a_p4))[0]
                beta_assimetry += beta_p34
                count+=1
                if beta_p34 > hight:
                    hight = beta_p34
                if beta_p34 < lower:
                    lower = beta_p34

            if 'T5' in cols and 'T6' in cols:
                a_t5 = data[:, raw.ch_names.index('T5')]
                a_t6 = data[:, raw.ch_names.index('T6')]
                beta_t56 = ((a_t5 - a_t6) / (a_t5 + a_t6))[0]
                beta_assimetry += beta_t56
                count+=1
                if beta_t56 > hight:
                    hight = beta_t56
                if beta_t56 < lower:
                    lower = beta_t56

            if 'O1' in cols and 'O2' in cols:
                a_o1 = data[:, raw.ch_names.index('O1')]
                a_o2 = data[:, raw.ch_names.index('O2')]
                beta_012 = ((a_o1 - a_o2) / (a_o1 + a_o2))[0]
                beta_assimetry += beta_012
                count+=1
                if beta_012 > hight:
                    hight = beta_012
                if beta_012 < lower:
                    lower = beta_012

    return (beta_assimetry/ count), hight, lower
#функции оценки состояния человека

#функции вычисления альфа,бета
def get_alpha_beta_average(raw, col):
    freq = get_freq(raw)
    band_names = np.array(['Delta', 'Theta', 'Alpha', 'Beta', 'Gamma', 'Slow', 'Low_beta'])
    filter_list = [[1, 3], [4, 7], [8, 12], [13, 30], [30, 43], [4, 13], [13, 17]]
    bands = []
    for filt in filter_list:
        pt = np.argwhere((freq >= filt[0]) & (freq <= filt[1])).reshape(-1)
        bands.append(pt)

    one_alpha = []
    one_beta = []

    for index, band in enumerate(bands):
        if (band_names[index] == 'Alpha'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            for x in col:
                item_1 = data[:, raw.ch_names.index(x)]
                one_alpha.append(item_1)
        if (band_names[index] == 'Beta'):
            power, freq = raw.compute_psd().get_data(return_freqs=True)
            data = power[::, band].mean(axis=1).reshape(1, -1)

            for x in col:
                item_2 = data[:, raw.ch_names.index(x)]
                one_beta.append(item_2)


    averge_dif_alpha_subject = sum(one_alpha) / len(one_alpha)
    averge_dif_bets_subject = sum(one_beta) / len(one_beta)
    return averge_dif_alpha_subject, averge_dif_bets_subject

#функции выявления паттернов выгорания

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    filename = "real_data\\Daniil Shirikov_close.edf"

    raw = input_eeg(filename, "edf")
    channels = raw.ch_names
    print(raw)

    alpha_frontal = get_alpha_assimetry(raw)
    print(alpha_frontal)

    averge_dif_alpha_subject, averge_dif_bets_subject = get_alpha_beta_average(raw,channels)
    print(averge_dif_alpha_subject, averge_dif_bets_subject)

    beta_assimenry, hight, lower = get_beta_assimetry(raw,channels)
    print(beta_assimenry, hight, lower )

# See PyCharm help at https://www.jetbrains.com/help/pycharm/